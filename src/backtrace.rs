//! Wavefront backtrace: reconstructs alignment CIGAR from stored wavefronts.
//!
//! The "high memory" backtrace walks backwards through all M-wavefronts,
//! using a piggyback encoding trick (offset + type packed into i64) to
//! determine which operation (mismatch, insertion, deletion) was taken
//! at each step.

use crate::bialign::ComponentType;
use crate::cigar::Cigar;
use crate::components::WavefrontComponents;
use crate::offset::{OFFSET_NULL, wavefront_h, wavefront_v};
use crate::penalties::{DistanceMetric, WavefrontPenalties};
use crate::slab::WavefrontSlab;

// Backtrace type constants (packed into low 4 bits of i64).
// Higher values win ties when offsets are equal.
const BACKTRACE_TYPE_BITS: u32 = 4;
const BACKTRACE_TYPE_MASK: i64 = 0xF;
const BT_I1_OPEN: i64 = 1;
const BT_I1_EXT: i64 = 2;
const BT_I2_OPEN: i64 = 3;
const BT_I2_EXT: i64 = 4;
const BT_D1_OPEN: i64 = 5;
const BT_D1_EXT: i64 = 6;
const BT_D2_OPEN: i64 = 7;
const BT_D2_EXT: i64 = 8;
const BT_M: i64 = 9;

#[inline(always)]
fn piggyback_set(offset: i32, bt_type: i64) -> i64 {
    ((offset as i64) << BACKTRACE_TYPE_BITS) | bt_type
}

#[inline(always)]
fn piggyback_get_offset(value: i64) -> i32 {
    (value >> BACKTRACE_TYPE_BITS) as i32
}

#[inline(always)]
fn piggyback_get_type(value: i64) -> i64 {
    value & BACKTRACE_TYPE_MASK
}

/// Look up mismatch predecessor: M-wavefront at (score, k), offset + 1.
#[inline(always)]
fn bt_misms<const N: usize>(wf_components: &WavefrontComponents<N>, score: i32, k: i32) -> i64 {
    if score < 0 {
        return OFFSET_NULL as i64;
    }
    let m_ptr = wf_components.get_m_ptr(score as usize);
    if m_ptr.is_null() {
        return OFFSET_NULL as i64;
    }
    let wf = unsafe { &*m_ptr };
    if wf.wf_elements_init_min <= k && k <= wf.wf_elements_init_max {
        piggyback_set(wf.get_offset(k) + 1, BT_M)
    } else {
        OFFSET_NULL as i64
    }
}

/// Look up insertion predecessor: M-wavefront at (score, k-1), offset + 1.
#[inline(always)]
fn bt_ins1_open<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> i64 {
    if score < 0 {
        return OFFSET_NULL as i64;
    }
    let m_ptr = wf_components.get_m_ptr(score as usize);
    if m_ptr.is_null() {
        return OFFSET_NULL as i64;
    }
    let wf = unsafe { &*m_ptr };
    if wf.wf_elements_init_min < k && k - 1 <= wf.wf_elements_init_max {
        piggyback_set(wf.get_offset(k - 1) + 1, BT_I1_OPEN)
    } else {
        OFFSET_NULL as i64
    }
}

/// Look up deletion predecessor: M-wavefront at (score, k+1), offset unchanged.
#[inline(always)]
fn bt_del1_open<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> i64 {
    if score < 0 {
        return OFFSET_NULL as i64;
    }
    let m_ptr = wf_components.get_m_ptr(score as usize);
    if m_ptr.is_null() {
        return OFFSET_NULL as i64;
    }
    let wf = unsafe { &*m_ptr };
    if wf.wf_elements_init_min <= k + 1 && k < wf.wf_elements_init_max {
        piggyback_set(wf.get_offset(k + 1), BT_D1_OPEN)
    } else {
        OFFSET_NULL as i64
    }
}

/// Write `num_matches` 'M' operations into the CIGAR (backwards).
#[inline(always)]
fn write_matches(cigar: &mut Cigar, num_matches: i32) {
    if num_matches < 0 {
        // Only reachable from a corrupted backtrace (X-drop with
        // steps_between_cutoffs > 1, a C bug kept for parity). C's
        // `wavefront_backtrace_matches` still does `begin_offset -= num_matches`,
        // discarding the last |num_matches| operations written; do the same.
        cigar.begin_offset += num_matches.unsigned_abs() as usize;
        return;
    }
    if num_matches == 0 {
        return;
    }
    let n = num_matches as usize;
    let start = cigar.begin_offset + 1 - n;
    cigar.operations[start..cigar.begin_offset + 1].fill(b'M');
    cigar.begin_offset = start - 1;
}

/// Perform backtrace for linear (edit/indel) distance, building a CIGAR.
///
/// Walks backwards through all stored M-wavefronts from the alignment endpoint
/// to reconstruct the full alignment. The CIGAR is built from end to beginning.
#[allow(clippy::too_many_arguments)]
pub fn backtrace_linear<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    _slab: &WavefrontSlab,
    penalties: &WavefrontPenalties,
    pattern_length: i32,
    text_length: i32,
    alignment_score: i32,
    alignment_k: i32,
    alignment_offset: i32,
    cigar: &mut Cigar,
) {
    let distance_metric = penalties.distance_metric;

    // Prepare cigar — write backwards from the end of the buffer
    let max_ops = cigar.operations.len();
    cigar.end_offset = max_ops;
    cigar.begin_offset = max_ops - 1;

    let mut score = alignment_score;
    let mut k = alignment_k;
    let mut offset = alignment_offset;
    let mut v = wavefront_v(k, offset);
    let mut h = wavefront_h(k, offset);

    // Account for ending insertions/deletions (ends-free)
    if v < pattern_length {
        for _ in 0..(pattern_length - v) {
            cigar.operations[cigar.begin_offset] = b'D';
            cigar.begin_offset -= 1;
        }
    }
    if h < text_length {
        for _ in 0..(text_length - h) {
            cigar.operations[cigar.begin_offset] = b'I';
            cigar.begin_offset -= 1;
        }
    }

    // Main backtrace loop
    while v > 0 && h > 0 && score > 0 {
        // Compute predecessor scores
        let mismatch_score = score - penalties.mismatch;
        let gap_score = score - penalties.gap_opening1;

        // Look up predecessor offsets (packed with type)
        let misms = if distance_metric != DistanceMetric::Indel {
            bt_misms(wf_components, mismatch_score, k)
        } else {
            OFFSET_NULL as i64
        };
        let ins = bt_ins1_open(wf_components, gap_score, k);
        let del = bt_del1_open(wf_components, gap_score, k);

        // Take max — higher type wins ties (M > D > I)
        let max_all = misms.max(ins.max(del));

        if max_all < 0 {
            break; // No valid predecessor
        }

        // Emit matches (gap between current offset and predecessor offset)
        let max_offset = piggyback_get_offset(max_all);
        let num_matches = offset - max_offset;
        write_matches(cigar, num_matches);
        offset = max_offset;

        // Check boundary
        v = wavefront_v(k, offset);
        h = wavefront_h(k, offset);
        if v <= 0 || h <= 0 {
            break;
        }

        // Emit the operation
        let bt_type = piggyback_get_type(max_all);
        match bt_type {
            BT_M => {
                score = mismatch_score;
                cigar.operations[cigar.begin_offset] = b'X';
                cigar.begin_offset -= 1;
                offset -= 1;
            }
            BT_I1_OPEN => {
                score = gap_score;
                cigar.operations[cigar.begin_offset] = b'I';
                cigar.begin_offset -= 1;
                k -= 1;
                offset -= 1;
            }
            BT_D1_OPEN => {
                score = gap_score;
                cigar.operations[cigar.begin_offset] = b'D';
                cigar.begin_offset -= 1;
                k += 1;
            }
            _ => panic!("Invalid backtrace type: {}", bt_type),
        }

        // Update coordinates
        v = wavefront_v(k, offset);
        h = wavefront_h(k, offset);
    }

    // Account for beginning matches
    if v > 0 && h > 0 {
        let num_matches = v.min(h);
        write_matches(cigar, num_matches);
        v -= num_matches;
        h -= num_matches;
    }

    // Account for beginning insertions/deletions
    while v > 0 {
        cigar.operations[cigar.begin_offset] = b'D';
        cigar.begin_offset -= 1;
        v -= 1;
    }
    while h > 0 {
        cigar.operations[cigar.begin_offset] = b'I';
        cigar.begin_offset -= 1;
        h -= 1;
    }

    // Adjust to point to first valid operation
    cigar.begin_offset += 1;
    cigar.score = alignment_score;
}

// --- Affine backtrace helpers ---

/// Look up both M→I1 open (k-1, offset+1) and M→D1 open (k+1, offset) from
/// the same M-wavefront at the given score. Returns (ins1_open, del1_open).
/// This avoids a redundant double-fetch.
#[inline(always)]
fn bt_open1_pair<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> (i64, i64) {
    if score < 0 {
        return (OFFSET_NULL as i64, OFFSET_NULL as i64);
    }
    let m_ptr = wf_components.get_m_ptr(score as usize);
    if m_ptr.is_null() {
        return (OFFSET_NULL as i64, OFFSET_NULL as i64);
    }
    let wf = unsafe { &*m_ptr };
    let min = wf.wf_elements_init_min;
    let max = wf.wf_elements_init_max;

    let ins = if min < k && k - 1 <= max {
        piggyback_set(wf.get_offset(k - 1) + 1, BT_I1_OPEN)
    } else {
        OFFSET_NULL as i64
    };
    let del = if min <= k + 1 && k < max {
        piggyback_set(wf.get_offset(k + 1), BT_D1_OPEN)
    } else {
        OFFSET_NULL as i64
    };
    (ins, del)
}

/// Affine backtrace: look up I1 extend at (score, k-1), offset + 1.
#[inline(always)]
fn bt_affine_ins1_ext<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> i64 {
    if score < 0 {
        return OFFSET_NULL as i64;
    }
    let i1_ptr = wf_components.get_i1_ptr(score as usize);
    if i1_ptr.is_null() {
        return OFFSET_NULL as i64;
    }
    let wf = unsafe { &*i1_ptr };
    if wf.wf_elements_init_min < k && k - 1 <= wf.wf_elements_init_max {
        piggyback_set(wf.get_offset(k - 1) + 1, BT_I1_EXT)
    } else {
        OFFSET_NULL as i64
    }
}

/// Affine backtrace: look up D1 extend at (score, k+1), offset unchanged.
#[inline(always)]
fn bt_affine_del1_ext<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> i64 {
    if score < 0 {
        return OFFSET_NULL as i64;
    }
    let d1_ptr = wf_components.get_d1_ptr(score as usize);
    if d1_ptr.is_null() {
        return OFFSET_NULL as i64;
    }
    let wf = unsafe { &*d1_ptr };
    if wf.wf_elements_init_min <= k + 1 && k < wf.wf_elements_init_max {
        piggyback_set(wf.get_offset(k + 1), BT_D1_EXT)
    } else {
        OFFSET_NULL as i64
    }
}

/// Perform backtrace for gap-affine distance, building a CIGAR.
///
/// Uses a state machine tracking current matrix type (M, I1, D1).
/// In M-state: predecessors are mismatch (M→M), gap open (M→I1/D1), or I1/D1→M.
/// In I1-state: predecessors are gap open (M→I1) or gap extend (I1→I1).
/// In D1-state: predecessors are gap open (M→D1) or gap extend (D1→D1).
#[allow(clippy::too_many_arguments)]
pub fn backtrace_affine<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    _slab: &WavefrontSlab,
    penalties: &WavefrontPenalties,
    pattern_length: i32,
    text_length: i32,
    alignment_score: i32,
    alignment_k: i32,
    alignment_offset: i32,
    component_begin: ComponentType,
    component_end: ComponentType,
    cigar: &mut Cigar,
) {
    let _ = component_begin; // C's backtrace does not use it
    backtrace_affine_c(
        wf_components,
        penalties,
        pattern_length,
        text_length,
        alignment_score,
        alignment_k,
        alignment_offset,
        component_end,
        cigar,
    );
}

// --- Affine2p backtrace helpers ---

/// Look up both M→I2 open (k-1, offset+1) and M→D2 open (k+1, offset) from
/// the same M-wavefront at the given score. Returns (ins2_open, del2_open).
#[inline(always)]
fn bt_open2_pair<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> (i64, i64) {
    if score < 0 {
        return (OFFSET_NULL as i64, OFFSET_NULL as i64);
    }
    let m_ptr = wf_components.get_m_ptr(score as usize);
    if m_ptr.is_null() {
        return (OFFSET_NULL as i64, OFFSET_NULL as i64);
    }
    let wf = unsafe { &*m_ptr };
    let min = wf.wf_elements_init_min;
    let max = wf.wf_elements_init_max;

    let ins = if min < k && k - 1 <= max {
        piggyback_set(wf.get_offset(k - 1) + 1, BT_I2_OPEN)
    } else {
        OFFSET_NULL as i64
    };
    let del = if min <= k + 1 && k < max {
        piggyback_set(wf.get_offset(k + 1), BT_D2_OPEN)
    } else {
        OFFSET_NULL as i64
    };
    (ins, del)
}

/// Affine2p backtrace: look up I2 extend at (score, k-1), offset + 1.
#[inline(always)]
fn bt_affine_ins2_ext<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> i64 {
    if score < 0 {
        return OFFSET_NULL as i64;
    }
    let i2_ptr = wf_components.get_i2_ptr(score as usize);
    if i2_ptr.is_null() {
        return OFFSET_NULL as i64;
    }
    let wf = unsafe { &*i2_ptr };
    if wf.wf_elements_init_min < k && k - 1 <= wf.wf_elements_init_max {
        piggyback_set(wf.get_offset(k - 1) + 1, BT_I2_EXT)
    } else {
        OFFSET_NULL as i64
    }
}

/// Affine2p backtrace: look up D2 extend at (score, k+1), offset unchanged.
#[inline(always)]
fn bt_affine_del2_ext<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    score: i32,
    k: i32,
) -> i64 {
    if score < 0 {
        return OFFSET_NULL as i64;
    }
    let d2_ptr = wf_components.get_d2_ptr(score as usize);
    if d2_ptr.is_null() {
        return OFFSET_NULL as i64;
    }
    let wf = unsafe { &*d2_ptr };
    if wf.wf_elements_init_min <= k + 1 && k < wf.wf_elements_init_max {
        piggyback_set(wf.get_offset(k + 1), BT_D2_EXT)
    } else {
        OFFSET_NULL as i64
    }
}

/// Perform backtrace for gap-affine 2-piece distance, building a CIGAR.
///
/// Uses a 5-state machine tracking current matrix type (M, I1, D1, I2, D2).
#[allow(clippy::too_many_arguments)]
pub fn backtrace_affine2p<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    _slab: &WavefrontSlab,
    penalties: &WavefrontPenalties,
    pattern_length: i32,
    text_length: i32,
    alignment_score: i32,
    alignment_k: i32,
    alignment_offset: i32,
    component_begin: ComponentType,
    component_end: ComponentType,
    cigar: &mut Cigar,
) {
    let _ = component_begin; // C's backtrace does not use it
    backtrace_affine_c(
        wf_components,
        penalties,
        pattern_length,
        text_length,
        alignment_score,
        alignment_k,
        alignment_offset,
        component_end,
        cigar,
    );
}

/// Backtrace for gap-affine and gap-affine 2-piece, ported from C's
/// `wavefront_backtrace_affine` so tie-breaking and boundary behaviour match.
///
/// In the M state every predecessor is read from its *source* wavefront:
/// mismatch (M at s-x), gap open (M at s-o-e, k∓1) and gap extend (I/D at
/// s-e, k∓1), taking the maximum packed (offset, type). Reading the stored
/// I/D cell at score s instead can miss sources when heuristics narrowed the
/// computed range. Leading matches/indels are added only if the trace ends
/// in M (as C).
#[allow(clippy::too_many_arguments)]
fn backtrace_affine_c<const N: usize>(
    wf_components: &WavefrontComponents<N>,
    penalties: &WavefrontPenalties,
    pattern_length: i32,
    text_length: i32,
    alignment_score: i32,
    alignment_k: i32,
    alignment_offset: i32,
    component_end: ComponentType,
    cigar: &mut Cigar,
) {
    let two_piece = penalties.distance_metric == DistanceMetric::GapAffine2p;
    let push = |cigar: &mut Cigar, op: u8| {
        cigar.operations[cigar.begin_offset] = op;
        cigar.begin_offset -= 1;
    };

    // Prepare cigar — write backwards from the end of the buffer
    let max_ops = cigar.operations.len();
    cigar.end_offset = max_ops;
    cigar.begin_offset = max_ops - 1;

    let mut matrix = component_end;
    let mut score = alignment_score;
    let mut k = alignment_k;
    let mut offset = alignment_offset;
    let mut v = wavefront_v(k, offset);
    let mut h = wavefront_h(k, offset);

    // Account for ending insertions/deletions (ends-free)
    if component_end == ComponentType::M {
        for _ in 0..(pattern_length - v).max(0) {
            push(cigar, b'D');
        }
        for _ in 0..(text_length - h).max(0) {
            push(cigar, b'I');
        }
    }

    while v > 0 && h > 0 && score > 0 {
        let mismatch = score - penalties.mismatch;
        let gap_open1 = score - penalties.gap_opening1 - penalties.gap_extension1;
        let gap_open2 = score - penalties.gap_opening2 - penalties.gap_extension2;
        let gap_extend1 = score - penalties.gap_extension1;
        let gap_extend2 = score - penalties.gap_extension2;
        let max_all = match matrix {
            ComponentType::M => {
                let misms = bt_misms(wf_components, mismatch, k);
                let (ins1_open, del1_open) = bt_open1_pair(wf_components, gap_open1, k);
                let max_ins1 = ins1_open.max(bt_affine_ins1_ext(wf_components, gap_extend1, k));
                let max_del1 = del1_open.max(bt_affine_del1_ext(wf_components, gap_extend1, k));
                if two_piece {
                    let (ins2_open, del2_open) = bt_open2_pair(wf_components, gap_open2, k);
                    let max_ins2 = ins2_open.max(bt_affine_ins2_ext(wf_components, gap_extend2, k));
                    let max_del2 = del2_open.max(bt_affine_del2_ext(wf_components, gap_extend2, k));
                    misms.max(max_ins1.max(max_ins2).max(max_del1.max(max_del2)))
                } else {
                    misms.max(max_ins1.max(max_del1))
                }
            }
            ComponentType::I1 => bt_open1_pair(wf_components, gap_open1, k)
                .0
                .max(bt_affine_ins1_ext(wf_components, gap_extend1, k)),
            ComponentType::D1 => bt_open1_pair(wf_components, gap_open1, k)
                .1
                .max(bt_affine_del1_ext(wf_components, gap_extend1, k)),
            ComponentType::I2 => bt_open2_pair(wf_components, gap_open2, k)
                .0
                .max(bt_affine_ins2_ext(wf_components, gap_extend2, k)),
            ComponentType::D2 => bt_open2_pair(wf_components, gap_open2, k)
                .1
                .max(bt_affine_del2_ext(wf_components, gap_extend2, k)),
        };
        if max_all < 0 {
            break; // No source
        }
        // Traceback matches
        if matrix == ComponentType::M {
            let max_offset = piggyback_get_offset(max_all);
            write_matches(cigar, offset - max_offset);
            offset = max_offset;
            v = wavefront_v(k, offset);
            h = wavefront_h(k, offset);
            if v <= 0 || h <= 0 {
                break;
            }
        }
        // Traceback operation
        let bt_type = piggyback_get_type(max_all);
        (score, matrix) = match bt_type {
            BT_M => (mismatch, ComponentType::M),
            BT_I1_OPEN | BT_D1_OPEN => (gap_open1, ComponentType::M),
            BT_I2_OPEN | BT_D2_OPEN => (gap_open2, ComponentType::M),
            BT_I1_EXT => (gap_extend1, ComponentType::I1),
            BT_I2_EXT => (gap_extend2, ComponentType::I2),
            BT_D1_EXT => (gap_extend1, ComponentType::D1),
            BT_D2_EXT => (gap_extend2, ComponentType::D2),
            _ => panic!("Invalid backtrace type: {bt_type}"),
        };
        match bt_type {
            BT_M => {
                push(cigar, b'X');
                offset -= 1;
            }
            BT_I1_OPEN | BT_I1_EXT | BT_I2_OPEN | BT_I2_EXT => {
                push(cigar, b'I');
                k -= 1;
                offset -= 1;
            }
            _ => {
                push(cigar, b'D');
                k += 1;
            }
        }
        v = wavefront_v(k, offset);
        h = wavefront_h(k, offset);
    }

    // Account for last operations (only when the trace ends in M)
    if matrix == ComponentType::M {
        if v > 0 && h > 0 {
            let num_matches = v.min(h);
            write_matches(cigar, num_matches);
            v -= num_matches;
            h -= num_matches;
        }
        for _ in 0..v {
            push(cigar, b'D');
        }
        for _ in 0..h {
            push(cigar, b'I');
        }
    } else {
        // C aborts here ("I?/D?-Beginning backtrace error").
        debug_assert!(
            v == 0 && h == 0 && (score == 0 || penalties.match_ != 0),
            "I?/D?-beginning backtrace error (v={v}, h={h}, score={score})"
        );
    }

    // Adjust to point to first valid operation
    cigar.begin_offset += 1;
    cigar.score = alignment_score;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aligner::{AlignStatus, AlignmentScope, EditAligner};

    #[test]
    fn test_backtrace_identical() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
        assert_eq!(aligner.cigar().to_string_rle(true), "4M");
        aligner.cigar().check_alignment(b"ACGT", b"ACGT").unwrap();
    }

    #[test]
    fn test_backtrace_one_mismatch() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 1);
        assert_eq!(aligner.cigar().to_string_rle(true), "2M1X1M");
        aligner.cigar().check_alignment(b"ACGT", b"ACTT").unwrap();
    }

    #[test]
    fn test_backtrace_one_insertion() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACT", b"ACGT");
        assert_eq!(score, 1);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"ACT", b"ACGT").unwrap();
        assert_eq!(cigar.score_edit(), 1);
    }

    #[test]
    fn test_backtrace_one_deletion() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACGT", b"ACT");
        assert_eq!(score, 1);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"ACGT", b"ACT").unwrap();
        assert_eq!(cigar.score_edit(), 1);
    }

    #[test]
    fn test_backtrace_multiple_edits() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"kitten", b"sitting");
        assert_eq!(score, 3);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"kitten", b"sitting").unwrap();
        assert_eq!(cigar.score_edit(), 3);
    }

    #[test]
    fn test_backtrace_empty_vs_nonempty() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"", b"ACGT");
        assert_eq!(score, 4);
        assert_eq!(aligner.cigar().to_string_rle(true), "4I");
    }

    #[test]
    fn test_backtrace_nonempty_vs_empty() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACGT", b"");
        assert_eq!(score, 4);
        assert_eq!(aligner.cigar().to_string_rle(true), "4D");
    }

    #[test]
    fn test_backtrace_indel_distance() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_indel());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 2);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"ACGT", b"ACTT").unwrap();
        // Indel distance doesn't allow mismatches, so no X in CIGAR
        let rle = cigar.to_string_rle(true);
        assert!(
            !rle.contains('X'),
            "Indel CIGAR should not contain X: {}",
            rle
        );
    }

    #[test]
    fn test_backtrace_completely_different() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"AAAA", b"TTTT");
        assert_eq!(score, 4);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"AAAA", b"TTTT").unwrap();
        assert_eq!(cigar.score_edit(), 4);
    }
}
