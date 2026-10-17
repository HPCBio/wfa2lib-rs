//! Heuristic strategies for wavefront pruning.
//!
//! Heuristics reduce computation by pruning diagonals from the wavefront
//! that are unlikely to contribute to the optimal alignment. This trades
//! exactness for speed on divergent sequences.

#![allow(clippy::too_many_arguments)]

use crate::offset::{DPMATRIX_DIAGONAL_NULL, OFFSET_NULL, WfOffset, wavefront_h, wavefront_v};
use crate::wavefront::Wavefront;

/// Heuristic strategy configuration.
#[derive(Debug, Clone)]
pub enum HeuristicStrategy {
    /// No heuristic pruning.
    None,
    /// WF-Adaptive: prune diagonals far behind the wavefront frontier.
    WfAdaptive {
        min_wavefront_length: i32,
        max_distance_threshold: i32,
        steps_between_cutoffs: i32,
    },
    /// X-Drop: prune diagonals whose Smith-Waterman score has dropped `xdrop`
    /// or more below the best score seen, every `steps_between_cutoffs` steps.
    ///
    /// **Known issue (parity with C WFA2-lib):** with `steps_between_cutoffs > 1`
    /// the CIGAR can be invalid (C's backtrace reads wavefronts the cutoff
    /// trimmed to empty, HPCBio/wfa2lib-rs#21). Scores are unaffected, and
    /// `steps_between_cutoffs == 1` produces valid CIGARs. Fixed by
    /// [`crate::aligner::Compat::Corrected`].
    XDrop {
        xdrop: i32,
        steps_between_cutoffs: i32,
    },
    /// Z-Drop: stop the alignment when the wavefront's best Smith-Waterman
    /// score falls more than `zdrop` below the best seen, checked every
    /// `steps_between_cutoffs` steps. The result is the partial alignment up to
    /// the best-scoring cell.
    ZDrop {
        zdrop: i32,
        steps_between_cutoffs: i32,
    },
    /// Static banding: clip wavefront to fixed diagonal range.
    BandedStatic { min_k: i32, max_k: i32 },
    /// Adaptive banding: move band dynamically based on diagonal progress.
    BandedAdaptive {
        min_k: i32,
        max_k: i32,
        steps_between_cutoffs: i32,
    },
}

/// Mutable heuristic state tracked across alignment steps.
#[derive(Debug, Clone)]
pub struct HeuristicState {
    pub steps_wait: i32,
    pub max_sw_score: i32,
    pub max_wf_score: i32,
    pub max_sw_score_k: i32,
    pub max_sw_score_offset: WfOffset,
}

impl Default for HeuristicState {
    fn default() -> Self {
        Self {
            steps_wait: 0,
            max_sw_score: 0,
            max_wf_score: 0,
            max_sw_score_k: DPMATRIX_DIAGONAL_NULL,
            max_sw_score_offset: OFFSET_NULL,
        }
    }
}

impl HeuristicState {
    pub fn reset(&mut self, strategy: &HeuristicStrategy) {
        *self = Self::default();
        self.steps_wait = match strategy {
            HeuristicStrategy::WfAdaptive {
                steps_between_cutoffs,
                ..
            } => *steps_between_cutoffs,
            HeuristicStrategy::BandedAdaptive {
                steps_between_cutoffs,
                ..
            }
            | HeuristicStrategy::XDrop {
                steps_between_cutoffs,
                ..
            }
            | HeuristicStrategy::ZDrop {
                steps_between_cutoffs,
                ..
            } => *steps_between_cutoffs,
            _ => 0,
        };
    }
}

/// Compute end-to-end distance for a diagonal: how far from completing alignment.
///
/// Returns `max(pattern_remaining, text_remaining)` or `i32::MAX` for null offsets.
#[inline(always)]
pub fn wf_distance_end2end(offset: WfOffset, k: i32, plen: i32, tlen: i32) -> i32 {
    if offset < 0 {
        return i32::MAX;
    }
    let v = wavefront_v(k, offset);
    let h = wavefront_h(k, offset);
    let left_v = plen - v;
    let left_h = tlen - h;
    left_v.max(left_h)
}

/// Apply heuristic cutoff to the M-wavefront.
///
/// Returns `true` if alignment should be terminated (z-drop only); the
/// best-scoring cell is then in `state.max_sw_score_{k,offset}` (reached at
/// score `state.max_wf_score`). Modifies `wf`'s lo/hi bounds in-place.
/// `wf_match` is the wavefront match penalty (`<= 0`).
pub fn heuristic_cutoff(
    strategy: &HeuristicStrategy,
    state: &mut HeuristicState,
    wf: &mut Wavefront,
    distances: &mut [WfOffset],
    distances_base_k: i32,
    pattern_length: i32,
    text_length: i32,
    wf_match: i32,
    score: i32,
) -> bool {
    if wf.lo > wf.hi {
        return false;
    }

    // Decrease wait steps
    state.steps_wait -= 1;

    match strategy {
        HeuristicStrategy::None => false,
        HeuristicStrategy::WfAdaptive {
            min_wavefront_length,
            max_distance_threshold,
            steps_between_cutoffs,
        } => {
            wfadaptive_cutoff(
                state,
                wf,
                distances,
                distances_base_k,
                pattern_length,
                text_length,
                *min_wavefront_length,
                *max_distance_threshold,
                *steps_between_cutoffs,
            );
            false
        }
        HeuristicStrategy::BandedStatic { min_k, max_k } => {
            banded_static_cutoff(wf, *min_k, *max_k);
            false
        }
        HeuristicStrategy::XDrop {
            xdrop,
            steps_between_cutoffs,
        } => {
            xdrop_cutoff(
                state,
                wf,
                distances,
                distances_base_k,
                wf_match,
                score,
                *xdrop,
                *steps_between_cutoffs,
            );
            false
        }
        HeuristicStrategy::ZDrop {
            zdrop,
            steps_between_cutoffs,
        } => zdrop_cutoff(
            state,
            wf,
            distances,
            distances_base_k,
            wf_match,
            score,
            *zdrop,
            *steps_between_cutoffs,
        ),
        HeuristicStrategy::BandedAdaptive {
            min_k,
            max_k,
            steps_between_cutoffs,
        } => {
            banded_adaptive_cutoff(
                state,
                wf,
                pattern_length,
                text_length,
                *min_k,
                *max_k,
                *steps_between_cutoffs,
            );
            false
        }
    }
}

/// WF-Adaptive heuristic: prune diagonals far from the wavefront frontier.
fn wfadaptive_cutoff(
    state: &mut HeuristicState,
    wf: &mut Wavefront,
    distances: &mut [WfOffset],
    distances_base_k: i32,
    pattern_length: i32,
    text_length: i32,
    min_wavefront_length: i32,
    max_distance_threshold: i32,
    steps_between_cutoffs: i32,
) {
    // Check steps
    if state.steps_wait > 0 {
        return;
    }

    // Check minimum wavefront length
    let base_lo = wf.lo;
    let base_hi = wf.hi;
    if (base_hi - base_lo + 1) < min_wavefront_length {
        return;
    }

    // Compute distance for each diagonal
    let offsets = wf.offsets_slice();
    let base_k = wf.base_k();
    let mut min_distance = pattern_length.max(text_length);

    for k in base_lo..=base_hi {
        let offset = offsets[(k - base_k) as usize];
        let distance = wf_distance_end2end(offset, k, pattern_length, text_length);
        distances[(k - distances_base_k) as usize] = distance;
        if distance < min_distance {
            min_distance = distance;
        }
    }

    // Reduce: preserve target diagonal
    let alignment_k = text_length - pattern_length;
    wfadaptive_reduce(
        wf,
        distances,
        distances_base_k,
        min_distance,
        max_distance_threshold,
        alignment_k,
        alignment_k,
    );

    // Reset wait steps
    state.steps_wait = steps_between_cutoffs;
}

/// Reduce wavefront bounds by removing diagonals far from minimum distance.
fn wfadaptive_reduce(
    wf: &mut Wavefront,
    distances: &[WfOffset],
    distances_base_k: i32,
    min_distance: i32,
    max_distance_threshold: i32,
    min_k: i32,
    max_k: i32,
) {
    // Reduce from bottom (lo)
    let top_limit = max_k.min(wf.hi);
    let mut lo_reduced = wf.lo;
    for k in wf.lo..top_limit {
        let d = distances[(k - distances_base_k) as usize];
        if d - min_distance <= max_distance_threshold {
            break;
        }
        lo_reduced += 1;
    }
    wf.lo = lo_reduced;

    // Reduce from top (hi)
    let bottom_limit = min_k.max(wf.lo);
    let mut hi_reduced = wf.hi;
    let mut k = wf.hi;
    while k > bottom_limit {
        let d = distances[(k - distances_base_k) as usize];
        if d - min_distance <= max_distance_threshold {
            break;
        }
        hi_reduced -= 1;
        k -= 1;
    }
    wf.hi = hi_reduced;
}

/// Static banding: clip wavefront to fixed diagonal range.
/// Adaptive band (C's `wavefront_heuristic_banded_adaptive`): every
/// `steps_between_cutoffs` steps, if the wavefront is wider than the band
/// (`max_k - min_k + 1`), keep a band-wide window, shifted towards the end of
/// the wavefront whose sampled diagonals are closer to the end of the
/// alignment. Wavefronts shorter than 4 are left alone (without resetting the
/// wait, as in C).
fn banded_adaptive_cutoff(
    state: &mut HeuristicState,
    wf: &mut Wavefront,
    pattern_length: i32,
    text_length: i32,
    min_k: i32,
    max_k: i32,
    steps_between_cutoffs: i32,
) {
    if state.steps_wait > 0 {
        return;
    }
    let (lo, hi) = (wf.lo, wf.hi);
    let wf_length = hi - lo + 1;
    if wf_length < 4 {
        return;
    }
    let max_wf_length = max_k - min_k + 1;
    if wf_length > max_wf_length {
        let leeway = (wf_length - max_wf_length) / 2;
        let quarter = wf_length / 4;
        let dist = |k: i32| wf_distance_end2end(wf.get_offset(k), k, pattern_length, text_length);
        let (d0, d1, d2, d3) = (
            dist(lo),
            dist(lo + quarter),
            dist(lo + 2 * quarter),
            dist(hi),
        );
        let mut new_lo = lo;
        if d0 > d3 {
            new_lo += leeway;
        }
        if d1 > d2 {
            new_lo += leeway;
        }
        wf.lo = new_lo.max(lo);
        wf.hi = (new_lo + max_wf_length - 1).min(hi);
    }
    state.steps_wait = steps_between_cutoffs;
}

fn banded_static_cutoff(wf: &mut Wavefront, min_k: i32, max_k: i32) {
    if wf.lo < min_k {
        wf.lo = min_k;
    }
    if wf.hi > max_k {
        wf.hi = max_k;
    }
}

/// X-Drop (C's `wavefront_heuristic_xdrop`): every `steps_between_cutoffs`
/// steps, trim diagonals from both ends of the wavefront while their
/// Smith-Waterman score is at least `xdrop` below the best seen so far (cells
/// with null offsets are skipped), then update the best. The first call only
/// records the best.
fn xdrop_cutoff(
    state: &mut HeuristicState,
    wf: &mut Wavefront,
    sw_scores: &mut [WfOffset],
    sw_base_k: i32,
    wf_match: i32,
    wf_score: i32,
    xdrop: i32,
    steps_between_cutoffs: i32,
) {
    if state.steps_wait > 0 {
        return;
    }
    let (cmax_sw_score, cmax_k, _) =
        compute_sw_scores(wf, sw_scores, sw_base_k, wf_match, wf_score);
    if state.max_sw_score_k != DPMATRIX_DIAGONAL_NULL {
        let max_sw_score = state.max_sw_score;
        let drops = |k: i32| max_sw_score - sw_scores[(k - sw_base_k) as usize] >= xdrop;
        // Reduce from bottom
        let mut k = wf.lo;
        while k <= wf.hi && (wf.get_offset(k) < 0 || drops(k)) {
            k += 1;
        }
        wf.lo = k;
        // Reduce from top
        let mut k = wf.hi;
        while k >= wf.lo && (wf.get_offset(k) < 0 || drops(k)) {
            k -= 1;
        }
        wf.hi = k;
        if cmax_sw_score > state.max_sw_score {
            state.max_sw_score = cmax_sw_score;
            state.max_sw_score_k = cmax_k;
        }
    } else {
        state.max_sw_score = cmax_sw_score;
        state.max_sw_score_k = cmax_k;
    }
    state.steps_wait = steps_between_cutoffs;
}

/// Z-Drop (C's `wavefront_heuristic_zdrop`): every `steps_between_cutoffs`
/// steps, compare the wavefront's best Smith-Waterman score with the best
/// seen. A new best is recorded (cell and the WFA score it was reached at);
/// otherwise, if it has dropped by more than `zdrop`, return `true`. The first
/// call only records the best.
fn zdrop_cutoff(
    state: &mut HeuristicState,
    wf: &mut Wavefront,
    sw_scores: &mut [WfOffset],
    sw_base_k: i32,
    wf_match: i32,
    wf_score: i32,
    zdrop: i32,
    steps_between_cutoffs: i32,
) -> bool {
    if state.steps_wait > 0 {
        return false;
    }
    let (cmax_sw_score, cmax_k, cmax_offset) =
        compute_sw_scores(wf, sw_scores, sw_base_k, wf_match, wf_score);
    if state.max_sw_score_k == DPMATRIX_DIAGONAL_NULL || cmax_sw_score > state.max_sw_score {
        state.max_sw_score = cmax_sw_score;
        state.max_wf_score = wf_score;
        state.max_sw_score_k = cmax_k;
        state.max_sw_score_offset = cmax_offset;
    } else if state.max_sw_score - cmax_sw_score > zdrop {
        return true; // Z-dropped
    }
    state.steps_wait = steps_between_cutoffs;
    false
}

/// Smith-Waterman score of each cell (C's `wf_heuristic_compute_sw_scores`):
/// `(swg_match * (v + h) - wf_score) / 2` with `swg_match = -wf_match`, or -1
/// when `wf_match == 0`. Cells with null offsets are skipped (left unset).
/// Returns the maximum score with its diagonal and offset.
fn compute_sw_scores(
    wf: &Wavefront,
    sw_scores: &mut [WfOffset],
    sw_base_k: i32,
    wf_match: i32,
    wf_score: i32,
) -> (i32, i32, WfOffset) {
    let swg_match = if wf_match != 0 { -wf_match } else { -1 };
    let (mut cmax_sw_score, mut cmax_k, mut cmax_offset) = (i32::MIN, 0, 0);
    for k in wf.lo..=wf.hi {
        let offset = wf.get_offset(k);
        if offset < 0 {
            continue;
        }
        let v = wavefront_v(k, offset);
        let h = wavefront_h(k, offset);
        let sw_score = (swg_match * (v + h) - wf_score) / 2;
        sw_scores[(k - sw_base_k) as usize] = sw_score;
        if cmax_sw_score < sw_score {
            cmax_sw_score = sw_score;
            cmax_k = k;
            cmax_offset = offset;
        }
    }
    (cmax_sw_score, cmax_k, cmax_offset)
}

/// Propagate M-wavefront's narrowed bounds to an I/D wavefront.
pub fn heuristic_equate(dst: &mut Wavefront, src: &Wavefront) {
    if src.lo > dst.lo {
        dst.lo = src.lo;
    }
    if src.hi < dst.hi {
        dst.hi = src.hi;
    }
    if dst.lo > dst.hi {
        dst.null = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wf_distance_end2end_basic() {
        // k=0, offset=5 → h=5, v=5. plen=10, tlen=10.
        // left_v = 5, left_h = 5 → distance = 5
        assert_eq!(wf_distance_end2end(5, 0, 10, 10), 5);
    }

    #[test]
    fn test_wf_distance_end2end_diagonal() {
        // k=2, offset=8 → h=8, v=6. plen=10, tlen=10.
        // left_v = 4, left_h = 2 → distance = 4
        assert_eq!(wf_distance_end2end(8, 2, 10, 10), 4);
    }

    #[test]
    fn test_wf_distance_end2end_null() {
        assert_eq!(wf_distance_end2end(OFFSET_NULL, 0, 10, 10), i32::MAX);
    }

    #[test]
    fn test_wf_distance_end2end_complete() {
        // k=0, offset=10 → h=10, v=10. plen=10, tlen=10.
        // left_v = 0, left_h = 0 → distance = 0
        assert_eq!(wf_distance_end2end(10, 0, 10, 10), 0);
    }

    #[test]
    fn test_wfadaptive_reduce_basic() {
        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(-3, 3);

        // Set distances in a buffer: middle diagonals close, edges far
        let mut distances = vec![0i32; 11];
        let base_k = -5;
        // k=-3: distance=100, k=-2: distance=5, k=-1..1: distance=2, k=2: distance=5, k=3: distance=100
        distances[(-3 - base_k) as usize] = 100;
        distances[(-2 - base_k) as usize] = 5;
        distances[(-1 - base_k) as usize] = 2;
        distances[(0 - base_k) as usize] = 2;
        distances[(1 - base_k) as usize] = 2;
        distances[(2 - base_k) as usize] = 5;
        distances[(3 - base_k) as usize] = 100;

        // min_distance=2, threshold=10, alignment_k=0
        wfadaptive_reduce(&mut wf, &distances, base_k, 2, 10, 0, 0);

        // k=-3 has distance 100, 100-2=98 > 10 → pruned
        // k=-2 has distance 5, 5-2=3 <= 10 → kept
        assert_eq!(wf.lo, -2);
        // k=3 has distance 100, 100-2=98 > 10 → pruned
        // k=2 has distance 5, 5-2=3 <= 10 → kept
        assert_eq!(wf.hi, 2);
    }

    #[test]
    fn test_wfadaptive_too_short() {
        // Wavefront length < min_wavefront_length → no reduction
        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(-1, 1); // length=3

        let mut distances = vec![0i32; 11];
        let mut state = HeuristicState::default();

        wfadaptive_cutoff(
            &mut state,
            &mut wf,
            &mut distances,
            -5,
            10,
            10,
            5, // min_wavefront_length=5, but wf length=3
            10,
            1,
        );

        assert_eq!(wf.lo, -1);
        assert_eq!(wf.hi, 1);
    }

    #[test]
    fn test_banded_static() {
        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(-3, 3);

        banded_static_cutoff(&mut wf, -2, 2);
        assert_eq!(wf.lo, -2);
        assert_eq!(wf.hi, 2);
    }

    #[test]
    fn test_compute_sw_scores_c_formula() {
        // C's wf_heuristic_compute_sw_scores: (swg_match * (v + h) - s) / 2 with
        // swg_match = -match, or -1 when match == 0; null offsets skipped.
        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(-1, 1);
        wf.set_offset(-1, 4); // v = 5, h = 4
        wf.set_offset(0, 6); // v = 6, h = 6
        wf.set_offset(1, OFFSET_NULL);
        let mut sw = vec![0; 11];
        let (max, k, off) = compute_sw_scores(&wf, &mut sw, -5, -3, 10);
        assert_eq!(sw[(-1 + 5) as usize], (3 * 9 - 10) / 2);
        assert_eq!((max, k, off), ((3 * 12 - 10) / 2, 0, 6));
        let (max, _, _) = compute_sw_scores(&wf, &mut sw, -5, 0, 10);
        assert_eq!(max, (-9 - 10) / 2); // match == 0: swg_match = -1
    }

    #[test]
    fn test_drop_state_reset_uses_steps() {
        let mut st = HeuristicState::default();
        st.reset(&HeuristicStrategy::XDrop {
            xdrop: 20,
            steps_between_cutoffs: 7,
        });
        assert_eq!(st.steps_wait, 7);
        assert_eq!(st.max_sw_score_k, DPMATRIX_DIAGONAL_NULL);
        st.reset(&HeuristicStrategy::ZDrop {
            zdrop: 20,
            steps_between_cutoffs: 3,
        });
        assert_eq!(st.steps_wait, 3);
    }

    #[test]
    fn test_banded_static_no_change() {
        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(-1, 1);

        banded_static_cutoff(&mut wf, -3, 3);
        assert_eq!(wf.lo, -1);
        assert_eq!(wf.hi, 1);
    }

    #[test]
    fn test_heuristic_equate() {
        let mut src = Wavefront::allocate(11, false);
        src.init(-5, 5);
        src.set_limits(-2, 2);

        let mut dst = Wavefront::allocate(11, false);
        dst.init(-5, 5);
        dst.set_limits(-4, 4);

        heuristic_equate(&mut dst, &src);
        assert_eq!(dst.lo, -2);
        assert_eq!(dst.hi, 2);
        assert!(!dst.null);
    }

    #[test]
    fn test_heuristic_equate_null() {
        let mut src = Wavefront::allocate(11, false);
        src.init(-5, 5);
        src.set_limits(2, 1); // lo > hi

        let mut dst = Wavefront::allocate(11, false);
        dst.init(-5, 5);
        dst.set_limits(-1, 1);

        heuristic_equate(&mut dst, &src);
        assert!(dst.null);
    }
}
