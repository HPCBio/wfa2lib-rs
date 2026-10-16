//! Wavefront extension: extends exact matches along diagonals.
//!
//! The extend kernel compares pattern and text characters starting from each
//! wavefront offset, advancing the offset while characters match. Sentinel
//! characters at the end of each sequence guarantee termination without
//! explicit bounds checks.

use crate::offset::{OFFSET_NULL, WfOffset};
use crate::sequences::WavefrontSequences;
use crate::termination;
use crate::wavefront::Wavefront;

/// Scalar extend kernel for a single diagonal.
#[inline(always)]
unsafe fn extend_scalar_one(
    pattern_base: *const u8,
    text_base: *const u8,
    k: i32,
    mut offset: WfOffset,
) -> WfOffset {
    unsafe {
        let v = (offset - k) as usize;
        let h = offset as usize;
        let mut p_ptr = pattern_base.add(v);
        let mut t_ptr = text_base.add(h);
        let mut cmp = (p_ptr as *const u64).read_unaligned()
            ^ (t_ptr as *const u64).read_unaligned();
        while cmp == 0 {
            offset += 8;
            p_ptr = p_ptr.add(8);
            t_ptr = t_ptr.add(8);
            cmp = (p_ptr as *const u64).read_unaligned()
                ^ (t_ptr as *const u64).read_unaligned();
        }
        offset += (cmp.trailing_zeros() / 8) as i32;
        offset
    }
}

/// Extend matches for end-to-end alignment using 64-bit blockwise comparison.
///
/// Scalar on every architecture. An AVX2 kernel that gathered 4 bytes for 8
/// diagonals at once (as C's) was removed: on x86 (Xeon Skylake-SP, EPYC Zen
/// 3) its gathers stalled, and long or divergent alignments ran up to 2× slower
/// than this loop (HPCBio/wfa2lib-rs#31).
#[inline]
pub fn extend_matches_packed_end2end(sequences: &WavefrontSequences, wf: &mut Wavefront) {
    let lo = wf.lo;
    let hi = wf.hi;
    let ptrs = sequences.ptrs();
    let pattern_base = ptrs.pattern;
    let text_base = ptrs.text;

    unsafe {
        let centered = wf.offsets_centered_mut_ptr();
        let mut off_ptr = centered.offset(lo as isize);
        let off_end = centered.offset(hi as isize).add(1);
        let mut k = lo;

        while off_ptr < off_end {
            let offset = *off_ptr;
            if offset >= 0 {
                *off_ptr = extend_scalar_one(pattern_base, text_base, k, offset);
            }
            off_ptr = off_ptr.add(1);
            k += 1;
        }
    }
}

/// Extend matches and return the maximum antidiagonal reached.
#[inline]
pub fn extend_matches_packed_end2end_max(
    sequences: &WavefrontSequences,
    wf: &mut Wavefront,
) -> WfOffset {
    let lo = wf.lo;
    let hi = wf.hi;
    let ptrs = sequences.ptrs();
    let pattern_base = ptrs.pattern;
    let text_base = ptrs.text;

    unsafe {
        let centered = wf.offsets_centered_mut_ptr();
        let mut off_ptr = centered.offset(lo as isize);
        let off_end = centered.offset(hi as isize).add(1);
        let mut k = lo;
        let mut max_antidiag: WfOffset = 0;

        while off_ptr < off_end {
            let offset = *off_ptr;
            if offset >= 0 {
                let extended = extend_scalar_one(pattern_base, text_base, k, offset);
                *off_ptr = extended;
                let antidiag = 2 * extended - k;
                if antidiag > max_antidiag {
                    max_antidiag = antidiag;
                }
            }
            off_ptr = off_ptr.add(1);
            k += 1;
        }

        max_antidiag
    }
}

/// Extend matches for ends-free alignment, checking termination at each diagonal.
///
/// Returns `Some((k, offset))` if a valid termination point is found, `None` otherwise.
#[inline]
pub fn extend_matches_packed_endsfree(
    sequences: &WavefrontSequences,
    wf: &mut Wavefront,
    pattern_end_free: i32,
    text_end_free: i32,
) -> Option<(i32, WfOffset)> {
    let lo = wf.lo;
    let hi = wf.hi;
    let ptrs = sequences.ptrs();
    let pattern_base = ptrs.pattern;
    let text_base = ptrs.text;
    let pattern_length = sequences.pattern_length;
    let text_length = sequences.text_length;

    // SAFETY: Same as extend_matches_packed_end2end.
    unsafe {
        let centered = wf.offsets_centered_mut_ptr();
        let mut off_ptr = centered.offset(lo as isize);
        let off_end = centered.offset(hi as isize).add(1);
        let mut k = lo;

        while off_ptr < off_end {
            let mut offset = *off_ptr;
            if offset != OFFSET_NULL {
                let v = (offset - k) as usize;
                let h = offset as usize;
                let mut p_ptr = pattern_base.add(v);
                let mut t_ptr = text_base.add(h);
                let mut cmp = (p_ptr as *const u64).read_unaligned()
                    ^ (t_ptr as *const u64).read_unaligned();
                while cmp == 0 {
                    offset += 8;
                    p_ptr = p_ptr.add(8);
                    t_ptr = t_ptr.add(8);
                    cmp = (p_ptr as *const u64).read_unaligned()
                        ^ (t_ptr as *const u64).read_unaligned();
                }
                offset += (cmp.trailing_zeros() / 8) as i32;
                *off_ptr = offset;

                // Check ends-free termination at this diagonal
                if termination::termination_endsfree(
                    pattern_length,
                    text_length,
                    pattern_end_free,
                    text_end_free,
                    k,
                    offset,
                ) {
                    return Some((k, offset));
                }
            }
            off_ptr = off_ptr.add(1);
            k += 1;
        }
    }

    None
}

/// Extend matches on every diagonal for ends-free alignment and return the
/// best termination point: the cell with the fewest bases left unaligned
/// (lowest k on ties), as `(k, offset, left)`. Unlike
/// [`extend_matches_packed_endsfree`], doesn't stop at the first one, so the
/// wavefront stays complete for later scores (corrected mode,
/// [`crate::aligner::Compat::Corrected`]).
#[inline]
pub fn extend_matches_packed_endsfree_best(
    sequences: &WavefrontSequences,
    wf: &mut Wavefront,
    pattern_end_free: i32,
    text_end_free: i32,
) -> Option<(i32, WfOffset, i32)> {
    let lo = wf.lo;
    let hi = wf.hi;
    let ptrs = sequences.ptrs();
    let pattern_base = ptrs.pattern;
    let text_base = ptrs.text;
    let pattern_length = sequences.pattern_length;
    let text_length = sequences.text_length;
    let mut best: Option<(i32, WfOffset, i32)> = None;

    // SAFETY: Same as extend_matches_packed_end2end.
    unsafe {
        let centered = wf.offsets_centered_mut_ptr();
        let mut off_ptr = centered.offset(lo as isize);
        let off_end = centered.offset(hi as isize).add(1);
        let mut k = lo;

        while off_ptr < off_end {
            let offset = *off_ptr;
            if offset != OFFSET_NULL {
                let offset = extend_scalar_one(pattern_base, text_base, k, offset);
                *off_ptr = offset;
                // Cheap check first: count the bases left only at an end
                if termination::termination_endsfree(
                    pattern_length,
                    text_length,
                    pattern_end_free,
                    text_end_free,
                    k,
                    offset,
                ) && let Some(left) = termination::endsfree_left(
                    pattern_length,
                    text_length,
                    pattern_end_free,
                    text_end_free,
                    k,
                    offset,
                ) && best.is_none_or(|(_, _, l)| left < l)
                {
                    best = Some((k, offset, left));
                }
            }
            off_ptr = off_ptr.add(1);
            k += 1;
        }
    }

    best
}

// --- Custom extend functions for Lambda mode (bounds-checked, no sentinels) ---

/// Extend matches for Lambda mode end-to-end alignment (bounds-checked).
#[inline]
pub fn extend_matches_custom_end2end(sequences: &WavefrontSequences, wf: &mut Wavefront) {
    let lo = wf.lo;
    let hi = wf.hi;
    let centered = unsafe { wf.offsets_centered_mut_ptr() };

    for k in lo..=hi {
        let mut offset = unsafe { *centered.offset(k as isize) };
        if offset == OFFSET_NULL {
            continue;
        }

        let mut v = offset - k;
        let mut h = offset;

        while sequences.cmp_lambda(v, h) {
            v += 1;
            h += 1;
            offset += 1;
        }

        unsafe { *centered.offset(k as isize) = offset; }
    }
}

/// Extend matches for Lambda mode and return the maximum antidiagonal reached.
#[inline]
pub fn extend_matches_custom_end2end_max(
    sequences: &WavefrontSequences,
    wf: &mut Wavefront,
) -> WfOffset {
    let lo = wf.lo;
    let hi = wf.hi;
    let centered = unsafe { wf.offsets_centered_mut_ptr() };
    let mut max_antidiag: WfOffset = 0;

    for k in lo..=hi {
        let mut offset = unsafe { *centered.offset(k as isize) };
        if offset == OFFSET_NULL {
            continue;
        }

        let mut v = offset - k;
        let mut h = offset;

        while sequences.cmp_lambda(v, h) {
            v += 1;
            h += 1;
            offset += 1;
        }

        unsafe { *centered.offset(k as isize) = offset; }

        let antidiag = 2 * offset - k;
        if antidiag > max_antidiag {
            max_antidiag = antidiag;
        }
    }

    max_antidiag
}

/// Extend matches for Lambda mode ends-free alignment, checking termination.
#[inline]
pub fn extend_matches_custom_endsfree(
    sequences: &WavefrontSequences,
    wf: &mut Wavefront,
    pattern_end_free: i32,
    text_end_free: i32,
) -> Option<(i32, WfOffset)> {
    let lo = wf.lo;
    let hi = wf.hi;
    let pattern_length = sequences.pattern_length;
    let text_length = sequences.text_length;
    let centered = unsafe { wf.offsets_centered_mut_ptr() };

    for k in lo..=hi {
        let mut offset = unsafe { *centered.offset(k as isize) };
        if offset == OFFSET_NULL {
            continue;
        }

        let mut v = offset - k;
        let mut h = offset;

        while sequences.cmp_lambda(v, h) {
            v += 1;
            h += 1;
            offset += 1;
        }

        unsafe { *centered.offset(k as isize) = offset; }

        if termination::termination_endsfree(
            pattern_length,
            text_length,
            pattern_end_free,
            text_end_free,
            k,
            offset,
        ) {
            return Some((k, offset));
        }
    }

    None
}

/// Lambda-mode [`extend_matches_packed_endsfree_best`].
#[inline]
pub fn extend_matches_custom_endsfree_best(
    sequences: &WavefrontSequences,
    wf: &mut Wavefront,
    pattern_end_free: i32,
    text_end_free: i32,
) -> Option<(i32, WfOffset, i32)> {
    let pattern_length = sequences.pattern_length;
    let text_length = sequences.text_length;
    let centered = unsafe { wf.offsets_centered_mut_ptr() };
    let mut best: Option<(i32, WfOffset, i32)> = None;

    for k in wf.lo..=wf.hi {
        let mut offset = unsafe { *centered.offset(k as isize) };
        if offset == OFFSET_NULL {
            continue;
        }

        while sequences.cmp_lambda(offset - k, offset) {
            offset += 1;
        }

        unsafe {
            *centered.offset(k as isize) = offset;
        }

        // Cheap check first: count the bases left only at an end
        if termination::termination_endsfree(
            pattern_length,
            text_length,
            pattern_end_free,
            text_end_free,
            k,
            offset,
        ) && let Some(left) = termination::endsfree_left(
            pattern_length,
            text_length,
            pattern_end_free,
            text_end_free,
            k,
            offset,
        ) && best.is_none_or(|(_, _, l)| left < l)
        {
            best = Some((k, offset, left));
        }
    }

    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extend_identical() {
        let mut seq = WavefrontSequences::new();
        seq.init_ascii(b"ACGT", b"ACGT", false);

        let mut wf = Wavefront::allocate(7, false);
        wf.init(-3, 3);
        wf.set_limits(0, 0);
        wf.set_offset(0, 0);

        extend_matches_packed_end2end(&seq, &mut wf);

        // Should extend all 4 characters
        assert_eq!(wf.get_offset(0), 4);
    }

    #[test]
    fn test_extend_no_match() {
        let mut seq = WavefrontSequences::new();
        seq.init_ascii(b"AAAA", b"TTTT", false);

        let mut wf = Wavefront::allocate(7, false);
        wf.init(-3, 3);
        wf.set_limits(0, 0);
        wf.set_offset(0, 0);

        extend_matches_packed_end2end(&seq, &mut wf);

        // No characters match at position 0
        assert_eq!(wf.get_offset(0), 0);
    }

    #[test]
    fn test_extend_partial_match() {
        let mut seq = WavefrontSequences::new();
        seq.init_ascii(b"ACGT", b"ACTT", false);

        let mut wf = Wavefront::allocate(7, false);
        wf.init(-3, 3);
        wf.set_limits(0, 0);
        wf.set_offset(0, 0);

        extend_matches_packed_end2end(&seq, &mut wf);

        // Matches: A, C then stops (G != T)
        assert_eq!(wf.get_offset(0), 2);
    }

    #[test]
    fn test_extend_null_offset_skipped() {
        let mut seq = WavefrontSequences::new();
        seq.init_ascii(b"ACGT", b"ACGT", false);

        let mut wf = Wavefront::allocate(7, false);
        wf.init(-3, 3);
        wf.set_limits(-1, 1);
        wf.set_offset(-1, OFFSET_NULL);
        wf.set_offset(0, 0);
        wf.set_offset(1, OFFSET_NULL);

        extend_matches_packed_end2end(&seq, &mut wf);

        // Only k=0 should be extended
        assert_eq!(wf.get_offset(0), 4);
        assert_eq!(wf.get_offset(-1), OFFSET_NULL);
        assert_eq!(wf.get_offset(1), OFFSET_NULL);
    }

    // --- Ends-free extend tests ---

    #[test]
    fn test_endsfree_terminates_on_text_end() {
        // pattern="ACGT" (4), text="ACGTXXXX" (8)
        // With text_end_free=4, once pattern is fully consumed (v=4),
        // remaining text (4) <= text_end_free → terminate
        let mut seq = WavefrontSequences::new();
        seq.init_ascii(b"ACGT", b"ACGTXXXX", false);

        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(0, 0);
        wf.set_offset(0, 0);

        // k=0, after extending: offset=4, v=4=pattern_length
        // text_left = 8 - 4 = 4 <= text_end_free=4 → terminate
        let result = extend_matches_packed_endsfree(&seq, &mut wf, 0, 4);
        assert!(result.is_some());
        let (k, offset) = result.unwrap();
        assert_eq!(k, 0);
        assert_eq!(offset, 4);
    }

    #[test]
    fn test_endsfree_best_fewest_left() {
        // pattern="AAAA", text="AAAAAA": diagonals 0, 1, 2 all consume the
        // pattern, leaving 2, 1 and 0 text bases. The first-hit kernel stops
        // at k=0; the best kernel extends all and picks k=2.
        let check = |seq: &WavefrontSequences, lambda: bool| {
            let mut wf = Wavefront::allocate(11, false);
            wf.init(-5, 5);
            wf.set_limits(0, 2);
            for k in 0..=2 {
                wf.set_offset(k, k);
            }
            let best = if lambda {
                extend_matches_custom_endsfree_best(seq, &mut wf, 0, 6)
            } else {
                extend_matches_packed_endsfree_best(seq, &mut wf, 0, 6)
            };
            assert_eq!(best, Some((2, 6, 0)));
            assert_eq!((wf.get_offset(0), wf.get_offset(1)), (4, 5));
        };
        let mut seq = WavefrontSequences::new();
        seq.init_ascii(b"AAAA", b"AAAAAA", false);
        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(0, 2);
        for k in 0..=2 {
            wf.set_offset(k, k);
        }
        assert_eq!(extend_matches_packed_endsfree(&seq, &mut wf, 0, 6), Some((0, 4)));
        check(&seq, false);
        let mut seq = WavefrontSequences::new();
        seq.init_lambda(Box::new(|_, _| true), 4, 6, false);
        check(&seq, true);
    }

    #[test]
    fn test_endsfree_no_termination() {
        // pattern="ACGTACGT" (8), text="ACGTXXXX" (8)
        // With pattern_end_free=0, text_end_free=0 (same as end2end)
        // After extending: offset=4 (matches ACGT then stops)
        // h=4 < 8, v=4 < 8 → no termination
        let mut seq = WavefrontSequences::new();
        seq.init_ascii(b"ACGTACGT", b"ACGTXXXX", false);

        let mut wf = Wavefront::allocate(11, false);
        wf.init(-5, 5);
        wf.set_limits(0, 0);
        wf.set_offset(0, 0);

        let result = extend_matches_packed_endsfree(&seq, &mut wf, 0, 0);
        assert!(result.is_none());
        assert_eq!(wf.get_offset(0), 4);
    }
}
