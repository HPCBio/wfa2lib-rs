//! BiWFA: Bidirectional Wavefront Alignment.
//!
//! Implements divide-and-conquer bidirectional alignment for O(s) memory usage.
//! Forward and reverse wavefronts expand toward each other until they overlap,
//! producing a breakpoint that splits the problem into two halves.

use crate::offset::{WfOffset, wavefront_k_inverse, wavefront_v};
use crate::wavefront::Wavefront;

/// DP matrix component type for BiWFA breakpoints.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComponentType {
    M,
    I1,
    D1,
    I2,
    D2,
}

/// Breakpoint where forward and reverse wavefronts meet.
#[derive(Debug, Clone)]
pub struct BiAlignBreakpoint {
    /// Best total score found so far (i32::MAX = no breakpoint yet).
    pub score: i32,
    pub score_forward: i32,
    pub score_reverse: i32,
    pub k_forward: i32,
    pub k_reverse: i32,
    pub offset_forward: WfOffset,
    pub offset_reverse: WfOffset,
    pub component: ComponentType,
}

impl Default for BiAlignBreakpoint {
    fn default() -> Self {
        Self::new()
    }
}

impl BiAlignBreakpoint {
    pub fn new() -> Self {
        Self {
            score: i32::MAX,
            score_forward: 0,
            score_reverse: 0,
            k_forward: 0,
            k_reverse: 0,
            offset_forward: 0,
            offset_reverse: 0,
            component: ComponentType::M,
        }
    }

    pub fn found(&self) -> bool {
        self.score < i32::MAX
    }
}

/// Get the gap opening penalty for overlap early-exit checks.
///
/// Returns the minimum gap opening penalty for the distance metric.
/// For edit/indel/linear: 0 (no gap opening concept or it's 0).
/// For affine: gap_opening1.
/// For affine2p: min(gap_opening1, gap_opening2).
pub fn overlap_gopen_adjust(gap_opening1: i32, gap_opening2: i32, has_affine2p: bool) -> i32 {
    if has_affine2p {
        gap_opening1.min(gap_opening2)
    } else {
        gap_opening1
    }
}

/// Check M-to-M wavefront overlap between forward (wf_0) and reverse (wf_1) wavefronts.
///
/// If `h_0 + h_1 >= text_length` for any overlapping diagonal, a breakpoint is found.
/// `breakpoint_forward`: if true, wf_0 is the forward aligner; if false, wf_0 is reverse.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub fn breakpoint_m2m(
    wf_0: &Wavefront,
    wf_1: &Wavefront,
    pattern_length: i32,
    text_length: i32,
    score_0: i32,
    score_1: i32,
    breakpoint_forward: bool,
    breakpoint: &mut BiAlignBreakpoint,
) {
    if score_0 + score_1 >= breakpoint.score {
        return;
    }

    let lo_0 = wf_0.lo;
    let hi_0 = wf_0.hi;

    // Convert reverse k-range to forward coordinates
    let lo_1_inv = wavefront_k_inverse(wf_1.hi, pattern_length, text_length);
    let hi_1_inv = wavefront_k_inverse(wf_1.lo, pattern_length, text_length);

    // Check k-range overlap
    if hi_1_inv < lo_0 || hi_0 < lo_1_inv {
        return;
    }

    let max_lo = lo_0.max(lo_1_inv);
    let min_hi = hi_0.min(hi_1_inv);

    // Use k-centered raw pointers: ptr[k] directly accesses diagonal k.
    // SAFETY: k_0 ∈ [max_lo, min_hi] ⊆ [lo_0, hi_0] ⊆ allocated range of wf_0.
    //         k_1 ∈ [wf_1.lo, wf_1.hi] ⊆ allocated range of wf_1 (by construction of lo_1_inv/hi_1_inv).
    let ptr_0 = unsafe { wf_0.offsets_centered_ptr() };
    let ptr_1 = unsafe { wf_1.offsets_centered_ptr() };

    // SAFETY: k_0 ∈ [max_lo, min_hi] ⊆ [lo_0, hi_0] and k_1 ∈ [wf_1.lo, wf_1.hi].
    unsafe {
        scan_overlap(
            ptr_0,
            ptr_1,
            max_lo,
            min_hi,
            pattern_length,
            text_length,
            |k_0, k_1, offset_0, offset_1| {
                if breakpoint_forward {
                    breakpoint.score_forward = score_0;
                    breakpoint.score_reverse = score_1;
                    breakpoint.k_forward = k_0;
                    breakpoint.k_reverse = k_1;
                    breakpoint.offset_forward = offset_0;
                    breakpoint.offset_reverse = offset_1;
                } else {
                    breakpoint.score_forward = score_1;
                    breakpoint.score_reverse = score_0;
                    breakpoint.k_forward = k_1;
                    breakpoint.k_reverse = k_0;
                    breakpoint.offset_forward = offset_1;
                    breakpoint.offset_reverse = offset_0;
                }
                breakpoint.score = score_0 + score_1;
                breakpoint.component = ComponentType::M;
                true
            },
        );
    }
}

/// Check I/D-to-I/D wavefront overlap between forward and reverse wavefronts.
///
/// Same as m2m but subtracts gap_opening from the combined score.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub fn breakpoint_indel2indel(
    wf_0: &Wavefront,
    wf_1: &Wavefront,
    pattern_length: i32,
    text_length: i32,
    score_0: i32,
    score_1: i32,
    gap_opening: i32,
    component: ComponentType,
    breakpoint_forward: bool,
    breakpoint: &mut BiAlignBreakpoint,
) {
    if score_0 + score_1 - gap_opening >= breakpoint.score {
        return;
    }

    let lo_0 = wf_0.lo;
    let hi_0 = wf_0.hi;

    let lo_1_inv = wavefront_k_inverse(wf_1.hi, pattern_length, text_length);
    let hi_1_inv = wavefront_k_inverse(wf_1.lo, pattern_length, text_length);

    if hi_1_inv < lo_0 || hi_0 < lo_1_inv {
        return;
    }

    let max_lo = lo_0.max(lo_1_inv);
    let min_hi = hi_0.min(hi_1_inv);

    // Use k-centered raw pointers: ptr[k] directly accesses diagonal k.
    // SAFETY: same argument as breakpoint_m2m — k_0 ∈ [lo_0, hi_0] and k_1 ∈ [wf_1.lo, wf_1.hi].
    let ptr_0 = unsafe { wf_0.offsets_centered_ptr() };
    let ptr_1 = unsafe { wf_1.offsets_centered_ptr() };

    unsafe {
        scan_overlap(
            ptr_0,
            ptr_1,
            max_lo,
            min_hi,
            pattern_length,
            text_length,
            |k_0, k_1, offset_0, offset_1| {
                // Extra bounds check for indel breakpoints
                let (check_k, check_h) = if breakpoint_forward {
                    (k_0, offset_0)
                } else {
                    (k_1, offset_1)
                };
                let v = wavefront_v(check_k, check_h);
                if v > pattern_length || check_h > text_length {
                    return false;
                }

                if breakpoint_forward {
                    breakpoint.score_forward = score_0;
                    breakpoint.score_reverse = score_1;
                    breakpoint.k_forward = k_0;
                    breakpoint.k_reverse = k_1;
                    breakpoint.offset_forward = offset_0;
                    breakpoint.offset_reverse = offset_1;
                } else {
                    breakpoint.score_forward = score_1;
                    breakpoint.score_reverse = score_0;
                    breakpoint.k_forward = k_1;
                    breakpoint.k_reverse = k_0;
                    breakpoint.offset_forward = offset_1;
                    breakpoint.offset_reverse = offset_0;
                }
                breakpoint.score = score_0 + score_1 - gap_opening;
                breakpoint.component = component;
                true
            },
        );
    }
}

/// Diagonals checked per vector step in [`scan_overlap`].
const SCAN_CHUNK: i32 = 8;

/// Scan the overlap `k_0 ∈ [lo, hi]` in ascending order for diagonals where
/// `ptr_0[k_0] + ptr_1[k_1] >= text_length` (`k_1` is `k_0` in the other
/// direction's coordinates), calling `accept` on each until it returns true.
///
/// Chunks of [`SCAN_CHUNK`] diagonals are first tested together with SIMD;
/// only a chunk with a candidate is re-checked one diagonal at a time. The
/// visiting order is unchanged, so the breakpoint chosen is the same as a
/// plain scalar scan (and as C's scalar and AVX2 versions).
///
/// # Safety
/// `ptr_0[k]` for `k ∈ [lo, hi]` and `ptr_1[k_inverse(k)]` must be readable.
#[inline(always)]
unsafe fn scan_overlap(
    ptr_0: *const WfOffset,
    ptr_1: *const WfOffset,
    lo: i32,
    hi: i32,
    pattern_length: i32,
    text_length: i32,
    mut accept: impl FnMut(i32, i32, WfOffset, WfOffset) -> bool,
) -> bool {
    let mut k_0 = lo;
    #[cfg(target_arch = "x86_64")]
    let use_avx2 = cfg!(target_feature = "avx2") || std::arch::is_x86_feature_detected!("avx2");
    while hi - k_0 + 1 >= SCAN_CHUNK {
        // The chunk's k_1 values run down from k_1_top: they sit at
        // ptr_1[k_1_top - 7 ..= k_1_top] in reverse order.
        let k_1_top = wavefront_k_inverse(k_0, pattern_length, text_length);
        let p_0 = unsafe { ptr_0.offset(k_0 as isize) };
        let p_1 = unsafe { ptr_1.offset((k_1_top - (SCAN_CHUNK - 1)) as isize) };
        #[cfg(target_arch = "aarch64")]
        let any = unsafe { chunk_any_neon(p_0, p_1, text_length) };
        #[cfg(target_arch = "x86_64")]
        let any = if use_avx2 {
            unsafe { chunk_any_avx2(p_0, p_1, text_length) }
        } else {
            unsafe { chunk_any_scalar(p_0, p_1, text_length) }
        };
        #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
        let any = unsafe { chunk_any_scalar(p_0, p_1, text_length) };
        if any {
            for i in 0..SCAN_CHUNK {
                let offset_0 = unsafe { *ptr_0.offset((k_0 + i) as isize) };
                let offset_1 = unsafe { *ptr_1.offset((k_1_top - i) as isize) };
                if offset_0 + offset_1 >= text_length
                    && accept(k_0 + i, k_1_top - i, offset_0, offset_1)
                {
                    return true;
                }
            }
        }
        k_0 += SCAN_CHUNK;
    }
    while k_0 <= hi {
        let k_1 = wavefront_k_inverse(k_0, pattern_length, text_length);
        let offset_0 = unsafe { *ptr_0.offset(k_0 as isize) };
        let offset_1 = unsafe { *ptr_1.offset(k_1 as isize) };
        if offset_0 + offset_1 >= text_length && accept(k_0, k_1, offset_0, offset_1) {
            return true;
        }
        k_0 += 1;
    }
    false
}

/// Whether any `p_0[i] + p_1[7 - i] >= text_length` for `i ∈ 0..8`.
#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
#[inline(always)]
unsafe fn chunk_any_scalar(p_0: *const WfOffset, p_1: *const WfOffset, text_length: i32) -> bool {
    let mut any = false;
    for i in 0..SCAN_CHUNK as usize {
        any |= unsafe { *p_0.add(i) + *p_1.add(SCAN_CHUNK as usize - 1 - i) } >= text_length;
    }
    any
}

#[cfg(target_arch = "aarch64")]
#[inline(always)]
unsafe fn chunk_any_neon(p_0: *const WfOffset, p_1: *const WfOffset, text_length: i32) -> bool {
    use std::arch::aarch64::*;
    unsafe {
        let v_tlen = vdupq_n_s32(text_length);
        // Reverse each 4-lane half of p_1; the high half pairs with the low half of p_0.
        let rev = |v: int32x4_t| {
            let r = vrev64q_s32(v);
            vcombine_s32(vget_high_s32(r), vget_low_s32(r))
        };
        let lo_1 = rev(vld1q_s32(p_1.add(4)));
        let hi_1 = rev(vld1q_s32(p_1));
        let m_lo = vcgeq_s32(vaddq_s32(vld1q_s32(p_0), lo_1), v_tlen);
        let m_hi = vcgeq_s32(vaddq_s32(vld1q_s32(p_0.add(4)), hi_1), v_tlen);
        vmaxvq_u32(vorrq_u32(m_lo, m_hi)) != 0
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn chunk_any_avx2(p_0: *const WfOffset, p_1: *const WfOffset, text_length: i32) -> bool {
    use std::arch::x86_64::*;
    unsafe {
        let rev = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
        let v_0 = _mm256_loadu_si256(p_0 as *const __m256i);
        let v_1 = _mm256_permutevar8x32_epi32(_mm256_loadu_si256(p_1 as *const __m256i), rev);
        let mask = _mm256_cmpgt_epi32(
            _mm256_add_epi32(v_0, v_1),
            _mm256_set1_epi32(text_length - 1),
        );
        _mm256_movemask_epi8(mask) != 0
    }
}

/// Wavefront access helper: get a wavefront from a raw pointer.
/// Returns None if the pointer is null or the wavefront is null.
fn get_wf_if_valid(ptr: *mut Wavefront) -> Option<&'static Wavefront> {
    if ptr.is_null() {
        return None;
    }
    let wf = unsafe { &*ptr };
    if wf.null {
        return None;
    }
    Some(wf)
}

/// Check all wavefront overlaps between two aligners at given scores.
///
/// `aligner_0` is at `score_0`, `aligner_1` is checked from `score_1` downward.
/// `breakpoint_forward`: true if aligner_0 is the forward aligner.
#[allow(clippy::too_many_arguments)]
pub fn bialign_overlap<const N: usize>(
    components_0: &crate::components::WavefrontComponents<N>,
    _slab_0: &crate::slab::WavefrontSlab,
    components_1: &crate::components::WavefrontComponents<N>,
    _slab_1: &crate::slab::WavefrontSlab,
    score_0: i32,
    score_1: i32,
    pattern_length: i32,
    text_length: i32,
    gap_opening1: i32,
    gap_opening2: i32,
    has_affine: bool,
    has_affine2p: bool,
    breakpoint_forward: bool,
    breakpoint: &mut BiAlignBreakpoint,
) {
    let max_score_scope = components_1.max_score_scope;

    // Get M-wavefront at score_0
    let mwf_0 = match get_wf_if_valid(components_0.get_m_ptr(score_0 as usize)) {
        Some(wf) => wf,
        None => return,
    };

    // Also get I/D wavefronts at score_0 (for indel-to-indel checks)
    let d1wf_0 = if has_affine {
        get_wf_if_valid(components_0.get_d1_ptr(score_0 as usize))
    } else {
        None
    };
    let i1wf_0 = if has_affine {
        get_wf_if_valid(components_0.get_i1_ptr(score_0 as usize))
    } else {
        None
    };
    let d2wf_0 = if has_affine2p {
        get_wf_if_valid(components_0.get_d2_ptr(score_0 as usize))
    } else {
        None
    };
    let i2wf_0 = if has_affine2p {
        get_wf_if_valid(components_0.get_i2_ptr(score_0 as usize))
    } else {
        None
    };

    // Iterate backwards through scores of aligner_1
    for i in 0..max_score_scope {
        let score_i = score_1 - i as i32;
        if score_i < 0 {
            break;
        }

        // Check D2/I2 breakpoints (affine2p only)
        if has_affine2p && score_0 + score_i - gap_opening2 < breakpoint.score {
            if let (Some(d2wf_0), Some(d2wf_1)) = (
                d2wf_0,
                get_wf_if_valid(components_1.get_d2_ptr(score_i as usize)),
            ) {
                breakpoint_indel2indel(
                    d2wf_0,
                    d2wf_1,
                    pattern_length,
                    text_length,
                    score_0,
                    score_i,
                    gap_opening2,
                    ComponentType::D2,
                    breakpoint_forward,
                    breakpoint,
                );
            }
            if let (Some(i2wf_0), Some(i2wf_1)) = (
                i2wf_0,
                get_wf_if_valid(components_1.get_i2_ptr(score_i as usize)),
            ) {
                breakpoint_indel2indel(
                    i2wf_0,
                    i2wf_1,
                    pattern_length,
                    text_length,
                    score_0,
                    score_i,
                    gap_opening2,
                    ComponentType::I2,
                    breakpoint_forward,
                    breakpoint,
                );
            }
        }

        // Check D1/I1 breakpoints (affine and above)
        if has_affine && score_0 + score_i - gap_opening1 < breakpoint.score {
            if let (Some(d1wf_0), Some(d1wf_1)) = (
                d1wf_0,
                get_wf_if_valid(components_1.get_d1_ptr(score_i as usize)),
            ) {
                breakpoint_indel2indel(
                    d1wf_0,
                    d1wf_1,
                    pattern_length,
                    text_length,
                    score_0,
                    score_i,
                    gap_opening1,
                    ComponentType::D1,
                    breakpoint_forward,
                    breakpoint,
                );
            }
            if let (Some(i1wf_0), Some(i1wf_1)) = (
                i1wf_0,
                get_wf_if_valid(components_1.get_i1_ptr(score_i as usize)),
            ) {
                breakpoint_indel2indel(
                    i1wf_0,
                    i1wf_1,
                    pattern_length,
                    text_length,
                    score_0,
                    score_i,
                    gap_opening1,
                    ComponentType::I1,
                    breakpoint_forward,
                    breakpoint,
                );
            }
        }

        // Check M-to-M breakpoints (all metrics)
        if score_0 + score_i >= breakpoint.score {
            continue;
        }
        if let Some(mwf_1) = get_wf_if_valid(components_1.get_m_ptr(score_i as usize)) {
            breakpoint_m2m(
                mwf_0,
                mwf_1,
                pattern_length,
                text_length,
                score_0,
                score_i,
                breakpoint_forward,
                breakpoint,
            );
        }
    }
}

/// Fallback score threshold: below this, use standard unidirectional WFA.
pub const BIALIGN_FALLBACK_MIN_SCORE: i32 = 250;

/// Sequences no longer than this are aligned with unidirectional WFA
/// (C's `WF_BIALIGN_FALLBACK_MIN_LENGTH`; top level only).
pub const BIALIGN_FALLBACK_MIN_LENGTH: i32 = 100;

/// When one direction reaches the end before the wavefronts overlap, fall
/// back to unidirectional WFA if the end score is at most this
/// (C's `WF_BIALIGN_RECOVERY_MIN_SCORE`); otherwise the end is unreachable.
pub const BIALIGN_RECOVERY_MIN_SCORE: i32 = 500;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wavefront::Wavefront;

    #[test]
    fn test_breakpoint_new() {
        let bp = BiAlignBreakpoint::new();
        assert!(!bp.found());
        assert_eq!(bp.score, i32::MAX);
    }

    #[test]
    fn test_breakpoint_m2m_found() {
        // Forward: k=0, offset=6, h=6. Reverse: k_inv=0, offset=5, h=5.
        // plen=10, tlen=10. h_0 + h_1 = 11 >= 10 → breakpoint
        let mut wf_fwd = Wavefront::allocate(7, false);
        wf_fwd.init(-3, 3);
        wf_fwd.set_limits(0, 0);
        wf_fwd.set_offset(0, 6);

        let mut wf_rev = Wavefront::allocate(7, false);
        wf_rev.init(-3, 3);
        wf_rev.set_limits(0, 0);
        wf_rev.set_offset(0, 5);

        let mut bp = BiAlignBreakpoint::new();
        breakpoint_m2m(&wf_fwd, &wf_rev, 10, 10, 3, 2, true, &mut bp);

        assert!(bp.found());
        assert_eq!(bp.score, 5);
        assert_eq!(bp.score_forward, 3);
        assert_eq!(bp.score_reverse, 2);
        assert_eq!(bp.k_forward, 0);
        assert_eq!(bp.offset_forward, 6);
        assert_eq!(bp.component, ComponentType::M);
    }

    #[test]
    fn test_breakpoint_m2m_not_found() {
        // h_0 + h_1 = 4 + 4 = 8 < 10 → no breakpoint
        let mut wf_fwd = Wavefront::allocate(7, false);
        wf_fwd.init(-3, 3);
        wf_fwd.set_limits(0, 0);
        wf_fwd.set_offset(0, 4);

        let mut wf_rev = Wavefront::allocate(7, false);
        wf_rev.init(-3, 3);
        wf_rev.set_limits(0, 0);
        wf_rev.set_offset(0, 4);

        let mut bp = BiAlignBreakpoint::new();
        breakpoint_m2m(&wf_fwd, &wf_rev, 10, 10, 3, 2, true, &mut bp);

        assert!(!bp.found());
    }

    #[test]
    fn test_breakpoint_m2m_different_lengths() {
        // plen=8, tlen=12 → alignment_k = 4
        // Forward: k=2, offset=9 → h=9
        // Reverse: k_inv for k=2: (12-8)-2 = 2. So reverse at k=2, offset=4 → h=4
        // h_0 + h_1 = 9 + 4 = 13 >= 12 → breakpoint
        let mut wf_fwd = Wavefront::allocate(11, false);
        wf_fwd.init(-5, 5);
        wf_fwd.set_limits(2, 2);
        wf_fwd.set_offset(2, 9);

        let mut wf_rev = Wavefront::allocate(11, false);
        wf_rev.init(-5, 5);
        wf_rev.set_limits(2, 2);
        wf_rev.set_offset(2, 4);

        let mut bp = BiAlignBreakpoint::new();
        breakpoint_m2m(&wf_fwd, &wf_rev, 8, 12, 5, 3, true, &mut bp);

        assert!(bp.found());
        assert_eq!(bp.score, 8);
    }

    #[test]
    fn test_breakpoint_m2m_no_k_overlap() {
        // Forward covers k=[0,2], reverse covers k=[-3,-1] → k_inv=[-1,1] for plen=tlen=10
        // Actually k_inv of k=[-3,-1] = (10-10)-k = [1,3]
        // Forward [0,2], reverse inv [1,3] → overlap [1,2]
        let mut wf_fwd = Wavefront::allocate(11, false);
        wf_fwd.init(-5, 5);
        wf_fwd.set_limits(4, 5); // k=[4,5]

        let mut wf_rev = Wavefront::allocate(11, false);
        wf_rev.init(-5, 5);
        wf_rev.set_limits(-5, -4); // k=[-5,-4], inv=[4,5] for plen=tlen=10
        wf_rev.set_offset(-5, 1);
        wf_rev.set_offset(-4, 1);

        wf_fwd.set_offset(4, 1);
        wf_fwd.set_offset(5, 1);

        let mut bp = BiAlignBreakpoint::new();
        breakpoint_m2m(&wf_fwd, &wf_rev, 10, 10, 3, 2, true, &mut bp);

        // h_0 + h_1 = 1 + 1 = 2 < 10 → no breakpoint
        assert!(!bp.found());
    }

    #[test]
    fn test_breakpoint_indel2indel_found() {
        // Same as m2m test but with gap_opening subtraction
        let mut wf_fwd = Wavefront::allocate(7, false);
        wf_fwd.init(-3, 3);
        wf_fwd.set_limits(0, 0);
        wf_fwd.set_offset(0, 6);

        let mut wf_rev = Wavefront::allocate(7, false);
        wf_rev.init(-3, 3);
        wf_rev.set_limits(0, 0);
        wf_rev.set_offset(0, 5);

        let mut bp = BiAlignBreakpoint::new();
        breakpoint_indel2indel(
            &wf_fwd,
            &wf_rev,
            10,
            10,
            3,
            2,
            6, // gap_opening
            ComponentType::D1,
            true,
            &mut bp,
        );

        assert!(bp.found());
        // score = 3 + 2 - 6 = -1
        assert_eq!(bp.score, -1);
        assert_eq!(bp.component, ComponentType::D1);
    }

    #[test]
    fn test_overlap_gopen_adjust() {
        assert_eq!(overlap_gopen_adjust(6, 24, false), 6);
        assert_eq!(overlap_gopen_adjust(6, 24, true), 6);
        assert_eq!(overlap_gopen_adjust(6, 3, true), 3);
    }

    /// Naive reference for `scan_overlap`: every hit in ascending order.
    #[allow(clippy::too_many_arguments)]
    fn naive_hits(
        a_0: &[i32],
        a_1: &[i32],
        c_0: i32,
        c_1: i32,
        lo: i32,
        hi: i32,
        plen: i32,
        tlen: i32,
    ) -> Vec<i32> {
        (lo..=hi)
            .filter(|&k| {
                let k_1 = wavefront_k_inverse(k, plen, tlen);
                a_0[(c_0 + k) as usize] + a_1[(c_1 + k_1) as usize] >= tlen
            })
            .collect()
    }

    #[test]
    fn test_scan_overlap_matches_scalar() {
        use crate::offset::OFFSET_NULL;
        // Deterministic pseudo-random offsets, some null, around tlen/2.
        let mut seed = 0x2545_f491u32;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 17;
            seed ^= seed << 5;
            seed
        };
        let (plen, tlen) = (90, 100);
        let n = 80usize;
        let (c_0, c_1) = (40i32, 40i32);
        for trial in 0..400 {
            let mut offsets = || -> Vec<i32> {
                (0..n)
                    .map(|_| {
                        let r = next();
                        if r % 7 == 0 {
                            OFFSET_NULL
                        } else {
                            30 + (r % 25) as i32 + (trial % 3) * 3
                        }
                    })
                    .collect()
            };
            let (a_0, a_1) = (offsets(), offsets());
            // Overlap lengths 1..=26 cover pure tails, whole chunks and both.
            let lo = -12 + (next() % 8) as i32;
            let hi = lo + (next() % 26) as i32;
            let hits = naive_hits(&a_0, &a_1, c_0, c_1, lo, hi, plen, tlen);
            // Reject the first `skip` hits, so a rejected hit inside a chunk
            // must not end or skip the rest of that chunk.
            for skip in 0..3 {
                let mut seen = Vec::new();
                let found = unsafe {
                    scan_overlap(
                        a_0.as_ptr().offset(c_0 as isize),
                        a_1.as_ptr().offset(c_1 as isize),
                        lo,
                        hi,
                        plen,
                        tlen,
                        |k_0, k_1, o_0, o_1| {
                            assert_eq!(k_1, wavefront_k_inverse(k_0, plen, tlen));
                            assert_eq!(o_0, a_0[(c_0 + k_0) as usize]);
                            assert_eq!(o_1, a_1[(c_1 + k_1) as usize]);
                            seen.push(k_0);
                            seen.len() > skip
                        },
                    )
                };
                let want: Vec<i32> = hits.iter().copied().take(skip + 1).collect();
                assert_eq!(seen, want, "trial {trial} lo {lo} hi {hi} skip {skip}");
                assert_eq!(found, hits.len() > skip);
            }
        }
    }

    #[test]
    fn test_chunk_any_simd_matches_scalar() {
        let a_0 = [10, 20, 30, 40, 50, 60, 70, 80];
        let tlen = 100;
        for lane in 0..8 {
            for delta in [-1, 0, 1] {
                // Only pair (lane, 7 - lane) can reach tlen.
                let mut a_1 = [0i32; 8];
                a_1[7 - lane] = tlen - a_0[lane] + delta;
                let want = unsafe { chunk_any_scalar(a_0.as_ptr(), a_1.as_ptr(), tlen) };
                assert_eq!(want, delta >= 0);
                #[cfg(target_arch = "aarch64")]
                assert_eq!(
                    unsafe { chunk_any_neon(a_0.as_ptr(), a_1.as_ptr(), tlen) },
                    want
                );
                #[cfg(target_arch = "x86_64")]
                if std::arch::is_x86_feature_detected!("avx2") {
                    assert_eq!(
                        unsafe { chunk_any_avx2(a_0.as_ptr(), a_1.as_ptr(), tlen) },
                        want
                    );
                }
            }
        }
    }
}
