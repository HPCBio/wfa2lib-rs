//! Compute kernels with piggyback backtrace (C's `*_idm_piggyback`), used
//! in the `med` and `low` memory modes.
//!
//! Alongside each offset, every cell carries a packed CIGAR (`pcigar`, up to
//! 16 operations) and the index of the backtrace block it continues
//! (`bt_prev`). Each kernel copies these from the winning predecessor and
//! pushes the operation taken, so the alignment can be recovered from the
//! end cell alone. Matches are not recorded; unpacking re-extends them.
//! Offsets are computed exactly as in the plain kernels, and ties are broken
//! as in C (later checks win).

#![allow(unsafe_op_in_unsafe_fn)]

use crate::bt_buffer::BacktraceBuffer;
use crate::offset::{OFFSET_NULL, WfOffset, wavefront_h, wavefront_v};
use crate::pcigar::{
    PCIGAR_ALMOST_FULL_MASK, PCIGAR_DELETION, PCIGAR_FULL_MASK, PCIGAR_HALF_FULL_MASK,
    PCIGAR_INSERTION, PCIGAR_MAX_LENGTH, PCIGAR_MISMATCH, Pcigar, pcigar_is_utilised,
};
use crate::wavefront::{BtBlockIdx, Wavefront};

/// k-centered offset and backtrace pointers of one wavefront (used by the
/// offload).
#[derive(Clone, Copy)]
struct Cells {
    off: *mut WfOffset,
    pcigar: *mut Pcigar,
    prev: *mut BtBlockIdx,
}

impl Cells {
    #[inline(always)]
    unsafe fn of(wf: &Wavefront) -> Self {
        let (pcigar, prev) = wf.bt_centered_ptrs();
        Self {
            off: wf.offsets_centered_ptr() as *mut WfOffset,
            pcigar,
            prev,
        }
    }
    #[inline(always)]
    unsafe fn off(&self, k: i32) -> WfOffset {
        *self.off.offset(k as isize)
    }
    #[inline(always)]
    unsafe fn pc(&self, k: i32) -> Pcigar {
        *self.pcigar.offset(k as isize)
    }
    #[inline(always)]
    unsafe fn pv(&self, k: i32) -> BtBlockIdx {
        *self.prev.offset(k as isize)
    }
}

/// `n` cells of an input wavefront starting at diagonal `k`.
struct In<'a> {
    off: &'a [WfOffset],
    pc: &'a [Pcigar],
    pv: &'a [BtBlockIdx],
}

impl<'a> In<'a> {
    #[inline(always)]
    unsafe fn at(wf: &'a Wavefront, k: i32, n: usize) -> Self {
        let c = Cells::of(wf);
        let k = k as isize;
        Self {
            off: std::slice::from_raw_parts(c.off.offset(k), n),
            pc: std::slice::from_raw_parts(c.pcigar.offset(k), n),
            pv: std::slice::from_raw_parts(c.prev.offset(k), n),
        }
    }
}

/// `n` cells of an output wavefront starting at diagonal `k`.
struct Out<'a> {
    off: &'a mut [WfOffset],
    pc: &'a mut [Pcigar],
    pv: &'a mut [BtBlockIdx],
}

impl<'a> Out<'a> {
    #[inline(always)]
    unsafe fn at(wf: &'a mut Wavefront, k: i32, n: usize) -> Self {
        let c = Cells::of(wf);
        let k = k as isize;
        Self {
            off: std::slice::from_raw_parts_mut(c.off.offset(k), n),
            pc: std::slice::from_raw_parts_mut(c.pcigar.offset(k), n),
            pv: std::slice::from_raw_parts_mut(c.prev.offset(k), n),
        }
    }
}

// Branch-free selects throughout (as C's predicated kernels). `push`
// appends a 2-bit operation.
#[inline(always)]
fn push(pcigar: Pcigar, op: Pcigar) -> Pcigar {
    (pcigar << 2) | op
}
#[inline(always)]
fn sel<T>(c: bool, a: T, b: T) -> T {
    if c { a } else { b }
}
/// All ones if `c`, else zero.
#[inline(always)]
fn mask(c: bool) -> u32 {
    (c as u32).wrapping_neg()
}
/// `a` where `m` is all ones, else `b`. Both values are loaded before the
/// choice: a plain `if` over two loads lets LLVM load from a selected
/// address instead, which is a gather and stops vectorization.
#[inline(always)]
fn blend(m: u32, a: u32, b: u32) -> u32 {
    (a & m) | (b & !m)
}
const DEL: Pcigar = PCIGAR_DELETION as Pcigar;
const MISMS: Pcigar = PCIGAR_MISMATCH as Pcigar;
const INS: Pcigar = PCIGAR_INSERTION as Pcigar;

/// NULL out an offset beyond the sequences (C does this after the max).
#[inline(always)]
fn bound(k: i32, offset: WfOffset, plen: i32, tlen: i32) -> WfOffset {
    let out = (wavefront_h(k, offset) as u32 > tlen as u32)
        | (wavefront_v(k, offset) as u32 > plen as u32);
    sel(out, OFFSET_NULL, offset)
}

/// Number of cells in `[lo, hi]`.
#[inline(always)]
fn len(lo: i32, hi: i32) -> usize {
    (hi - lo + 1).max(0) as usize
}

/// Reslice to `n` cells, so the kernel loops have no bounds checks.
macro_rules! cut {
    ($n:expr; $($s:ident),+) => { $(let $s = &$s[..$n];)+ };
}
macro_rules! cut_mut {
    ($n:expr; $($s:ident),+) => { $(let $s = &mut $s[..$n];)+ };
}

/// One gap cell (C's I or D update): from the M-open input or the
/// extension input, the extension winning ties.
#[inline(always)]
fn gap(
    open: WfOffset,
    open_pc: Pcigar,
    open_pv: BtBlockIdx,
    ext: WfOffset,
    ext_pc: Pcigar,
    ext_pv: BtBlockIdx,
) -> (WfOffset, Pcigar, BtBlockIdx) {
    let m = mask(ext >= open);
    (
        blend(m, ext as u32, open as u32) as WfOffset,
        blend(m, ext_pc, open_pc),
        blend(m, ext_pv, open_pv),
    )
}

// The kernels below take every array as a separate slice parameter, so the
// compiler knows outputs don't alias inputs and vectorizes them (C gets the
// same with `PRAGMA_LOOP_VECTORIZE`). They are not inlined into their
// callers, which would lose that.

/// Edit distance (C's `wavefront_compute_edit_idm_piggyback`).
///
/// # Safety
/// `prev` must cover `[lo - 1, hi + 1]` and `curr` `[lo, hi]`, both with
/// backtrace arrays.
pub unsafe fn compute_edit_idm(
    plen: i32,
    tlen: i32,
    prev: &Wavefront,
    curr: &mut Wavefront,
    lo: i32,
    hi: i32,
) {
    let n = len(lo, hi);
    let (pi, px, pd) = (
        In::at(prev, lo - 1, n),
        In::at(prev, lo, n),
        In::at(prev, lo + 1, n),
    );
    let o = Out::at(curr, lo, n);
    edit_kernel(
        lo, plen, tlen, pi.off, pi.pc, pi.pv, px.off, px.pc, px.pv, pd.off, pd.pc, pd.pv, o.off,
        o.pc, o.pv,
    );
}

#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn edit_kernel(
    lo: i32,
    plen: i32,
    tlen: i32,
    i_off: &[WfOffset],
    i_pc: &[Pcigar],
    i_pv: &[BtBlockIdx],
    x_off: &[WfOffset],
    x_pc: &[Pcigar],
    x_pv: &[BtBlockIdx],
    d_off: &[WfOffset],
    d_pc: &[Pcigar],
    d_pv: &[BtBlockIdx],
    o_off: &mut [WfOffset],
    o_pc: &mut [Pcigar],
    o_pv: &mut [BtBlockIdx],
) {
    let n = o_off.len();
    cut!(n; i_off, i_pc, i_pv, x_off, x_pc, x_pv, d_off, d_pc, d_pv);
    cut_mut!(n; o_pc, o_pv);
    for i in 0..n {
        let ins = i_off[i] + 1;
        let del = d_off[i];
        let misms = x_off[i] + 1;
        let max = del.max(ins.max(misms));
        // Later checks win: ins, then del, then misms.
        let (d, x) = (mask(max == del), mask(max == misms));
        let pc = blend(d, push(d_pc[i], DEL), push(i_pc[i], INS));
        o_pc[i] = blend(x, push(x_pc[i], MISMS), pc);
        o_pv[i] = blend(x, x_pv[i], blend(d, d_pv[i], i_pv[i]));
        o_off[i] = bound(lo + i as i32, max, plen, tlen);
    }
}

/// Indel distance (C's `wavefront_compute_indel_idm_piggyback`).
///
/// # Safety
/// As [`compute_edit_idm`].
pub unsafe fn compute_indel_idm(
    plen: i32,
    tlen: i32,
    prev: &Wavefront,
    curr: &mut Wavefront,
    lo: i32,
    hi: i32,
) {
    let n = len(lo, hi);
    let (pi, pd) = (In::at(prev, lo - 1, n), In::at(prev, lo + 1, n));
    let o = Out::at(curr, lo, n);
    indel_kernel(
        lo, plen, tlen, pi.off, pi.pc, pi.pv, pd.off, pd.pc, pd.pv, o.off, o.pc, o.pv,
    );
}

#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn indel_kernel(
    lo: i32,
    plen: i32,
    tlen: i32,
    i_off: &[WfOffset],
    i_pc: &[Pcigar],
    i_pv: &[BtBlockIdx],
    d_off: &[WfOffset],
    d_pc: &[Pcigar],
    d_pv: &[BtBlockIdx],
    o_off: &mut [WfOffset],
    o_pc: &mut [Pcigar],
    o_pv: &mut [BtBlockIdx],
) {
    let n = o_off.len();
    cut!(n; i_off, i_pc, i_pv, d_off, d_pc, d_pv);
    cut_mut!(n; o_pc, o_pv);
    for i in 0..n {
        let ins = i_off[i] + 1;
        let del = d_off[i];
        let max = del.max(ins);
        let d = mask(max == del);
        o_pc[i] = blend(d, push(d_pc[i], DEL), push(i_pc[i], INS));
        o_pv[i] = blend(d, d_pv[i], i_pv[i]);
        o_off[i] = bound(lo + i as i32, max, plen, tlen);
    }
}

/// Gap-linear (C's `wavefront_compute_linear_idm_piggyback`).
///
/// # Safety
/// Inputs must cover `[lo - 1, hi + 1]` and the output `[lo, hi]`, all with
/// backtrace arrays; the output must not be an input.
pub unsafe fn compute_linear_idm(
    plen: i32,
    tlen: i32,
    m_misms: &Wavefront,
    m_open: &Wavefront,
    m_curr: &mut Wavefront,
    lo: i32,
    hi: i32,
) {
    let n = len(lo, hi);
    let mx = In::at(m_misms, lo, n);
    let (oi, od) = (In::at(m_open, lo - 1, n), In::at(m_open, lo + 1, n));
    let o = Out::at(m_curr, lo, n);
    // Same cell update as edit, with misms and the gaps from their inputs
    edit_kernel(
        lo, plen, tlen, oi.off, oi.pc, oi.pv, mx.off, mx.pc, mx.pv, od.off, od.pc, od.pv, o.off,
        o.pc, o.pv,
    );
}

/// Gap-affine (C's `wavefront_compute_affine_idm_piggyback`): I1, D1 and M
/// in one pass, as C's fused loop.
///
/// # Safety
/// As [`compute_linear_idm`].
#[allow(clippy::too_many_arguments)]
pub unsafe fn compute_affine_idm(
    plen: i32,
    tlen: i32,
    m_misms: &Wavefront,
    m_open: &Wavefront,
    i1_ext: &Wavefront,
    d1_ext: &Wavefront,
    m_curr: &mut Wavefront,
    i1_curr: &mut Wavefront,
    d1_curr: &mut Wavefront,
    lo: i32,
    hi: i32,
) {
    let n = len(lo, hi);
    let mx = In::at(m_misms, lo, n);
    let (io, ie) = (In::at(m_open, lo - 1, n), In::at(i1_ext, lo - 1, n));
    let (dop, de) = (In::at(m_open, lo + 1, n), In::at(d1_ext, lo + 1, n));
    let (m, i1, d1) = (
        Out::at(m_curr, lo, n),
        Out::at(i1_curr, lo, n),
        Out::at(d1_curr, lo, n),
    );
    affine_kernel(
        lo, plen, tlen, mx.off, mx.pc, mx.pv, io.off, io.pc, io.pv, ie.off, ie.pc, ie.pv, dop.off,
        dop.pc, dop.pv, de.off, de.pc, de.pv, m.off, m.pc, m.pv, i1.off, i1.pc, i1.pv, d1.off,
        d1.pc, d1.pv,
    );
}

#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn affine_kernel(
    lo: i32,
    plen: i32,
    tlen: i32,
    mx_off: &[WfOffset],
    mx_pc: &[Pcigar],
    mx_pv: &[BtBlockIdx],
    io_off: &[WfOffset],
    io_pc: &[Pcigar],
    io_pv: &[BtBlockIdx],
    ie_off: &[WfOffset],
    ie_pc: &[Pcigar],
    ie_pv: &[BtBlockIdx],
    do_off: &[WfOffset],
    do_pc: &[Pcigar],
    do_pv: &[BtBlockIdx],
    de_off: &[WfOffset],
    de_pc: &[Pcigar],
    de_pv: &[BtBlockIdx],
    m_off: &mut [WfOffset],
    m_pc: &mut [Pcigar],
    m_pv: &mut [BtBlockIdx],
    i_off: &mut [WfOffset],
    i_pc: &mut [Pcigar],
    i_pv: &mut [BtBlockIdx],
    d_off: &mut [WfOffset],
    d_pc: &mut [Pcigar],
    d_pv: &mut [BtBlockIdx],
) {
    let n = m_off.len();
    cut!(n; mx_off, mx_pc, mx_pv, io_off, io_pc, io_pv, ie_off, ie_pc, ie_pv);
    cut!(n; do_off, do_pc, do_pv, de_off, de_pc, de_pv);
    cut_mut!(n; m_pc, m_pv, i_off, i_pc, i_pv, d_off, d_pc, d_pv);
    for i in 0..n {
        let (ins, ins_pc, ins_pv) =
            gap(io_off[i], io_pc[i], io_pv[i], ie_off[i], ie_pc[i], ie_pv[i]);
        let (ins, ins_pc) = (ins + 1, push(ins_pc, INS));
        let (del, del_pc, del_pv) =
            gap(do_off[i], do_pc[i], do_pv[i], de_off[i], de_pc[i], de_pv[i]);
        let del_pc = push(del_pc, DEL);
        (i_off[i], i_pc[i], i_pv[i]) = (ins, ins_pc, ins_pv);
        (d_off[i], d_pc[i], d_pv[i]) = (del, del_pc, del_pv);
        // M: I1, then D1 (later wins ties), then misms. From I/D the X is a
        // fake gap-close; from M it is a real mismatch.
        let misms = mx_off[i] + 1;
        let max = del.max(misms.max(ins));
        let (d, x) = (mask(max == del), mask(max == misms));
        let pc = blend(d, del_pc, ins_pc);
        m_pc[i] = push(blend(x, mx_pc[i], pc), MISMS);
        m_pv[i] = blend(x, mx_pv[i], blend(d, del_pv, ins_pv));
        m_off[i] = bound(lo + i as i32, max, plen, tlen);
    }
}

/// Gap-affine 2-piece (C's `wavefront_compute_affine2p_idm_piggyback`): I1,
/// I2, D1, D2 and M in one pass, as C's fused loop.
///
/// # Safety
/// As [`compute_linear_idm`].
#[allow(clippy::too_many_arguments)]
pub unsafe fn compute_affine2p_idm(
    plen: i32,
    tlen: i32,
    m_misms: &Wavefront,
    m_open1: &Wavefront,
    m_open2: &Wavefront,
    i1_ext: &Wavefront,
    i2_ext: &Wavefront,
    d1_ext: &Wavefront,
    d2_ext: &Wavefront,
    m_curr: &mut Wavefront,
    i1_curr: &mut Wavefront,
    i2_curr: &mut Wavefront,
    d1_curr: &mut Wavefront,
    d2_curr: &mut Wavefront,
    lo: i32,
    hi: i32,
) {
    let n = len(lo, hi);
    let mx = In::at(m_misms, lo, n);
    let (io1, ie1) = (In::at(m_open1, lo - 1, n), In::at(i1_ext, lo - 1, n));
    let (io2, ie2) = (In::at(m_open2, lo - 1, n), In::at(i2_ext, lo - 1, n));
    let (do1, de1) = (In::at(m_open1, lo + 1, n), In::at(d1_ext, lo + 1, n));
    let (do2, de2) = (In::at(m_open2, lo + 1, n), In::at(d2_ext, lo + 1, n));
    let m = Out::at(m_curr, lo, n);
    let (i1, i2) = (Out::at(i1_curr, lo, n), Out::at(i2_curr, lo, n));
    let (d1, d2) = (Out::at(d1_curr, lo, n), Out::at(d2_curr, lo, n));
    affine2p_kernel(
        lo, plen, tlen, mx.off, mx.pc, mx.pv, io1.off, io1.pc, io1.pv, ie1.off, ie1.pc, ie1.pv,
        io2.off, io2.pc, io2.pv, ie2.off, ie2.pc, ie2.pv, do1.off, do1.pc, do1.pv, de1.off, de1.pc,
        de1.pv, do2.off, do2.pc, do2.pv, de2.off, de2.pc, de2.pv, m.off, m.pc, m.pv, i1.off, i1.pc,
        i1.pv, i2.off, i2.pc, i2.pv, d1.off, d1.pc, d1.pv, d2.off, d2.pc, d2.pv,
    );
}

#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn affine2p_kernel(
    lo: i32,
    plen: i32,
    tlen: i32,
    mx_off: &[WfOffset],
    mx_pc: &[Pcigar],
    mx_pv: &[BtBlockIdx],
    io1_off: &[WfOffset],
    io1_pc: &[Pcigar],
    io1_pv: &[BtBlockIdx],
    ie1_off: &[WfOffset],
    ie1_pc: &[Pcigar],
    ie1_pv: &[BtBlockIdx],
    io2_off: &[WfOffset],
    io2_pc: &[Pcigar],
    io2_pv: &[BtBlockIdx],
    ie2_off: &[WfOffset],
    ie2_pc: &[Pcigar],
    ie2_pv: &[BtBlockIdx],
    do1_off: &[WfOffset],
    do1_pc: &[Pcigar],
    do1_pv: &[BtBlockIdx],
    de1_off: &[WfOffset],
    de1_pc: &[Pcigar],
    de1_pv: &[BtBlockIdx],
    do2_off: &[WfOffset],
    do2_pc: &[Pcigar],
    do2_pv: &[BtBlockIdx],
    de2_off: &[WfOffset],
    de2_pc: &[Pcigar],
    de2_pv: &[BtBlockIdx],
    m_off: &mut [WfOffset],
    m_pc: &mut [Pcigar],
    m_pv: &mut [BtBlockIdx],
    i1_off: &mut [WfOffset],
    i1_pc: &mut [Pcigar],
    i1_pv: &mut [BtBlockIdx],
    i2_off: &mut [WfOffset],
    i2_pc: &mut [Pcigar],
    i2_pv: &mut [BtBlockIdx],
    d1_off: &mut [WfOffset],
    d1_pc: &mut [Pcigar],
    d1_pv: &mut [BtBlockIdx],
    d2_off: &mut [WfOffset],
    d2_pc: &mut [Pcigar],
    d2_pv: &mut [BtBlockIdx],
) {
    let n = m_off.len();
    cut!(n; mx_off, mx_pc, mx_pv, io1_off, io1_pc, io1_pv, ie1_off, ie1_pc, ie1_pv);
    cut!(n; io2_off, io2_pc, io2_pv, ie2_off, ie2_pc, ie2_pv);
    cut!(n; do1_off, do1_pc, do1_pv, de1_off, de1_pc, de1_pv);
    cut!(n; do2_off, do2_pc, do2_pv, de2_off, de2_pc, de2_pv);
    cut_mut!(n; m_pc, m_pv, i1_off, i1_pc, i1_pv, i2_off, i2_pc, i2_pv);
    cut_mut!(n; d1_off, d1_pc, d1_pv, d2_off, d2_pc, d2_pv);
    for i in 0..n {
        let (ins1, ins1_pc, ins1_pv) = gap(
            io1_off[i], io1_pc[i], io1_pv[i], ie1_off[i], ie1_pc[i], ie1_pv[i],
        );
        let (ins1, ins1_pc) = (ins1 + 1, push(ins1_pc, INS));
        let (ins2, ins2_pc, ins2_pv) = gap(
            io2_off[i], io2_pc[i], io2_pv[i], ie2_off[i], ie2_pc[i], ie2_pv[i],
        );
        let (ins2, ins2_pc) = (ins2 + 1, push(ins2_pc, INS));
        let (del1, del1_pc, del1_pv) = gap(
            do1_off[i], do1_pc[i], do1_pv[i], de1_off[i], de1_pc[i], de1_pv[i],
        );
        let del1_pc = push(del1_pc, DEL);
        let (del2, del2_pc, del2_pv) = gap(
            do2_off[i], do2_pc[i], do2_pv[i], de2_off[i], de2_pc[i], de2_pv[i],
        );
        let del2_pc = push(del2_pc, DEL);
        (i1_off[i], i1_pc[i], i1_pv[i]) = (ins1, ins1_pc, ins1_pv);
        (i2_off[i], i2_pc[i], i2_pv[i]) = (ins2, ins2_pc, ins2_pv);
        (d1_off[i], d1_pc[i], d1_pv[i]) = (del1, del1_pc, del1_pv);
        (d2_off[i], d2_pc[i], d2_pv[i]) = (del2, del2_pc, del2_pv);
        // M: I1, I2, D1, D2 (later wins ties), then misms
        let misms = mx_off[i] + 1;
        let max = ins1.max(ins2).max(del1).max(del2).max(misms);
        let (s2, t1, t2, x) = (
            mask(max == ins2),
            mask(max == del1),
            mask(max == del2),
            mask(max == misms),
        );
        let pc = blend(t2, del2_pc, blend(t1, del1_pc, blend(s2, ins2_pc, ins1_pc)));
        let pv = blend(t2, del2_pv, blend(t1, del1_pv, blend(s2, ins2_pv, ins1_pv)));
        m_pc[i] = push(blend(x, mx_pc[i], pc), MISMS);
        m_pv[i] = blend(x, mx_pv[i], pv);
        m_off[i] = bound(lo + i as i32, max, plen, tlen);
    }
}

// --- Offloading (C's wavefront_backtrace_offload.c) ---

/// C's `wavefront_backtrace_offload_blocks_selective`: move every pcigar of
/// a live cell (offset >= 0) at or above `mask` into a buffer block, and
/// restart the cell from that block.
fn offload_selective(
    wf: &mut Wavefront,
    lo: i32,
    hi: i32,
    mask: Pcigar,
    buf: &mut BacktraceBuffer,
) {
    // SAFETY: [lo, hi] is within the output wavefront's allocation.
    unsafe {
        let c = Cells::of(wf);
        for k in lo..=hi {
            let (pcigar, prev) = (c.pc(k), c.pv(k));
            if c.off(k) >= 0 && pcigar_is_utilised(pcigar, mask) {
                let idx = buf.store_block(pcigar, prev);
                *c.pcigar.offset(k as isize) = 0;
                *c.prev.offset(k as isize) = idx;
            }
        }
    }
}

/// Offload for edit/indel/linear (C's
/// `wavefront_backtrace_offload_blocks_linear`); returns the occupancy left.
pub fn offload_blocks_linear(
    wf: &mut Wavefront,
    lo: i32,
    hi: i32,
    low: bool,
    buf: &mut BacktraceBuffer,
) -> i32 {
    if low {
        offload_selective(wf, lo, hi, PCIGAR_FULL_MASK, buf);
        PCIGAR_MAX_LENGTH as i32 - 1 // At least 1 slot free
    } else {
        offload_selective(wf, lo, hi, PCIGAR_HALF_FULL_MASK, buf);
        PCIGAR_MAX_LENGTH as i32 / 2 // Half occupancy
    }
}

/// Offload for affine/affine2p (C's
/// `wavefront_backtrace_offload_blocks_affine`); returns the occupancy left.
fn offload_blocks_affine(
    wf: &mut Wavefront,
    lo: i32,
    hi: i32,
    low: bool,
    buf: &mut BacktraceBuffer,
) -> i32 {
    if low {
        offload_selective(wf, lo, hi, PCIGAR_ALMOST_FULL_MASK, buf);
        PCIGAR_MAX_LENGTH as i32 - 2 // At least 2 slots free
    } else {
        offload_selective(wf, lo, hi, PCIGAR_HALF_FULL_MASK, buf);
        PCIGAR_MAX_LENGTH as i32 / 2 // Half occupancy
    }
}

/// Occupancy contributed by an input wavefront (0 if null).
#[inline(always)]
fn occ(wf: &Wavefront) -> i32 {
    if wf.null { 0 } else { wf.bt_occupancy_max }
}

/// C's `wavefront_backtrace_offload_linear` (gap-linear).
pub fn offload_linear(
    m_misms: &Wavefront,
    m_open: &Wavefront,
    out_m: &mut Wavefront,
    lo: i32,
    hi: i32,
    low: bool,
    buf: &mut BacktraceBuffer,
) {
    let occ_max = occ(m_open).max(occ(m_misms)) + 1;
    out_m.bt_occupancy_max = occ_max;
    if !out_m.null && occ_max >= PCIGAR_MAX_LENGTH as i32 {
        out_m.bt_occupancy_max = offload_blocks_linear(out_m, lo, hi, low, buf);
    }
}

/// Occupancy of a gap component fed by M-open and its own extension
/// (C: `m_open->bt_occupancy_max + 1`, then the max with `ext + 1`).
#[inline(always)]
fn occ_gap(m_open: &Wavefront, ext: &Wavefront) -> i32 {
    let mut o = if m_open.null {
        0
    } else {
        m_open.bt_occupancy_max + 1
    };
    if !ext.null {
        o = o.max(ext.bt_occupancy_max + 1);
    }
    o
}

/// Offload one affine output component if its occupancy is near full.
#[inline(always)]
fn offload_affine_component(
    wf: &mut Wavefront,
    lo: i32,
    hi: i32,
    low: bool,
    buf: &mut BacktraceBuffer,
) {
    if !wf.null && wf.bt_occupancy_max >= PCIGAR_MAX_LENGTH as i32 - 1 {
        wf.bt_occupancy_max = offload_blocks_affine(wf, lo, hi, low, buf);
    }
}

/// C's `wavefront_backtrace_offload_affine` for gap-affine (M, I1, D1).
#[allow(clippy::too_many_arguments)]
pub fn offload_affine(
    m_misms: &Wavefront,
    m_open: &Wavefront,
    i1_ext: &Wavefront,
    d1_ext: &Wavefront,
    out_m: &mut Wavefront,
    out_i1: &mut Wavefront,
    out_d1: &mut Wavefront,
    lo: i32,
    hi: i32,
    low: bool,
    buf: &mut BacktraceBuffer,
) {
    let occ_i1 = occ_gap(m_open, i1_ext);
    let occ_d1 = occ_gap(m_open, d1_ext);
    let occ_m = occ(m_misms).max(occ_i1).max(occ_d1) + 1;
    out_i1.bt_occupancy_max = occ_i1;
    out_d1.bt_occupancy_max = occ_d1;
    out_m.bt_occupancy_max = occ_m;
    offload_affine_component(out_m, lo, hi, low, buf);
    offload_affine_component(out_i1, lo, hi, low, buf);
    offload_affine_component(out_d1, lo, hi, low, buf);
}

/// C's `wavefront_backtrace_offload_affine` for gap-affine 2-piece.
#[allow(clippy::too_many_arguments)]
pub fn offload_affine2p(
    m_misms: &Wavefront,
    m_open1: &Wavefront,
    m_open2: &Wavefront,
    i1_ext: &Wavefront,
    i2_ext: &Wavefront,
    d1_ext: &Wavefront,
    d2_ext: &Wavefront,
    out_m: &mut Wavefront,
    out_i1: &mut Wavefront,
    out_i2: &mut Wavefront,
    out_d1: &mut Wavefront,
    out_d2: &mut Wavefront,
    lo: i32,
    hi: i32,
    low: bool,
    buf: &mut BacktraceBuffer,
) {
    let occ_i1 = occ_gap(m_open1, i1_ext);
    let occ_d1 = occ_gap(m_open1, d1_ext);
    let occ_i2 = occ_gap(m_open2, i2_ext);
    let occ_d2 = occ_gap(m_open2, d2_ext);
    let occ_m = occ(m_misms).max(occ_i1).max(occ_i2).max(occ_d1).max(occ_d2) + 1;
    out_i1.bt_occupancy_max = occ_i1;
    out_i2.bt_occupancy_max = occ_i2;
    out_d1.bt_occupancy_max = occ_d1;
    out_d2.bt_occupancy_max = occ_d2;
    out_m.bt_occupancy_max = occ_m;
    offload_affine_component(out_m, lo, hi, low, buf);
    offload_affine_component(out_i1, lo, hi, low, buf);
    offload_affine_component(out_d1, lo, hi, low, buf);
    offload_affine_component(out_i2, lo, hi, low, buf);
    offload_affine_component(out_d2, lo, hi, low, buf);
}
