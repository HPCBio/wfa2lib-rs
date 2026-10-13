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

/// `n` cells of an input wavefront starting at diagonal `k`. The kernels
/// work on slices so the compiler knows outputs don't alias inputs and can
/// vectorize (C gets the same with `PRAGMA_LOOP_VECTORIZE`).
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
    for i in 0..n {
        let k = lo + i as i32;
        let ins = pi.off[i] + 1;
        let del = pd.off[i];
        let misms = px.off[i] + 1;
        let max = del.max(ins.max(misms));
        // Later checks win: ins, then del, then misms.
        let (d, x) = (max == del, max == misms);
        let pc = sel(
            x,
            push(px.pc[i], MISMS),
            sel(d, push(pd.pc[i], DEL), push(pi.pc[i], INS)),
        );
        let pv = sel(x, px.pv[i], sel(d, pd.pv[i], pi.pv[i]));
        o.off[i] = bound(k, max, plen, tlen);
        o.pc[i] = pc;
        o.pv[i] = pv;
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
    for i in 0..n {
        let k = lo + i as i32;
        let ins = pi.off[i] + 1;
        let del = pd.off[i];
        let max = del.max(ins);
        let d = max == del;
        o.off[i] = bound(k, max, plen, tlen);
        o.pc[i] = sel(d, push(pd.pc[i], DEL), push(pi.pc[i], INS));
        o.pv[i] = sel(d, pd.pv[i], pi.pv[i]);
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
    for i in 0..n {
        let k = lo + i as i32;
        let ins = oi.off[i] + 1;
        let del = od.off[i];
        let misms = mx.off[i] + 1;
        let max = del.max(misms.max(ins));
        let (d, x) = (max == del, max == misms);
        let pc = sel(
            x,
            push(mx.pc[i], MISMS),
            sel(d, push(od.pc[i], DEL), push(oi.pc[i], INS)),
        );
        let pv = sel(x, mx.pv[i], sel(d, od.pv[i], oi.pv[i]));
        o.off[i] = bound(k, max, plen, tlen);
        o.pc[i] = pc;
        o.pv[i] = pv;
    }
}

/// One gap component (C's I or D update): from the M-open input or the
/// extension input (extension wins ties), plus `inc` (1 for I, 0 for D),
/// pushing `op`. All slices are cells `[lo, hi]` of their component, taken
/// at the right diagonal shift. Slices as parameters tell the compiler
/// outputs don't alias inputs, so this vectorizes.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn gap_pass(
    open_off: &[WfOffset],
    open_pc: &[Pcigar],
    open_pv: &[BtBlockIdx],
    ext_off: &[WfOffset],
    ext_pc: &[Pcigar],
    ext_pv: &[BtBlockIdx],
    out_off: &mut [WfOffset],
    out_pc: &mut [Pcigar],
    out_pv: &mut [BtBlockIdx],
    inc: WfOffset,
    op: Pcigar,
) {
    let n = out_off.len();
    let (open_off, open_pc, open_pv) = (&open_off[..n], &open_pc[..n], &open_pv[..n]);
    let (ext_off, ext_pc, ext_pv) = (&ext_off[..n], &ext_pc[..n], &ext_pv[..n]);
    let (out_pc, out_pv) = (&mut out_pc[..n], &mut out_pv[..n]);
    for i in 0..n {
        let x = ext_off[i] >= open_off[i];
        out_off[i] = sel(x, ext_off[i], open_off[i]) + inc;
        out_pc[i] = push(sel(x, ext_pc[i], open_pc[i]), op);
        out_pv[i] = sel(x, ext_pv[i], open_pv[i]);
    }
}

/// Fold a gap component into the M candidate (C's `if (max == ins/del)`):
/// where it equals `max`, it wins over earlier candidates.
#[inline(never)]
fn m_fold(
    max: &[WfOffset],
    g_off: &[WfOffset],
    g_pc: &[Pcigar],
    g_pv: &[BtBlockIdx],
    pc: &mut [Pcigar],
    pv: &mut [BtBlockIdx],
) {
    let n = max.len();
    let (g_off, g_pc, g_pv) = (&g_off[..n], &g_pc[..n], &g_pv[..n]);
    let (pc, pv) = (&mut pc[..n], &mut pv[..n]);
    for i in 0..n {
        let c = max[i] == g_off[i];
        pc[i] = sel(c, g_pc[i], pc[i]);
        pv[i] = sel(c, g_pv[i], pv[i]);
    }
}

/// M = max(misms, gaps...) with C's tie order (gaps in order, then misms
/// last), the X push and the out-of-bounds NULL. `max` holds the max over
/// the gap components on entry; `pc`/`pv` the winning gap's backtrace.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn m_finish(
    lo: i32,
    plen: i32,
    tlen: i32,
    mx_off: &[WfOffset],
    mx_pc: &[Pcigar],
    mx_pv: &[BtBlockIdx],
    max: &mut [WfOffset],
    pc: &mut [Pcigar],
    pv: &mut [BtBlockIdx],
) {
    let n = max.len();
    let (mx_off, mx_pc, mx_pv) = (&mx_off[..n], &mx_pc[..n], &mx_pv[..n]);
    let (pc, pv) = (&mut pc[..n], &mut pv[..n]);
    for i in 0..n {
        let misms = mx_off[i] + 1;
        let m = misms.max(max[i]);
        let x = m == misms;
        // From I/D the X is a fake gap-close; from M it is a real mismatch.
        pc[i] = push(sel(x, mx_pc[i], pc[i]), MISMS);
        pv[i] = sel(x, mx_pv[i], pv[i]);
        max[i] = bound(lo + i as i32, m, plen, tlen);
    }
}

/// Elementwise max of gap offsets into `max`.
#[inline(never)]
fn max_into(max: &mut [WfOffset], g: &[WfOffset]) {
    let n = max.len();
    let g = &g[..n];
    for i in 0..n {
        max[i] = max[i].max(g[i]);
    }
}

/// Run a gap pass from inputs `open`/`ext` (at diagonal `k_in`) into `out`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn gap_into(
    open: &Wavefront,
    ext: &Wavefront,
    k_in: i32,
    out: &mut Wavefront,
    lo: i32,
    n: usize,
    inc: WfOffset,
    op: Pcigar,
) {
    let (o, e) = (In::at(open, k_in, n), In::at(ext, k_in, n));
    let w = Out::at(out, lo, n);
    gap_pass(
        o.off, o.pc, o.pv, e.off, e.pc, e.pv, w.off, w.pc, w.pv, inc, op,
    );
}

/// Gap-affine (C's `wavefront_compute_affine_idm_piggyback`), computed as
/// separate I, D and M passes over the cells (same result as C's fused loop,
/// whose M update only reads this cell's I/D).
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
    gap_into(m_open, i1_ext, lo - 1, i1_curr, lo, n, 1, INS);
    gap_into(m_open, d1_ext, lo + 1, d1_curr, lo, n, 0, DEL);
    let (i1, d1) = (In::at(i1_curr, lo, n), In::at(d1_curr, lo, n));
    let m = Out::at(m_curr, lo, n);
    // M candidates in C's order: I1, then D1 (later wins ties), then misms.
    m.off.copy_from_slice(i1.off);
    m.pc.copy_from_slice(i1.pc);
    m.pv.copy_from_slice(i1.pv);
    max_into(m.off, d1.off);
    m_fold(m.off, d1.off, d1.pc, d1.pv, m.pc, m.pv);
    let mx = In::at(m_misms, lo, n);
    m_finish(lo, plen, tlen, mx.off, mx.pc, mx.pv, m.off, m.pc, m.pv);
}

/// Gap-affine 2-piece (C's `wavefront_compute_affine2p_idm_piggyback`), as
/// separate passes like [`compute_affine_idm`].
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
    gap_into(m_open1, i1_ext, lo - 1, i1_curr, lo, n, 1, INS);
    gap_into(m_open2, i2_ext, lo - 1, i2_curr, lo, n, 1, INS);
    gap_into(m_open1, d1_ext, lo + 1, d1_curr, lo, n, 0, DEL);
    gap_into(m_open2, d2_ext, lo + 1, d2_curr, lo, n, 0, DEL);
    let m = Out::at(m_curr, lo, n);
    // Max over all gaps first; then C's tie order I1, I2, D1, D2 (later
    // wins), then misms.
    let i1 = In::at(i1_curr, lo, n);
    m.off.copy_from_slice(i1.off);
    m.pc.copy_from_slice(i1.pc);
    m.pv.copy_from_slice(i1.pv);
    let gaps = [
        In::at(i2_curr, lo, n),
        In::at(d1_curr, lo, n),
        In::at(d2_curr, lo, n),
    ];
    for g in &gaps {
        max_into(m.off, g.off);
    }
    // `m.off` is now the max over gaps; I1 wins only where it equals it.
    for g in &gaps {
        m_fold(m.off, g.off, g.pc, g.pv, m.pc, m.pv);
    }
    let mx = In::at(m_misms, lo, n);
    m_finish(lo, plen, tlen, mx.off, mx.pc, mx.pv, m.off, m.pc, m.pv);
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
