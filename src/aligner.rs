//! WavefrontAligner: main state machine and public API for wavefront alignment.

use crate::backtrace;
use crate::bialign::{
    self, BIALIGN_FALLBACK_MIN_LENGTH, BIALIGN_FALLBACK_MIN_SCORE, BIALIGN_RECOVERY_MIN_SCORE,
    BiAlignBreakpoint, ComponentType,
};
use crate::bt_buffer::{BT_BLOCK_IDX_NULL, BacktraceBuffer};
use crate::cigar::Cigar;
use crate::components::WavefrontComponents;
use crate::compute;
use crate::compute::{affine, affine2p, edit, linear, piggyback};
use crate::extend;
use crate::heuristic::{self, HeuristicState, HeuristicStrategy};
use crate::offset::{DPMATRIX_DIAGONAL_NULL, OFFSET_NULL, WfOffset, wavefront_h, wavefront_v};
use crate::pcigar::{PCIGAR_MAX_LENGTH, Pcigar};
use crate::penalties::{DistanceMetric, WavefrontPenalties};
use crate::sequences::{SequenceMode, WavefrontSequences};
use crate::components::WF_PTR_NONE;
use crate::slab::{SlabMode, WavefrontSlab};
use crate::termination;
use crate::wavefront::{BtBlockIdx, WavefrontPos, WavefrontStatus};

/// Alignment computation scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentScope {
    /// Only compute the alignment score.
    ComputeScore,
    /// Compute score and full alignment (CIGAR).
    ComputeAlignment,
}

/// Wavefront memory mode (C's `wavefront_memory_t`, without `ultralow`,
/// which is BiWFA: [`WavefrontAligner::align_biwfa`]).
///
/// Only affects alignments (CIGAR); score-only always keeps just the
/// wavefronts it still needs. `Med` and `Low` keep only those too, and record
/// the alignment as packed CIGARs alongside the wavefronts ("piggyback"
/// backtrace), stored in a backtrace buffer; they differ in how full a packed
/// CIGAR gets before it is offloaded. Between co-optimal alignments they may
/// choose differently from `High`, exactly as C does.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MemoryMode {
    /// Keep every wavefront and backtrace through them (the default).
    #[default]
    High,
    /// Piggyback backtrace, offloading packed CIGARs at half occupancy.
    Med,
    /// Piggyback backtrace, offloading packed CIGARs only when (nearly) full.
    Low,
}

/// Whether to reproduce known C WFA2-lib bugs (HPCBio/wfa2lib-rs#32).
///
/// By default the aligner's output is identical to C WFA2-lib, including
/// results that C gets wrong. [`Compat::Corrected`] fixes these known
/// issues, so output can differ from C in exactly those cases:
///
/// - **Ends-free with a match reward** (`match < 0`): C ends the alignment
///   at the first cell to reach a free end, but with a match reward the
///   score of an ends-free alignment depends on how many bases its free end
///   gaps skip, so that cell need not be the best. The result can be
///   suboptimal, from a trailing-end miss of a point or two
///   (Marco-Sola-Lab/WFA2-lib#102) to a shift by a repeat unit (HPCBio/wfa2lib-rs#37)
///   or the all-gap, no-overlap alignment when the free ends span whole
///   sequences (HPCBio/wfa2lib-rs#38). Corrected mode charges skipped
///   trailing bases like the begin-free ones and returns the optimal
///   alignment. Its score, as returned by
///   [`WavefrontAligner::align_endsfree`], is then the WFA score of the whole
///   alignment, free end gaps included, so
///   [`WavefrontPenalties::wf_score_to_sw_score`] with the full sequence
///   lengths gives its real score. Also applies to extension alignment.
///   It can take more steps than C (up to the optimal score).
///
/// - **Edit/indel ends-free with begin gaps** (HPCBio/wfa2lib-rs#18, a C
///   regression from `380eb31b`, reported as Marco-Sola-Lab/WFA2-lib#142): C can
///   return an invalid CIGAR (a mismatching column as `M`; the score is
///   right). Corrected mode marks the begin-free seeds as initialized, the
///   fix proposed upstream, so the backtrace can reach them.
///
/// Gap-linear/affine/affine2p with `match == 0` give identical results in
/// both modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Compat {
    /// Output identical to C WFA2-lib, including its known bugs (the default).
    #[default]
    C,
    /// Fix the known C bugs listed in [`Compat`].
    Corrected,
}

/// Corrected-mode ends-free end candidate: the best cell reaching a free
/// end so far, and its total score (see [`Compat::Corrected`]).
#[derive(Debug, Clone, Copy)]
struct EndsfreeCandidate {
    /// Score at the cell plus `-match` per skipped trailing base.
    total: i32,
    /// The cell (score, k, offset).
    pos: WavefrontPos,
    /// The cell's piggyback backtrace `(pcigar, prev)` (med/low), kept from
    /// compaction while it is the candidate.
    bt: (Pcigar, BtBlockIdx),
}

/// Alignment span mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentSpan {
    /// End-to-end (global) alignment.
    End2End,
    /// Ends-free (semi-global) alignment.
    EndsFree,
}

/// Alignment result status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignStatus {
    /// Initial state / in progress.
    Ok,
    /// Alignment target reached (internal, before finalization).
    EndReached,
    /// Target unreachable (heuristic cutoff).
    EndUnreachable,
    /// Maximum alignment steps exceeded.
    MaxStepsReached,
    /// Out of memory.
    Oom,
    /// Alignment completed successfully.
    Completed,
    /// Partial alignment (extension or heuristic).
    Partial,
}

/// Main wavefront aligner.
///
/// The const generic `N` is the number of component slots per score:
/// - N=1: Edit, Indel, GapLinear (only M wavefront)
/// - N=3: GapAffine (M + I1 + D1)
/// - N=5: GapAffine2p (M + I1 + D1 + I2 + D2)
///
/// Use the type aliases [`EditAligner`], [`AffineAligner`], [`Affine2pAligner`]
/// for convenient construction.
pub struct WavefrontAligner<const N: usize> {
    // Sequences
    pub sequences: WavefrontSequences,

    // Wavefront data
    pub wf_components: WavefrontComponents<N>,
    pub wavefront_slab: WavefrontSlab,

    // Configuration
    pub penalties: WavefrontPenalties,
    pub alignment_scope: AlignmentScope,
    pub alignment_span: AlignmentSpan,

    // Alignment state
    pub cigar: Cigar,
    pub alignment_end_pos: WavefrontPos,
    pub status: AlignStatus,
    pub num_null_steps: i32,

    // Ends-free configuration
    pub pattern_begin_free: i32,
    pub pattern_end_free: i32,
    pub text_begin_free: i32,
    pub text_end_free: i32,
    pub extension: bool,

    /// C compatibility or corrected behaviour (see [`Compat`]).
    pub compat: Compat,
    /// Corrected-mode ends-free end candidate (see [`Compat::Corrected`]).
    endsfree_best: Option<EndsfreeCandidate>,

    // Heuristic
    pub heuristic: HeuristicStrategy,
    pub heuristic_state: HeuristicState,

    /// Memory mode for alignments (see [`MemoryMode`]).
    pub memory_mode: MemoryMode,
    /// Backtrace-buffer size (bytes) above which med/low compact it.
    pub max_memory_compact: u64,

    // System limits
    pub max_alignment_steps: i32,

    // BiWFA sub-aligners (lazily initialized, reused across align_biwfa calls)
    biwfa_fwd: Option<Box<WavefrontAligner<N>>>,
    biwfa_rev: Option<Box<WavefrontAligner<N>>>,
    biwfa_base: Option<Box<WavefrontAligner<N>>>,
}

/// Type alias for edit/indel/gap-linear distance (only M wavefront needed).
pub type EditAligner = WavefrontAligner<1>;
/// Type alias for gap-affine distance (M + I1 + D1 wavefronts).
pub type AffineAligner = WavefrontAligner<3>;
/// Type alias for gap-affine 2-piece distance (M + I1 + D1 + I2 + D2 wavefronts).
pub type Affine2pAligner = WavefrontAligner<5>;

/// Start a backtrace chain at cell `k` of `wf`, which is (v, h): an empty
/// pcigar pointing to a new initial block (C's
/// `wf_backtrace_buffer_init_block`).
fn set_bt_init(
    wf: &mut crate::wavefront::Wavefront,
    buf: &mut BacktraceBuffer,
    k: i32,
    v: i32,
    h: i32,
) {
    wf.set_bt_pcigar(k, 0);
    wf.set_bt_prev(k, buf.init_block(v, h));
}

/// C's `wf_backtrace_buffer_unpack_cigar_{linear,affine}`: build the CIGAR
/// from packed blocks (`packed` is end-first), starting at `begin` (preceded
/// by its begin-free I/D) and ending at `end`, then padded with matches and
/// end I/D to the full sequences. With `affine`, matches are only
/// re-extended in M, and an X after I/D (gap close) is dropped.
fn unpack_piggyback(
    packed: &[crate::pcigar::Pcigar],
    sequences: &WavefrontSequences,
    begin: (i32, i32),
    end: (i32, i32),
    affine: bool,
    cigar: &mut Cigar,
) {
    use crate::pcigar::{PCIGAR_FULL_MASK, pcigar_extract, pcigar_free_slots};
    let (plen, tlen) = (sequences.pattern_length, sequences.text_length);
    let lambda = sequences.mode == SequenceMode::Lambda;
    let matches = |v: i32, h: i32| -> bool {
        if lambda {
            sequences.cmp_lambda(v, h)
        } else {
            v < plen && h < tlen && sequences.cmp(v, h)
        }
    };
    let ops = &mut cigar.operations;
    let mut n = 0usize;
    let mut put = |op: u8| {
        ops[n] = op;
        n += 1;
    };
    let (mut v, mut h) = begin;
    (0..h).for_each(|_| put(b'I'));
    (0..v).for_each(|_| put(b'D'));
    let mut in_m = true;
    for &block in packed.iter().rev() {
        let mut pcigar = block;
        let mut length = PCIGAR_MAX_LENGTH;
        if pcigar < PCIGAR_FULL_MASK {
            let free = pcigar_free_slots(pcigar);
            length -= free;
            pcigar = pcigar.checked_shl(2 * free as u32).unwrap_or(0);
        }
        for _ in 0..length {
            if in_m || !affine {
                while matches(v, h) {
                    put(b'M');
                    v += 1;
                    h += 1;
                }
            }
            let op = pcigar_extract(pcigar);
            pcigar <<= 2;
            let (ch, dv, dh, m) = match op {
                1 => (b'D', 1, 0, false),
                2 => (b'X', 1, 1, true),
                3 => (b'I', 0, 1, false),
                _ => (b'?', 0, 0, true),
            };
            if affine && !in_m && ch == b'X' {
                in_m = true;
                continue;
            }
            put(ch);
            v += dv;
            h += dh;
            in_m = m;
        }
    }
    // As C, a negative count (end behind the unpacked position) moves back.
    let num_matches = (end.0 - v).min(end.1 - h);
    for _ in 0..num_matches {
        put(b'M');
    }
    v += num_matches;
    h += num_matches;
    while h < tlen {
        put(b'I');
        h += 1;
    }
    while v < plen {
        put(b'D');
        v += 1;
    }
    cigar.begin_offset = 0;
    cigar.end_offset = n;
}

/// Outcome of a BiWFA breakpoint search or sub-alignment (C's `WF_STATUS_*`
/// values inside `wavefront_bialign`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BiStatus {
    Ok,
    /// One direction reached the end before the wavefronts overlapped (score).
    EndReached(i32),
    Unreachable,
    MaxSteps,
    /// A base-case sub-alignment did not complete.
    Unattainable,
}

/// The heuristic's band, for banded-static / banded-adaptive.
fn heuristic_band(h: &HeuristicStrategy) -> Option<(i32, i32)> {
    match *h {
        HeuristicStrategy::BandedStatic { min_k, max_k }
        | HeuristicStrategy::BandedAdaptive { min_k, max_k, .. } => Some((min_k, max_k)),
        _ => None,
    }
}

/// `h` with its band replaced by [min_k, max_k] (no-op for other heuristics).
fn with_band(h: &HeuristicStrategy, min_k: i32, max_k: i32) -> HeuristicStrategy {
    match *h {
        HeuristicStrategy::BandedStatic { .. } => HeuristicStrategy::BandedStatic { min_k, max_k },
        HeuristicStrategy::BandedAdaptive {
            steps_between_cutoffs,
            ..
        } => HeuristicStrategy::BandedAdaptive {
            min_k,
            max_k,
            steps_between_cutoffs,
        },
        ref other => other.clone(),
    }
}

/// Initial sequence lengths the components are sized for, as C's
/// `wavefront_aligner_new` (`PATTERN_LENGTH_INIT`/`TEXT_LENGTH_INIT`). Arrays
/// only grow afterwards, so their capacity (the score limit in full-memory
/// mode) is at least that for 1000 x 1000 -- heuristics can push an
/// alignment past the bound computed for the actual lengths.
/// C's default `max_memory_compact` (`BUFFER_SIZE_256M`).
const MAX_MEMORY_COMPACT: u64 = 256 << 20;
/// C's `probe_interval_compact` (6000; checked within the 3000 global
/// interval, so every 6000 scores).
const PROBE_INTERVAL_COMPACT: i32 = 6000;
const PATTERN_LENGTH_INIT: i32 = 1000;
const TEXT_LENGTH_INIT: i32 = 1000;

/// CIGAR buffer size for a `plen` × `tlen` alignment. As C, twice the
/// sequence lengths (partial Z-drop backtraces can exceed `plen + tlen`
/// before max-trimming), plus one slot: the backtraces write backwards and
/// leave `begin_offset` one below the first op, which would underflow for
/// an empty buffer (empty pattern and text, HPCBio/wfa2lib-rs#41) or a full one.
#[inline]
fn cigar_buffer_len(plen: i32, tlen: i32) -> usize {
    (2 * (plen + tlen)) as usize + 1
}

/// True if a component slot holds no wavefront or one trimmed to empty
/// (C's `wavefront->null`).
#[inline(always)]
fn is_null_wf(ptr: *mut crate::wavefront::Wavefront) -> bool {
    ptr.is_null() || unsafe { (*ptr).null }
}

/// `ptr` if it holds a non-empty wavefront, else `WF_PTR_NONE` (read as the
/// null wavefront, as C's `wavefront_compute_get_*wavefront` does).
#[inline(always)]
fn live_or_none(ptr: *mut crate::wavefront::Wavefront) -> *mut crate::wavefront::Wavefront {
    if is_null_wf(ptr) { WF_PTR_NONE } else { ptr }
}

/// lo/hi of an input wavefront; a missing one has the null wavefront's
/// bounds (lo = 1, hi = -1), which C includes in the input range union.
#[inline(always)]
fn null_wf_bounds(ptr: *mut crate::wavefront::Wavefront) -> (i32, i32) {
    if ptr.is_null() {
        (1, -1)
    } else {
        unsafe { ((*ptr).lo, (*ptr).hi) }
    }
}

impl<const N: usize> WavefrontAligner<N> {
    /// Create a new aligner with the given penalties.
    pub fn new(penalties: WavefrontPenalties) -> Self {
        let wf_components = WavefrontComponents::<N>::new(
            PATTERN_LENGTH_INIT,
            TEXT_LENGTH_INIT,
            &penalties,
            false,
            false,
        );
        // C's `wavefront_slab_new(1000, ...)`.
        let init_wf_length = 1000;

        Self {
            sequences: WavefrontSequences::new(),
            wf_components,
            wavefront_slab: {
                let mut slab = WavefrontSlab::new(init_wf_length, false, SlabMode::Reuse);
                // Enable cache coloring for 5-component (affine2p) aligners
                // to reduce L1 cache set conflicts between the 12+ input arrays.
                if N == 5 {
                    slab.set_cache_coloring(24);
                }
                slab
            },
            penalties,
            alignment_scope: AlignmentScope::ComputeScore,
            alignment_span: AlignmentSpan::End2End,
            cigar: Cigar::new(0),
            alignment_end_pos: WavefrontPos::default(),
            status: AlignStatus::Ok,
            num_null_steps: 0,
            pattern_begin_free: 0,
            pattern_end_free: 0,
            text_begin_free: 0,
            text_end_free: 0,
            extension: false,
            compat: Compat::C,
            endsfree_best: None,
            heuristic: HeuristicStrategy::None,
            heuristic_state: HeuristicState::default(),
            memory_mode: MemoryMode::High,
            max_memory_compact: MAX_MEMORY_COMPACT,
            max_alignment_steps: i32::MAX,
            biwfa_fwd: None,
            biwfa_rev: None,
            biwfa_base: None,
        }
    }

    /// Set the memory mode used for alignments (see [`MemoryMode`]).
    pub fn set_memory_mode(&mut self, mode: MemoryMode) {
        self.memory_mode = mode;
    }

    /// Set the backtrace-buffer size (bytes) above which the med/low modes
    /// compact it (default 256 MB, fixed in C). Compaction only frees
    /// memory; it doesn't change results.
    pub fn set_max_memory_compact(&mut self, bytes: u64) {
        self.max_memory_compact = bytes;
    }

    /// C's backtrace-buffer check in `wavefront_unialign_reached_limits`,
    /// after computing `score`: every `PROBE_INTERVAL_COMPACT` scores,
    /// compact the buffer if it exceeds `max_memory_compact`. After
    /// `max_partial_compacts` compactions (4 for med, 1 for low) the next
    /// one starts over from the first block.
    #[inline(always)]
    fn probe_bt_buffer(&mut self, score: i32) {
        if !self.wf_components.bt_piggyback || score % PROBE_INTERVAL_COMPACT != 0 {
            return;
        }
        self.compact_bt_buffer_if_needed(score);
    }

    #[cold]
    fn compact_bt_buffer_if_needed(&mut self, score: i32) {
        let Some(buf) = self.wf_components.bt_buffer.as_ref() else {
            return;
        };
        if buf.get_size_used() <= self.max_memory_compact {
            return;
        }
        let mut pinned = self.endsfree_best.map_or(BT_BLOCK_IDX_NULL, |b| b.bt.1);
        self.wf_components
            .compact_bt_buffer(score as usize, &mut pinned);
        if let Some(b) = self.endsfree_best.as_mut() {
            b.bt.1 = pinned;
        }
        let max_partial_compacts = if self.memory_mode == MemoryMode::Low {
            1
        } else {
            4
        };
        let buf = self.wf_components.bt_buffer.as_mut().unwrap();
        if buf.num_compactions >= max_partial_compacts {
            buf.reset_compaction();
        }
    }

    /// Rebuild the components and slab with or without backtrace arrays and
    /// a backtrace buffer (piggyback backtrace for the med/low modes).
    fn set_piggyback(&mut self, piggyback: bool) {
        self.wf_components = WavefrontComponents::<N>::new(
            PATTERN_LENGTH_INIT,
            TEXT_LENGTH_INIT,
            &self.penalties,
            piggyback,
            piggyback,
        );
        let mut slab = WavefrontSlab::new(1000, piggyback, SlabMode::Reuse);
        if N == 5 {
            slab.set_cache_coloring(24);
        }
        self.wavefront_slab = slab;
    }

    /// Set the maximum alignment steps (score limit).
    pub fn set_max_alignment_steps(&mut self, max_steps: i32) {
        self.max_alignment_steps = max_steps;
    }

    /// Set the alignment scope.
    ///
    /// [`AlignmentScope::ComputeScore`] (the default) computes only the
    /// alignment score; [`AlignmentScope::ComputeAlignment`] additionally
    /// produces the full CIGAR (retrievable via [`Self::cigar`]).
    pub fn set_alignment_scope(&mut self, scope: AlignmentScope) {
        self.alignment_scope = scope;
    }

    /// Configure ends-free alignment with the given free-gap parameters.
    ///
    /// Up to `pattern_begin_free` / `text_begin_free` leading and
    /// `pattern_end_free` / `text_end_free` trailing bases of each sequence
    /// may be left unaligned at no cost. Use [`Self::align_endsfree`].
    ///
    /// # Known issues (parity with C WFA2-lib)
    ///
    /// Output is identical to C WFA2-lib, including these upstream issues:
    ///
    /// - **Edit/indel with begin gaps** (`*_begin_free > 0`) can return an
    ///   invalid CIGAR: a mismatching column reported as `M`. The score is
    ///   correct, and gap-linear/affine/affine2p are unaffected. This is a
    ///   regression in C WFA2-lib (from commit `380eb31b`); tracked in
    ///   <https://github.com/HPCBio/wfa2lib-rs/issues/18>. Fixed by
    ///   [`Compat::Corrected`].
    /// - **With a match reward** (`match < 0`), the chosen ends-free alignment
    ///   can be suboptimal (<https://github.com/Marco-Sola-Lab/WFA2-lib/issues/102>):
    ///   shifted by a repeat unit
    ///   (<https://github.com/HPCBio/wfa2lib-rs/issues/37>), or, when the free
    ///   ends span whole sequences, the all-gap alignment
    ///   (<https://github.com/HPCBio/wfa2lib-rs/issues/38>). Fixed by
    ///   [`Compat::Corrected`] (see [`Self::set_compat`]).
    pub fn set_alignment_free_ends(
        &mut self,
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_begin_free: i32,
        text_end_free: i32,
    ) {
        self.alignment_span = AlignmentSpan::EndsFree;
        self.extension = false;
        self.pattern_begin_free = pattern_begin_free;
        self.pattern_end_free = pattern_end_free;
        self.text_begin_free = text_begin_free;
        self.text_end_free = text_end_free;
    }

    /// Configure extension alignment mode.
    ///
    /// Extension alignment uses ends-free with full end-free parameters
    /// and applies maxtrim post-processing to find the maximal-scoring prefix.
    pub fn set_alignment_extension(&mut self) {
        self.alignment_span = AlignmentSpan::EndsFree;
        self.extension = true;
    }

    /// Choose C-compatible output (the default) or corrected behaviour for
    /// known C WFA2-lib bugs (see [`Compat`]).
    pub fn set_compat(&mut self, compat: Compat) {
        self.compat = compat;
    }

    /// Set the heuristic strategy for wavefront pruning.
    pub fn set_heuristic(&mut self, strategy: HeuristicStrategy) {
        self.heuristic = strategy;
    }

    /// Align two sequences end-to-end and return the alignment score.
    ///
    /// For edit distance, the score equals the edit distance (Levenshtein distance).
    /// For indel distance, the score equals the indel distance (LCS-based).
    /// For gap-linear/affine, the score equals the WFA score (sum of penalties).
    pub fn align_end2end(&mut self, pattern: &[u8], text: &[u8]) -> i32 {
        // Init sequences
        self.sequences.init_ascii(pattern, text, false);
        self.align_end2end_inner()
    }

    /// Shared end-to-end alignment loop (sequences must already be initialized).
    fn align_end2end_inner(&mut self) -> i32 {
        // Init alignment
        self.init_alignment();

        // Main alignment loop
        let mut score = 0;
        loop {
            // Extend
            if self.extend_end2end(score) {
                break;
            }

            // Next score
            score += 1;

            // Stop before writing past the score-indexed component arrays
            if !self.wf_components.has_slot(score as usize) {
                self.status = AlignStatus::MaxStepsReached;
                break;
            }

            // Compute next wavefront based on distance metric
            self.compute_step(score);

            // Check limits
            if score >= self.max_alignment_steps {
                self.status = AlignStatus::MaxStepsReached;
                break;
            }
            self.probe_bt_buffer(score);
        }

        // Record how many flat slots were used (for efficient clear on next alignment).
        self.wf_components.prev_score_count = (score + 1) as usize;

        // Finalize
        if self.status == AlignStatus::EndReached {
            // Backtrace to produce CIGAR if requested
            if self.alignment_scope == AlignmentScope::ComputeAlignment {
                self.backtrace_from(
                    score,
                    self.alignment_end_pos.k,
                    self.alignment_end_pos.offset,
                );
            }
            self.status = AlignStatus::Completed;
        } else if self.status == AlignStatus::EndUnreachable
            && self.alignment_end_pos.offset != OFFSET_NULL
        {
            // Z-dropped: partial alignment up to the best cell, max-trimmed
            // (C's `wavefront_unialign_terminate`, which backtraces with the
            // current score).
            if self.alignment_scope == AlignmentScope::ComputeAlignment {
                self.backtrace_from(
                    score,
                    self.alignment_end_pos.k,
                    self.alignment_end_pos.offset,
                );
                self.maxtrim_cigar();
            }
            self.status = AlignStatus::Partial;
        }

        score
    }

    /// Align two sequences with ends-free (semi-global) alignment.
    ///
    /// Ends-free alignment allows free gaps at the beginning and end of each
    /// sequence, as set by [`Self::set_alignment_free_ends`] (see its known
    /// issues) or [`Self::set_alignment_extension`].
    ///
    /// For extension mode, the CIGAR is trimmed to the maximal-scoring prefix.
    pub fn align_endsfree(&mut self, pattern: &[u8], text: &[u8]) -> i32 {
        // Init sequences
        self.sequences.init_ascii(pattern, text, false);

        // Extension: no begin gaps, free ends up to the full sequence lengths
        // (C's `wavefront_align_presets__checks`)
        if self.extension {
            self.pattern_begin_free = 0;
            self.pattern_end_free = self.sequences.pattern_length;
            self.text_begin_free = 0;
            self.text_end_free = self.sequences.text_length;
        }

        // Init alignment
        self.init_alignment();

        // Main alignment loop
        let mut score = 0;
        loop {
            // Extend with ends-free termination
            if self.extend_endsfree(score) {
                break;
            }

            // Next score
            score += 1;

            // Stop before writing past the score-indexed component arrays
            if !self.wf_components.has_slot(score as usize) {
                self.status = AlignStatus::MaxStepsReached;
                break;
            }

            // Compute next wavefront based on distance metric
            match self.penalties.distance_metric {
                DistanceMetric::Edit | DistanceMetric::Indel => {
                    self.compute_edit_step(score);
                }
                DistanceMetric::GapLinear => {
                    self.compute_linear_step(score);
                }
                DistanceMetric::GapAffine => {
                    self.compute_affine_step(score);
                }
                DistanceMetric::GapAffine2p => {
                    self.compute_affine2p_step(score);
                }
            }

            // Check limits
            if score >= self.max_alignment_steps {
                self.status = AlignStatus::MaxStepsReached;
                break;
            }
            self.probe_bt_buffer(score);
        }

        // Record how many flat slots were used (for efficient clear on next alignment).
        self.wf_components.prev_score_count = (score + 1) as usize;

        // Finalize
        let mut result = score;
        if self.status == AlignStatus::EndReached || self.status == AlignStatus::EndUnreachable {
            let unreachable = self.status == AlignStatus::EndUnreachable;

            // Corrected mode: the end candidate's cell can be from an earlier
            // score, and its total is the alignment's score.
            let candidate = self.endsfree_best.filter(|_| !unreachable);
            if let Some(best) = candidate {
                result = best.total;
            }

            // Backtrace to produce CIGAR if requested
            if self.alignment_scope == AlignmentScope::ComputeAlignment
                && (self.status == AlignStatus::EndReached
                    || self.alignment_end_pos.offset != OFFSET_NULL)
            {
                if let Some(best) = candidate {
                    let WavefrontPos { score, k, offset } = best.pos;
                    if self.wf_components.bt_piggyback {
                        self.backtrace_piggyback_bt(best.bt, score, k, offset);
                    } else {
                        self.backtrace_from(score, k, offset);
                    }
                    self.cigar.score = best.total;
                } else {
                    self.backtrace_from(
                        score,
                        self.alignment_end_pos.k,
                        self.alignment_end_pos.offset,
                    );
                }

                // Post-processing: extension mode or unreachable → maxtrim
                if self.extension || unreachable {
                    let trimmed = self.maxtrim_cigar();
                    self.status = if trimmed || unreachable {
                        AlignStatus::Partial
                    } else {
                        AlignStatus::Completed
                    };
                } else {
                    self.status = AlignStatus::Completed;
                }
            } else {
                self.status = if unreachable {
                    AlignStatus::Partial
                } else {
                    AlignStatus::Completed
                };
            }
        }

        result
    }

    /// Get the alignment status after the last alignment.
    pub fn status(&self) -> AlignStatus {
        self.status
    }

    /// Get the CIGAR from the last alignment (only valid after ComputeAlignment).
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    // --- Internal methods ---

    /// Initialize for a new alignment.
    fn init_alignment(&mut self) {
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;

        // Resize CIGAR buffer if computing alignment
        if self.alignment_scope == AlignmentScope::ComputeAlignment {
            self.cigar.resize(cigar_buffer_len(plen, tlen));
        }

        // Enable modular wavefront storage for score-only mode (like C reference).
        // This keeps only max_score_scope wavefronts alive, dramatically reducing
        // memory usage and improving cache performance. Med/low alignments
        // are modular too, with piggyback backtrace (C's
        // `wavefront_aligner_init_alignment`).
        let score_only = self.alignment_scope == AlignmentScope::ComputeScore;
        let piggyback = !score_only && self.memory_mode != MemoryMode::High;
        if piggyback != self.wf_components.bt_piggyback {
            self.set_piggyback(piggyback);
        }
        let score_only = score_only || piggyback;
        self.wf_components.memory_modular = score_only;

        // Resize components for sequence dimensions
        self.wf_components.resize(plen, tlen, &self.penalties);

        // As C: modular (score-only) wavefronts all have the slab's current
        // length and are reused; non-modular ones are sized to their own
        // range ("tight"), so memory follows the wavefronts' actual spread.
        // Score-only (modular): full clear; circular window already handles reuse.
        self.wavefront_slab.set_mode(if score_only { SlabMode::Reuse } else { SlabMode::Tight });
        if score_only {
            self.wavefront_slab.clear();
        } else {
            self.wavefront_slab.reap_cigar();
        }

        // Pre-size slab to avoid Vec reallocations that would invalidate stored *mut Wavefront ptrs.
        // In CIGAR mode: reserve for all wavefronts across all scores.
        // In score-only (modular) mode: reserve for the circular window (max_score_scope * num_types).
        {
            let num_types = match self.penalties.distance_metric {
                DistanceMetric::Edit | DistanceMetric::Indel => 1,
                DistanceMetric::GapLinear => 1,
                DistanceMetric::GapAffine => 3,
                DistanceMetric::GapAffine2p => 5,
            };
            let capacity = if score_only {
                // Modular: at most max_score_scope * num_types live wavefronts at once
                self.wf_components.max_score_scope * num_types + num_types
            } else {
                self.wf_components.num_wavefronts * num_types
            };
            self.wavefront_slab.reserve(capacity);
        }

        // Reset state
        self.num_null_steps = 0;
        self.status = AlignStatus::Ok;
        // No end position until the end is reached or a Z-drop sets one
        // (initialized as in C's `wavefront_aligner_init`).
        self.alignment_end_pos = WavefrontPos {
            score: -1,
            k: DPMATRIX_DIAGONAL_NULL,
            offset: OFFSET_NULL,
        };
        self.heuristic_state.reset(&self.heuristic);
        self.endsfree_best = None;

        // Compute allocation limits based on sequence lengths.
        // k ranges from -plen to +tlen, so we need wavefronts that cover this
        // full range. Add padding for the kernel's k±1 reads and score scope.
        let max_score_scope = self.wf_components.max_score_scope as i32;
        let eff_lo = -(plen + max_score_scope + 2);
        let eff_hi = tlen + max_score_scope + 2;

        // Limits the computed lo/hi are clamped to.
        self.wf_components.historic_min_lo = eff_lo;
        self.wf_components.historic_max_hi = eff_hi;

        // Allocate initial M-wavefront at score 0, over C's initial range
        // (the begin-free diagonals with match == 0).
        let (lo0, hi0) =
            if self.penalties.match_ == 0 && self.alignment_span == AlignmentSpan::EndsFree {
                (-self.pattern_begin_free, self.text_begin_free)
            } else {
                (0, 0)
            };
        let (alloc_lo, alloc_hi) = self.output_alloc_range(0, lo0, hi0);
        let wf_ptr = self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi);
        unsafe {
            (*wf_ptr).set_offset(0, 0);
            (*wf_ptr).set_limits(0, 0);
        }
        if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
            buf.clear();
            unsafe { set_bt_init(&mut *wf_ptr, buf, 0, 0, 0) };
        }
        // Ends-free with match == 0: the begin-free cells (h,0) and (0,v) are
        // all reachable at score 0 (C's `wavefront_aligner_init_wf_m`, with
        // C #126's restriction to ends-free span). For match != 0 they are
        // added during compute instead (`endsfree_begin_k`).
        if self.alignment_span == AlignmentSpan::EndsFree && self.penalties.match_ == 0 {
            let (t0, p0) = (self.text_begin_free, self.pattern_begin_free);
            unsafe {
                for h in 1..=t0 {
                    (*wf_ptr).set_offset(h, h); // k = h, offset = h
                    if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
                        set_bt_init(&mut *wf_ptr, buf, h, 0, h);
                    }
                }
                for v in 1..=p0 {
                    (*wf_ptr).set_offset(-v, 0); // k = -v, offset = h = 0
                    if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
                        set_bt_init(&mut *wf_ptr, buf, -v, v, 0);
                    }
                }
                // C sets lo/hi only; the initialized range stays [0, 0]
                // (widened later by the gap-metric `init_ends`, never by edit).
                (*wf_ptr).lo = -p0;
                (*wf_ptr).hi = t0;
                // Corrected mode: the seeds are initialized too, so the
                // edit/indel backtrace can step onto them (HPCBio/wfa2lib-rs#18,
                // the fix proposed in Marco-Sola-Lab/WFA2-lib#142).
                if self.compat == Compat::Corrected {
                    (*wf_ptr).wf_elements_init_min = -p0;
                    (*wf_ptr).wf_elements_init_max = t0;
                }
                // The seeds' backtraces are valid (read by med/low when the
                // alignment ends on a seed's diagonal at score 0).
                if self.wf_components.bt_buffer.is_some() {
                    (*wf_ptr).mark_written(-p0, t0);
                }
            }
        }
        self.wf_components.set_m_ptr(0, wf_ptr);
    }

    /// Extend the M-wavefront at the given score, check termination.
    /// Returns true if alignment is finished.
    fn extend_end2end(&mut self, score: i32) -> bool {
        let score_idx = score as usize;
        let m_ptr = self.wf_components.get_m_ptr(score_idx);

        if m_ptr.is_null() {
            // Null steps are counted by the compute steps (as in C).
            if self.num_null_steps > self.wf_components.max_score_scope as i32 {
                self.status = AlignStatus::EndUnreachable;
                return true;
            }
            return false;
        }

        // Extend matches
        let wf = unsafe { &mut *m_ptr };
        if self.sequences.mode == SequenceMode::Lambda {
            extend::extend_matches_custom_end2end(&self.sequences, wf);
        } else {
            extend::extend_matches_packed_end2end(&self.sequences, wf);
        }

        // Check termination
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;
        let wf = unsafe { &*m_ptr };
        if termination::termination_end2end(plen, tlen, wf) {
            let alignment_k = tlen - plen;
            self.alignment_end_pos = WavefrontPos {
                score,
                k: alignment_k,
                offset: tlen,
            };
            self.status = AlignStatus::EndReached;
            return true;
        }

        // Apply heuristic cutoff
        if !matches!(self.heuristic, HeuristicStrategy::None)
            && self.apply_heuristic_cutoff(score, m_ptr)
        {
            self.status = AlignStatus::EndUnreachable;
            return true;
        }

        false
    }

    /// Extend the M-wavefront at the given score with ends-free termination.
    /// Returns true if alignment is finished.
    fn extend_endsfree(&mut self, score: i32) -> bool {
        if self.endsfree_corrected() {
            return self.extend_endsfree_corrected(score);
        }
        let score_idx = score as usize;
        let m_ptr = self.wf_components.get_m_ptr(score_idx);

        if m_ptr.is_null() {
            // Null steps are counted by the compute steps (as in C).
            if self.num_null_steps > self.wf_components.max_score_scope as i32 {
                self.status = AlignStatus::EndUnreachable;
                return true;
            }
            return false;
        }

        // Extend matches with ends-free termination check
        let wf = unsafe { &mut *m_ptr };
        let result = if self.sequences.mode == SequenceMode::Lambda {
            extend::extend_matches_custom_endsfree(
                &self.sequences,
                wf,
                self.pattern_end_free,
                self.text_end_free,
            )
        } else {
            extend::extend_matches_packed_endsfree(
                &self.sequences,
                wf,
                self.pattern_end_free,
                self.text_end_free,
            )
        };

        if let Some((k, offset)) = result {
            self.alignment_end_pos = WavefrontPos { score, k, offset };
            self.status = AlignStatus::EndReached;
            return true;
        }

        // Apply heuristic cutoff
        if !matches!(self.heuristic, HeuristicStrategy::None)
            && self.apply_heuristic_cutoff(score, m_ptr)
        {
            self.status = AlignStatus::EndUnreachable;
            return true;
        }

        false
    }

    /// Whether ends-free termination uses the corrected rule (see
    /// [`Compat::Corrected`]); with `match == 0` both rules agree.
    #[inline(always)]
    fn endsfree_corrected(&self) -> bool {
        self.compat == Compat::Corrected && self.penalties.match_ < 0
    }

    /// Corrected-mode [`Self::extend_endsfree`] (see [`Compat::Corrected`]).
    ///
    /// With a match reward the WFA score of an ends-free alignment is
    /// `2 * penalty + (-match) * (aligned bases)`, so on its own it doesn't
    /// rank alignments that leave different numbers of bases to free end
    /// gaps. Counting each free-gap base as `-match`, as C already does for
    /// begin-free ones (`endsfree_begin_k`), makes the total
    /// `2 * penalty + (-match) * (|P| + |T|)` for every alignment. A cell
    /// reaching a free end at `score` with `left` bases to skip is then a
    /// complete alignment of total `score + (-match) * left`. Every later
    /// cell totals more than the current score, so the best candidate is
    /// optimal once the current score reaches its total.
    fn extend_endsfree_corrected(&mut self, score: i32) -> bool {
        let m_ptr = self.wf_components.get_m_ptr(score as usize);
        let prune = matches!(
            self.heuristic,
            HeuristicStrategy::None | HeuristicStrategy::BandedStatic { .. }
        );
        if prune
            && let Some(best) = self.endsfree_best
            && best.total > score
        {
            self.endsfree_prune(score, best.total);
        }

        if !m_ptr.is_null() {
            let wf = unsafe { &mut *m_ptr };
            let found = if self.sequences.mode == SequenceMode::Lambda {
                extend::extend_matches_custom_endsfree_best(
                    &self.sequences,
                    wf,
                    self.pattern_end_free,
                    self.text_end_free,
                )
            } else {
                extend::extend_matches_packed_endsfree_best(
                    &self.sequences,
                    wf,
                    self.pattern_end_free,
                    self.text_end_free,
                )
            };
            if let Some((k, offset, left)) = found {
                let total = score + (-self.penalties.match_) * left;
                if self.endsfree_best.is_none_or(|b| total < b.total) {
                    let bt = if self.wf_components.bt_piggyback {
                        wf.bt_or_zero(k)
                    } else {
                        (0, BT_BLOCK_IDX_NULL)
                    };
                    self.endsfree_best = Some(EndsfreeCandidate {
                        total,
                        pos: WavefrontPos { score, k, offset },
                        bt,
                    });
                    if prune && total > score {
                        self.endsfree_prune(score, total);
                    }
                }
            }
        }

        // No later cell can beat the candidate
        if let Some(best) = self.endsfree_best
            && best.total <= score
        {
            self.alignment_end_pos = best.pos;
            self.status = AlignStatus::EndReached;
            return true;
        }

        let done = if m_ptr.is_null() {
            // Null steps are counted by the compute steps (as in C). With a
            // candidate the search is already bounded by its total, and null
            // steps don't end it: pruning can empty every wavefront while
            // begin-free seeds are still to come.
            self.endsfree_best.is_none()
                && self.num_null_steps > self.wf_components.max_score_scope as i32
        } else {
            !matches!(self.heuristic, HeuristicStrategy::None)
                && self.apply_heuristic_cutoff(score, m_ptr)
        };
        if done {
            // Nothing left to search (or pruned by the heuristic): the
            // candidate, if any, is the alignment.
            if let Some(best) = self.endsfree_best {
                self.alignment_end_pos = best.pos;
                self.status = AlignStatus::EndReached;
            } else {
                self.status = AlignStatus::EndUnreachable;
            }
            return true;
        }

        false
    }

    /// Corrected-mode pruning: drop the cells at `score` that can't complete
    /// to a total below the candidate's `total`.
    ///
    /// From cell (v, h) on diagonal k, the bases left are `|P| - v` and
    /// `|T| - h`, which differ by `k - k_end` (`k_end = |T| - |P|`) all along
    /// the diagonal. Matches cost nothing, but each base of that difference
    /// becomes a gap base (at least `-match`) or a free-end base (`-match`),
    /// so any completion totals at least `score + (-match) * |k - k_end|`.
    /// Diagonals where that exceeds `total` can't give a better end (ties are
    /// kept), and neither can any cell computed only from them, so the
    /// result doesn't change. The kept diagonals form a band around `k_end`
    /// that narrows as the score grows.
    ///
    /// Only with no heuristic or `BandedStatic`: the other heuristics read
    /// the wavefronts, so pruning could change what they cut.
    fn endsfree_prune(&mut self, score: i32, total: i32) {
        let radius = (total - score) / -self.penalties.match_;
        let k_end = self.sequences.text_length - self.sequences.pattern_length;
        let (lo, hi) = (k_end - radius, k_end + radius);
        let score_idx = score as usize;
        let getters: &[fn(&WavefrontComponents<N>, usize) -> *mut crate::wavefront::Wavefront] =
            match self.penalties.distance_metric {
                DistanceMetric::GapAffine => &[
                    WavefrontComponents::get_m_ptr,
                    WavefrontComponents::get_i1_ptr,
                    WavefrontComponents::get_d1_ptr,
                ],
                DistanceMetric::GapAffine2p => &[
                    WavefrontComponents::get_m_ptr,
                    WavefrontComponents::get_i1_ptr,
                    WavefrontComponents::get_d1_ptr,
                    WavefrontComponents::get_i2_ptr,
                    WavefrontComponents::get_d2_ptr,
                ],
                _ => &[WavefrontComponents::get_m_ptr],
            };
        for get_ptr in getters {
            let ptr = get_ptr(&self.wf_components, score_idx);
            if !ptr.is_null() {
                unsafe { (*ptr).trim_nulled(lo, hi) };
            }
        }
    }

    /// Backtrace the CIGAR from M-wavefront cell (`k`, `offset`), starting at
    /// `score` (C's `wavefront_unialign_terminate`).
    fn backtrace_from(&mut self, score: i32, k: i32, offset: WfOffset) {
        if self.wf_components.bt_piggyback {
            self.backtrace_piggyback(score, k, offset);
        } else {
            self.backtrace_components(score, k, offset, ComponentType::M, ComponentType::M);
        }
    }

    /// C's `wavefront_backtrace_pcigar` (med/low): follow the chain of
    /// packed CIGARs from cell `k` of the M-wavefront at `score` back to its
    /// initial position, then unpack it, re-extending matches. As in C, the
    /// cell is read at the current score even for a Z-drop end found earlier.
    fn backtrace_piggyback(&mut self, score: i32, k: i32, offset: WfOffset) {
        let m = self.wf_components.get_m_ptr(score as usize);
        // C reads the cell unchecked. After a Z-drop it may lie outside
        // [lo, hi]: cells the kernel computed before trimming are as in C,
        // but never-written ones hold leftovers of earlier alignments in C
        // (zeroes in fresh memory), so C's result depends on what it aligned
        // before (HPCBio/wfa2lib-rs#22). Rust reads them as zeroes, so each
        // result is deterministic. A NULL wavefront crashes C; Rust returns
        // an empty CIGAR.
        if m.is_null() || !unsafe { (*m).is_allocated(k) } {
            self.cigar.begin_offset = 0;
            self.cigar.end_offset = 0;
            return;
        }
        let bt = unsafe { (*m).bt_or_zero(k) };
        self.backtrace_piggyback_bt(bt, score, k, offset);
    }

    /// Unpack the packed-CIGAR chain `(pcigar, prev)` of M-wavefront cell
    /// (`k`, `offset`) at `score` into the CIGAR.
    fn backtrace_piggyback_bt(
        &mut self,
        (pcigar, prev): (Pcigar, BtBlockIdx),
        score: i32,
        k: i32,
        offset: WfOffset,
    ) {
        let buf = self.wf_components.bt_buffer.as_mut().unwrap();
        let init = buf.traceback_from(pcigar, prev);
        let init_pos = buf.alignment_init_pos[init.pcigar as usize];
        let affine = matches!(
            self.penalties.distance_metric,
            DistanceMetric::GapAffine | DistanceMetric::GapAffine2p
        );
        unpack_piggyback(
            &buf.alignment_packed,
            &self.sequences,
            (init_pos.v, init_pos.h),
            (wavefront_v(k, offset), wavefront_h(k, offset)),
            affine,
            &mut self.cigar,
        );
        self.cigar.score = score;
    }

    /// `backtrace_from` between the given begin/end components (BiWFA halves).
    fn backtrace_components(
        &mut self,
        score: i32,
        k: i32,
        offset: WfOffset,
        component_begin: ComponentType,
        component_end: ComponentType,
    ) {
        let (plen, tlen) = (self.sequences.pattern_length, self.sequences.text_length);
        match self.penalties.distance_metric {
            DistanceMetric::GapAffine => backtrace::backtrace_affine(
                &self.wf_components,
                &self.wavefront_slab,
                &self.penalties,
                plen,
                tlen,
                score,
                k,
                offset,
                component_begin,
                component_end,
                &mut self.cigar,
            ),
            DistanceMetric::GapAffine2p => backtrace::backtrace_affine2p(
                &self.wf_components,
                &self.wavefront_slab,
                &self.penalties,
                plen,
                tlen,
                score,
                k,
                offset,
                component_begin,
                component_end,
                &mut self.cigar,
            ),
            _ => backtrace::backtrace_linear(
                &self.wf_components,
                &self.wavefront_slab,
                &self.penalties,
                plen,
                tlen,
                score,
                k,
                offset,
                &mut self.cigar,
            ),
        }
    }

    /// Apply maxtrim to the CIGAR based on the distance metric.
    /// Returns true if the CIGAR was trimmed.
    fn maxtrim_cigar(&mut self) -> bool {
        match self.penalties.distance_metric {
            DistanceMetric::GapLinear => {
                let p = self.penalties.linear_penalties.as_ref().unwrap();
                self.cigar.maxtrim_gap_linear(p)
            }
            DistanceMetric::GapAffine => {
                let p = self.penalties.affine_penalties.as_ref().unwrap();
                self.cigar.maxtrim_gap_affine(p)
            }
            DistanceMetric::GapAffine2p => {
                let p = self.penalties.affine2p_penalties.as_ref().unwrap();
                self.cigar.maxtrim_gap_affine2p(p)
            }
            _ => false, // Edit/indel don't support maxtrim
        }
    }

    /// Apply heuristic cutoff to the M-wavefront (and I/D wavefronts if affine).
    /// Returns true if alignment should be terminated (z-drop).
    fn apply_heuristic_cutoff(&mut self, score: i32, m_ptr: *mut crate::wavefront::Wavefront) -> bool {
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;

        // Use the victim wavefront's offsets as scratch buffer for distances
        let distances_base_k = self.wf_components.wavefront_victim.base_k();
        let distances = self.wf_components.wavefront_victim.offsets_slice_mut();

        let wf = unsafe { &mut *m_ptr };
        let (lo_base, hi_base) = (wf.lo, wf.hi);
        let zdropped = heuristic::heuristic_cutoff(
            &self.heuristic,
            &mut self.heuristic_state,
            wf,
            distances,
            distances_base_k,
            plen,
            tlen,
            self.penalties.match_,
            score,
        );

        if zdropped {
            // C's `wavefront_heuristic_zdrop`: the alignment ends at the best
            // cell seen; it is backtraced as a partial alignment.
            let st = &self.heuristic_state;
            self.alignment_end_pos = WavefrontPos {
                score: st.max_wf_score,
                k: st.max_sw_score_k,
                offset: st.max_sw_score_offset,
            };
            return true;
        }

        // As C's `wavefront_heuristic_cufoff`: an empty M is left alone, and I/D
        // are only equated when the cutoff actually narrowed M. (Equating
        // unconditionally empties I/D whenever M is empty, ending alignments
        // that C completes through the I/D wavefronts.)
        if lo_base > hi_base || (wf.lo == lo_base && wf.hi == hi_base) {
            return false;
        }
        if wf.lo > wf.hi {
            wf.null = true;
        }

        // Propagate narrowed M-wavefront bounds to I/D wavefronts
        let score_idx = score as usize;
        match self.penalties.distance_metric {
            DistanceMetric::GapAffine => {
                let (m_lo, m_hi) = unsafe { ((*m_ptr).lo, (*m_ptr).hi) };

                let i1_ptr = self.wf_components.get_i1_ptr(score_idx);
                if !i1_ptr.is_null() {
                    unsafe {
                        if m_lo > (*i1_ptr).lo { (*i1_ptr).lo = m_lo; }
                        if m_hi < (*i1_ptr).hi { (*i1_ptr).hi = m_hi; }
                        if (*i1_ptr).lo > (*i1_ptr).hi { (*i1_ptr).null = true; }
                    }
                }

                let d1_ptr = self.wf_components.get_d1_ptr(score_idx);
                if !d1_ptr.is_null() {
                    unsafe {
                        if m_lo > (*d1_ptr).lo { (*d1_ptr).lo = m_lo; }
                        if m_hi < (*d1_ptr).hi { (*d1_ptr).hi = m_hi; }
                        if (*d1_ptr).lo > (*d1_ptr).hi { (*d1_ptr).null = true; }
                    }
                }
            }
            DistanceMetric::GapAffine2p => {
                let (m_lo, m_hi) = unsafe { ((*m_ptr).lo, (*m_ptr).hi) };

                for get_ptr_fn in [
                    WavefrontComponents::get_i1_ptr,
                    WavefrontComponents::get_d1_ptr,
                    WavefrontComponents::get_i2_ptr,
                    WavefrontComponents::get_d2_ptr,
                ] {
                    let ptr = get_ptr_fn(&self.wf_components, score_idx);
                    if !ptr.is_null() {
                        unsafe {
                            if m_lo > (*ptr).lo { (*ptr).lo = m_lo; }
                            if m_hi < (*ptr).hi { (*ptr).hi = m_hi; }
                            if (*ptr).lo > (*ptr).hi { (*ptr).null = true; }
                        }
                    }
                }
            }
            _ => {} // Edit/indel/linear only have M-wavefront
        }

        false
    }

    /// Compute the next wavefront for edit/indel distance.
    fn compute_edit_step(&mut self, score: i32) {
        let score_prev = (score - 1) as usize;
        let score_curr = score as usize;

        let prev_ptr = self.wf_components.get_m_ptr(score_prev);
        // A missing or trimmed-empty previous wavefront yields a null output, so
        // the extend step counts null steps toward END_UNREACHABLE. (C computes
        // from it anyway and never terminates; here the range would drift by one
        // diagonal per step past the allocation.)
        if is_null_wf(prev_ptr) {
            self.num_null_steps += 1;
            self.set_null_output_m(score_curr);
            return;
        }
        self.num_null_steps = 0;
        // Get prev wavefront dimensions
        let (prev_lo, prev_hi) = unsafe { ((*prev_ptr).lo, (*prev_ptr).hi) };

        // Clamp lo/hi to leave room for kernel's k±1 reads within allocated range.
        let hist_lo = self.wf_components.historic_min_lo;
        let hist_hi = self.wf_components.historic_max_hi;
        let lo = (prev_lo - 1).max(hist_lo + 1);
        let hi = (prev_hi + 1).min(hist_hi - 1);
        if lo > hi {
            self.set_null_output_m(score_curr);
            return;
        }

        // Null the two cells beyond each end of prev that the kernel reads
        // (prev[k-1], prev[k], prev[k+1] for k in [prev.lo - 1, prev.hi + 1]).
        // As in C's `wavefront_compute_edit`, these are written directly and
        // prev's initialized range is left unchanged: the backtrace bounds
        // checks use it, so widening it would let the backtrace step onto
        // M0's begin-free diagonals, which C's does not.
        unsafe {
            for k in [prev_lo - 2, prev_lo - 1, prev_hi + 1, prev_hi + 2] {
                (*prev_ptr).set_null_cell_if_allocated(k);
            }
        }

        // Allocate or reuse output wavefront
        let (alloc_lo, alloc_hi) = self.output_alloc_range(score, lo, hi);
        let old = if self.wf_components.memory_modular {
            self.wf_components.get_m_ptr(score_curr)
        } else {
            WF_PTR_NONE
        };
        let curr_ptr = self.reuse_or_allocate(old, alloc_lo, alloc_hi);
        unsafe { (*curr_ptr).set_limits(lo, hi); }

        // Compute kernel
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;
        let indel = match self.penalties.distance_metric {
            DistanceMetric::Edit => false,
            DistanceMetric::Indel => true,
            _ => unreachable!("compute_edit_step called for non-edit metric"),
        };
        unsafe {
            match (self.wf_components.bt_buffer.as_mut(), indel) {
                (None, false) => {
                    edit::compute_edit_idm(plen, tlen, &*prev_ptr, &mut *curr_ptr, lo, hi)
                }
                (None, true) => {
                    edit::compute_indel_idm(plen, tlen, &*prev_ptr, &mut *curr_ptr, lo, hi)
                }
                (Some(buf), _) => {
                    if indel {
                        piggyback::compute_indel_idm(
                            plen,
                            tlen,
                            &*prev_ptr,
                            &mut *curr_ptr,
                            lo,
                            hi,
                        );
                    } else {
                        piggyback::compute_edit_idm(plen, tlen, &*prev_ptr, &mut *curr_ptr, lo, hi);
                    }
                    // C offloads edit/indel every PCIGAR_MAX_LENGTH scores.
                    if score % PCIGAR_MAX_LENGTH as i32 == 0 {
                        let low = self.memory_mode == MemoryMode::Low;
                        piggyback::offload_blocks_linear(&mut *curr_ptr, lo, hi, low, buf);
                    }
                }
            }
        }

        // Store in components
        self.wf_components.set_m_ptr(score_curr, curr_ptr);

        // Trim ends
        {
            let wf_curr = unsafe { &mut *curr_ptr };
            compute::trim_ends(plen, tlen, wf_curr);
        }
    }
    /// Compute the next wavefront for gap-linear distance.
    fn compute_linear_step(&mut self, score: i32) {
        let score_curr = score as usize;
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;

        // Fetch input wavefronts
        let misms_score = score - self.penalties.mismatch;
        let gap_score = score - self.penalties.gap_opening1;

        let misms_ptr = if misms_score >= 0 {
            self.wf_components.get_m_ptr(misms_score as usize)
        } else {
            WF_PTR_NONE
        };
        let gap_ptr = if gap_score >= 0 {
            self.wf_components.get_m_ptr(gap_score as usize)
        } else {
            WF_PTR_NONE
        };

        // All inputs missing or trimmed-empty: null output (C's
        // `wavefront_compute_allocate_output_null`), counted as a null step by
        // extend. Computing from them would widen the range every step.
        if is_null_wf(misms_ptr) && is_null_wf(gap_ptr) {
            self.num_null_steps += 1;
            self.set_null_output_m(score_curr);
            if let Some(k) = self.endsfree_begin_k(score) {
                let wf = self.endsfree_null_m(k);
                self.wf_components.set_m_ptr(score_curr, wf);
            }
            return;
        }
        self.num_null_steps = 0;

        // Trimmed-empty inputs are read as the null wavefront, like missing ones,
        // whose bounds (lo = 1, hi = -1) still enter the range union as in C.
        let misms_ptr = live_or_none(misms_ptr);
        let gap_ptr = live_or_none(gap_ptr);
        let (misms_lo, misms_hi) = null_wf_bounds(misms_ptr);
        let (gap_lo, gap_hi) = null_wf_bounds(gap_ptr);

        // Union of ranges (gap expands by ±1 for ins/del k-shifts), clamped to
        // leave room for the kernel's k±1 reads within the allocated range.
        let hist_lo = self.wf_components.historic_min_lo;
        let hist_hi = self.wf_components.historic_max_hi;
        let (lo, hi) = self.clamp_to_band(misms_lo.min(gap_lo - 1), misms_hi.max(gap_hi + 1));
        let lo = lo.max(hist_lo + 1);
        let hi = hi.min(hist_hi - 1);

        if lo > hi {
            self.set_null_output_m(score_curr);
            return;
        }

        // Allocation range (also resizes the null/victim wavefronts)
        let (alloc_lo, alloc_hi) = self.output_alloc_range(score, lo, hi);

        // Initialize input wavefront ends (like C's wavefront_compute_init_ends)
        // misms needs [lo, hi]
        if !misms_ptr.is_null() {
            unsafe {
                (*misms_ptr).init_ends_higher(hi);
                (*misms_ptr).init_ends_lower(lo);
            }
        }
        // gap needs [lo-1, hi+1] (kernel reads k-1 and k+1)
        if !gap_ptr.is_null() {
            unsafe {
                (*gap_ptr).init_ends_higher(hi + 1);
                (*gap_ptr).init_ends_lower(lo - 1);
            }
        }

        // Allocate or reuse output wavefront
        let old = if self.wf_components.memory_modular {
            self.wf_components.get_m_ptr(score_curr)
        } else {
            WF_PTR_NONE
        };
        let curr_ptr = self.reuse_or_allocate(old, alloc_lo, alloc_hi);
        unsafe { (*curr_ptr).set_limits(lo, hi); }

        // Compute kernel using raw pointers for multi-borrow
        unsafe {
            let null_wf = &self.wf_components.wavefront_null as *const _;

            let wf_misms_cptr = if !misms_ptr.is_null() { misms_ptr as *const _ } else { null_wf };
            let wf_gap_cptr = if !gap_ptr.is_null() { gap_ptr as *const _ } else { null_wf };

            if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
                piggyback::compute_linear_idm(
                    plen,
                    tlen,
                    &*wf_misms_cptr,
                    &*wf_gap_cptr,
                    &mut *curr_ptr,
                    lo,
                    hi,
                );
                let low = self.memory_mode == MemoryMode::Low;
                piggyback::offload_linear(
                    &*wf_misms_cptr,
                    &*wf_gap_cptr,
                    &mut *curr_ptr,
                    lo,
                    hi,
                    low,
                    buf,
                );
            } else {
                linear::compute_linear_idm(
                    plen,
                    tlen,
                    &*wf_misms_cptr,
                    &*wf_gap_cptr,
                    &mut *curr_ptr,
                    lo,
                    hi,
                );
            }
        }

        // Store in components
        self.wf_components.set_m_ptr(score_curr, curr_ptr);

        // Trim ends
        {
            let wf_curr = unsafe { &mut *curr_ptr };
            if let Some(k) = self.endsfree_begin_k(score) {
                self.endsfree_init(wf_curr, k);
            }
            compute::trim_ends(plen, tlen, wf_curr);
        }
    }

    /// Compute the next wavefront for gap-affine distance.
    fn compute_affine_step(&mut self, score: i32) {
        let score_curr = score as usize;
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;

        let x = self.penalties.mismatch;
        let o_plus_e = self.penalties.gap_opening1 + self.penalties.gap_extension1;
        let e = self.penalties.gap_extension1;

        // Pre-compute output slot once (avoids redundant modular divisions).
        let out_slot = self.wf_components.slot_of(score_curr);

        // Fetch input wavefront pointers directly from flat array.
        // Compute input slots once per distinct score to avoid redundant divs.
        let ptr_m_misms = if score - x >= 0 { self.wf_components.get_m_ptr((score - x) as usize) } else { WF_PTR_NONE };
        let ptr_m_open  = if score - o_plus_e >= 0 { self.wf_components.get_m_ptr((score - o_plus_e) as usize) } else { WF_PTR_NONE };
        let (ptr_i1_ext, ptr_d1_ext) = if score - e >= 0 {
            let s = self.wf_components.slot_of((score - e) as usize);
            (self.wf_components.get_i1_ptr_slot(s), self.wf_components.get_d1_ptr_slot(s))
        } else {
            (WF_PTR_NONE, WF_PTR_NONE)
        };

        // All inputs missing or trimmed-empty: null output (C's
        // `wavefront_compute_allocate_output_null`), counted as a null step by extend.
        if is_null_wf(ptr_m_misms)
            && is_null_wf(ptr_m_open)
            && is_null_wf(ptr_i1_ext)
            && is_null_wf(ptr_d1_ext)
        {
            self.num_null_steps += 1;
            self.free_output_wavefronts_affine(out_slot);
            if let Some(k) = self.endsfree_begin_k(score) {
                let wf = self.endsfree_null_m(k);
                self.wf_components.set_m_ptr_slot(out_slot, wf);
            }
            return;
        }
        self.num_null_steps = 0;

        // Trimmed-empty inputs are read as the null wavefront, like missing ones.
        let ptr_m_misms = live_or_none(ptr_m_misms);
        let ptr_m_open = live_or_none(ptr_m_open);
        let ptr_i1_ext = live_or_none(ptr_i1_ext);
        let ptr_d1_ext = live_or_none(ptr_d1_ext);

        // Effective lo/hi over all inputs. As in C's `wavefront_compute_limits_input`,
        // null inputs contribute the null wavefront's bounds (lo = 1, hi = -1).
        let (misms_lo, misms_hi) = null_wf_bounds(ptr_m_misms);
        let (open_lo, open_hi) = null_wf_bounds(ptr_m_open);
        let (i1_lo, i1_hi) = null_wf_bounds(ptr_i1_ext);
        let (d1_lo, d1_hi) = null_wf_bounds(ptr_d1_ext);
        let (mut lo, mut hi) = self.clamp_to_band(
            misms_lo.min(open_lo - 1).min(i1_lo + 1).min(d1_lo - 1),
            misms_hi.max(open_hi + 1).max(i1_hi + 1).max(d1_hi - 1),
        );

        if lo > hi {
            if self.wf_components.memory_modular {
                self.free_output_wavefronts_affine(out_slot);
            }
            return;
        }

        // Clamp lo/hi to leave room for kernel's k±1 reads within allocated range.
        let hist_lo = self.wf_components.historic_min_lo;
        let hist_hi = self.wf_components.historic_max_hi;
        lo = lo.max(hist_lo + 1);
        hi = hi.min(hist_hi - 1);

        if lo > hi {
            if self.wf_components.memory_modular {
                self.free_output_wavefronts_affine(out_slot);
            }
            return;
        }

        // Allocation range (also resizes the null/victim wavefronts)
        let (alloc_lo, alloc_hi) = self.output_alloc_range(score, lo, hi);

        // Prefetch input wavefront offset arrays into L1 cache.
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            macro_rules! prefetch_wf {
                ($ptr:expr, $lo:expr) => {
                    if !$ptr.is_null() {
                        let base = (*$ptr).centered_offsets.offset($lo as isize) as *const i8;
                        _mm_prefetch(base, _MM_HINT_T0);
                        _mm_prefetch(base.add(64), _MM_HINT_T0);
                    }
                }
            }
            prefetch_wf!(ptr_m_misms, lo);
            prefetch_wf!(ptr_m_open, lo);
            prefetch_wf!(ptr_i1_ext, lo);
            prefetch_wf!(ptr_d1_ext, lo);
        }

        // Initialize input wavefront ends using pointers (no slab re-lookup)
        unsafe {
            if !ptr_m_misms.is_null() { (*ptr_m_misms).init_ends_higher(hi);     (*ptr_m_misms).init_ends_lower(lo);     }
            if !ptr_m_open.is_null()  { (*ptr_m_open).init_ends_higher(hi + 1);  (*ptr_m_open).init_ends_lower(lo - 1);  }
            if !ptr_i1_ext.is_null()  { (*ptr_i1_ext).init_ends_higher(hi);      (*ptr_i1_ext).init_ends_lower(lo - 1);  }
            if !ptr_d1_ext.is_null()  { (*ptr_d1_ext).init_ends_higher(hi + 1);  (*ptr_d1_ext).init_ends_lower(lo);      }
        }

        // Allocate or reuse 3 output wavefronts using pre-computed slot.
        let (m_curr_ptr, i1_curr_ptr, d1_curr_ptr) = if self.wf_components.memory_modular {
            macro_rules! reuse_or_alloc_slot {
                ($comp:ident, $slot:expr) => {{
                    let old = self.wf_components.$comp($slot);
                    self.reuse_or_allocate(old, alloc_lo, alloc_hi)
                }}
            }
            (
                reuse_or_alloc_slot!(get_m_ptr_slot,  out_slot),
                reuse_or_alloc_slot!(get_i1_ptr_slot, out_slot),
                reuse_or_alloc_slot!(get_d1_ptr_slot, out_slot),
            )
        } else {
            (
                self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi),
                self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi),
                self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi),
            )
        };

        // Compute kernel using raw pointers for multi-borrow
        unsafe {
            let null_wf = &self.wf_components.wavefront_null as *const _;
            let wf_m_misms = if !ptr_m_misms.is_null() { ptr_m_misms as *const _ } else { null_wf };
            let wf_m_open  = if !ptr_m_open.is_null()  { ptr_m_open  as *const _ } else { null_wf };
            let wf_i1_ext  = if !ptr_i1_ext.is_null()  { ptr_i1_ext  as *const _ } else { null_wf };
            let wf_d1_ext  = if !ptr_d1_ext.is_null()  { ptr_d1_ext  as *const _ } else { null_wf };

            (*m_curr_ptr).set_limits(lo, hi);
            (*i1_curr_ptr).set_limits(lo, hi);
            (*d1_curr_ptr).set_limits(lo, hi);

            if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
                piggyback::compute_affine_idm(
                    plen,
                    tlen,
                    &*wf_m_misms,
                    &*wf_m_open,
                    &*wf_i1_ext,
                    &*wf_d1_ext,
                    &mut *m_curr_ptr,
                    &mut *i1_curr_ptr,
                    &mut *d1_curr_ptr,
                    lo,
                    hi,
                );
                let low = self.memory_mode == MemoryMode::Low;
                piggyback::offload_affine(
                    &*wf_m_misms,
                    &*wf_m_open,
                    &*wf_i1_ext,
                    &*wf_d1_ext,
                    &mut *m_curr_ptr,
                    &mut *i1_curr_ptr,
                    &mut *d1_curr_ptr,
                    lo,
                    hi,
                    low,
                    buf,
                );
            } else {
                affine::compute_affine_idm(
                    plen,
                    tlen,
                    &*wf_m_misms,
                    &*wf_m_open,
                    &*wf_i1_ext,
                    &*wf_d1_ext,
                    &mut *m_curr_ptr,
                    &mut *i1_curr_ptr,
                    &mut *d1_curr_ptr,
                    lo,
                    hi,
                );
            }
        }

        // Store in components using pre-computed slot
        self.wf_components.set_m_ptr_slot(out_slot, m_curr_ptr);
        self.wf_components.set_i1_ptr_slot(out_slot, i1_curr_ptr);
        self.wf_components.set_d1_ptr_slot(out_slot, d1_curr_ptr);

        // Trim all output wavefronts independently, as C's `process_ends` does.
        // I/D may keep valid diagonals outside M's trimmed range.
        if let Some(k) = self.endsfree_begin_k(score) {
            self.endsfree_init(unsafe { &mut *m_curr_ptr }, k);
        }
        unsafe {
            compute::trim_ends(plen, tlen, &mut *m_curr_ptr);
            compute::trim_ends(plen, tlen, &mut *i1_curr_ptr);
            compute::trim_ends(plen, tlen, &mut *d1_curr_ptr);
        }
    }

    /// Compute the next wavefront for gap-affine 2-piece distance.
    ///
    /// Matches C's conditional allocation and dispatcher: when all O2/E2
    /// inputs are null, delegates to the affine (3-component) kernel and
    /// skips I2/D2 allocation entirely. For I1/D1, allocation is also
    /// skipped when both inputs are null (using victim wavefront for output).
    fn compute_affine2p_step(&mut self, score: i32) {
        let score_curr = score as usize;
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;

        let x = self.penalties.mismatch;
        let o1_plus_e1 = self.penalties.gap_opening1 + self.penalties.gap_extension1;
        let e1 = self.penalties.gap_extension1;
        let o2_plus_e2 = self.penalties.gap_opening2 + self.penalties.gap_extension2;
        let e2 = self.penalties.gap_extension2;

        // Pre-compute output slot once (avoids redundant modular divisions).
        let out_slot = self.wf_components.slot_of(score_curr);

        // Fetch input wavefront pointers, sharing slot computation for same-score lookups.
        let ptr_m_misms = if score - x >= 0 { self.wf_components.get_m_ptr((score - x) as usize) } else { WF_PTR_NONE };
        let ptr_m_open1 = if score - o1_plus_e1 >= 0 { self.wf_components.get_m_ptr((score - o1_plus_e1) as usize) } else { WF_PTR_NONE };
        let ptr_m_open2 = if score - o2_plus_e2 >= 0 { self.wf_components.get_m_ptr((score - o2_plus_e2) as usize) } else { WF_PTR_NONE };
        let (ptr_i1_ext, ptr_d1_ext) = if score - e1 >= 0 {
            let s = self.wf_components.slot_of((score - e1) as usize);
            (self.wf_components.get_i1_ptr_slot(s), self.wf_components.get_d1_ptr_slot(s))
        } else {
            (WF_PTR_NONE, WF_PTR_NONE)
        };
        let (ptr_i2_ext, ptr_d2_ext) = if score - e2 >= 0 {
            let s = self.wf_components.slot_of((score - e2) as usize);
            (self.wf_components.get_i2_ptr_slot(s), self.wf_components.get_d2_ptr_slot(s))
        } else {
            (WF_PTR_NONE, WF_PTR_NONE)
        };

        // All inputs missing or trimmed-empty: null output (C's
        // `wavefront_compute_allocate_output_null`), counted as a null step by extend.
        if is_null_wf(ptr_m_misms)
            && is_null_wf(ptr_m_open1)
            && is_null_wf(ptr_m_open2)
            && is_null_wf(ptr_i1_ext)
            && is_null_wf(ptr_d1_ext)
            && is_null_wf(ptr_i2_ext)
            && is_null_wf(ptr_d2_ext)
        {
            self.num_null_steps += 1;
            self.free_output_wavefronts_affine2p(out_slot);
            if let Some(k) = self.endsfree_begin_k(score) {
                let wf = self.endsfree_null_m(k);
                self.wf_components.set_m_ptr_slot(out_slot, wf);
            }
            return;
        }
        self.num_null_steps = 0;

        // Trimmed-empty inputs are read as the null wavefront, like missing ones.
        let ptr_m_misms = live_or_none(ptr_m_misms);
        let ptr_m_open1 = live_or_none(ptr_m_open1);
        let ptr_m_open2 = live_or_none(ptr_m_open2);
        let ptr_i1_ext = live_or_none(ptr_i1_ext);
        let ptr_d1_ext = live_or_none(ptr_d1_ext);
        let ptr_i2_ext = live_or_none(ptr_i2_ext);
        let ptr_d2_ext = live_or_none(ptr_d2_ext);

        // When all O2/E2 inputs are null, dispatch to the affine kernel (3 components).
        let i2_active = !ptr_m_open2.is_null() || !ptr_i2_ext.is_null();
        let d2_active = !ptr_m_open2.is_null() || !ptr_d2_ext.is_null();
        let use_affine_kernel = !i2_active && !d2_active;

        // Effective lo/hi over all inputs. As in C's `wavefront_compute_limits_input`,
        // null inputs contribute the null wavefront's bounds (lo = 1, hi = -1).
        let (misms_lo, misms_hi) = null_wf_bounds(ptr_m_misms);
        let (open1_lo, open1_hi) = null_wf_bounds(ptr_m_open1);
        let (open2_lo, open2_hi) = null_wf_bounds(ptr_m_open2);
        let (i1_lo, i1_hi) = null_wf_bounds(ptr_i1_ext);
        let (d1_lo, d1_hi) = null_wf_bounds(ptr_d1_ext);
        let (i2_lo, i2_hi) = null_wf_bounds(ptr_i2_ext);
        let (d2_lo, d2_hi) = null_wf_bounds(ptr_d2_ext);
        let mut lo = misms_lo
            .min(open1_lo - 1)
            .min(open2_lo - 1)
            .min(i1_lo + 1)
            .min(i2_lo + 1)
            .min(d1_lo - 1)
            .min(d2_lo - 1);
        let mut hi = misms_hi
            .max(open1_hi + 1)
            .max(open2_hi + 1)
            .max(i1_hi + 1)
            .max(i2_hi + 1)
            .max(d1_hi - 1)
            .max(d2_hi - 1);
        (lo, hi) = self.clamp_to_band(lo, hi);

        if lo > hi {
            if self.wf_components.memory_modular {
                self.free_output_wavefronts_affine2p(out_slot);
            }
            return;
        }

        // Clamp lo/hi to leave room for kernel's k±1 reads
        let hist_lo = self.wf_components.historic_min_lo;
        let hist_hi = self.wf_components.historic_max_hi;
        lo = lo.max(hist_lo + 1);
        hi = hi.min(hist_hi - 1);

        if lo > hi {
            if self.wf_components.memory_modular {
                self.free_output_wavefronts_affine2p(out_slot);
            }
            return;
        }

        // Allocation range (also resizes the null/victim wavefronts)
        let (alloc_lo, alloc_hi) = self.output_alloc_range(score, lo, hi);

        // Prefetch input wavefront offset arrays into L1 cache.
        // These prefetches overlap with the init_ends work below, hiding L1 miss latency.
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            macro_rules! prefetch_wf {
                ($ptr:expr, $lo:expr) => {
                    if !$ptr.is_null() {
                        let base = (*$ptr).centered_offsets.offset($lo as isize) as *const i8;
                        _mm_prefetch(base, _MM_HINT_T0);
                        _mm_prefetch(base.add(64), _MM_HINT_T0);
                    }
                }
            }
            prefetch_wf!(ptr_m_misms, lo);
            prefetch_wf!(ptr_m_open1, lo);
            prefetch_wf!(ptr_m_open2, lo);
            prefetch_wf!(ptr_i1_ext, lo);
            prefetch_wf!(ptr_i2_ext, lo);
            prefetch_wf!(ptr_d1_ext, lo);
            prefetch_wf!(ptr_d2_ext, lo);
        }

        // Initialize input wavefront ends using pointers (no slab re-lookup)
        unsafe {
            if !ptr_m_misms.is_null() { (*ptr_m_misms).init_ends_higher(hi);     (*ptr_m_misms).init_ends_lower(lo);     }
            if !ptr_m_open1.is_null() { (*ptr_m_open1).init_ends_higher(hi + 1); (*ptr_m_open1).init_ends_lower(lo - 1); }
            if !ptr_m_open2.is_null() { (*ptr_m_open2).init_ends_higher(hi + 1); (*ptr_m_open2).init_ends_lower(lo - 1); }
            if !ptr_i1_ext.is_null()  { (*ptr_i1_ext).init_ends_higher(hi);      (*ptr_i1_ext).init_ends_lower(lo - 1);  }
            if !ptr_i2_ext.is_null()  { (*ptr_i2_ext).init_ends_higher(hi);      (*ptr_i2_ext).init_ends_lower(lo - 1);  }
            if !ptr_d1_ext.is_null()  { (*ptr_d1_ext).init_ends_higher(hi + 1);  (*ptr_d1_ext).init_ends_lower(lo);      }
            if !ptr_d2_ext.is_null()  { (*ptr_d2_ext).init_ends_higher(hi + 1);  (*ptr_d2_ext).init_ends_lower(lo);      }
        }

        // Allocate or reuse output wavefronts using pre-computed slot.
        macro_rules! reuse_or_alloc_slot {
            ($comp:ident, $slot:expr) => {{
                let old = self.wf_components.$comp($slot);
                self.reuse_or_allocate(old, alloc_lo, alloc_hi)
            }}
        }

        let victim_ptr = &raw mut self.wf_components.wavefront_victim;

        let (m_curr_ptr, i1_curr_ptr, d1_curr_ptr) =
            if self.wf_components.memory_modular {
                (
                    reuse_or_alloc_slot!(get_m_ptr_slot,  out_slot),
                    reuse_or_alloc_slot!(get_i1_ptr_slot, out_slot),
                    reuse_or_alloc_slot!(get_d1_ptr_slot, out_slot),
                )
            } else {
                (
                    self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi),
                    self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi),
                    self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi),
                )
            };

        let (i2_curr_ptr, i2_stored) = if !use_affine_kernel {
            let p = if self.wf_components.memory_modular {
                reuse_or_alloc_slot!(get_i2_ptr_slot, out_slot)
            } else {
                self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi)
            };
            (p, true)
        } else {
            (victim_ptr, false)
        };

        let (d2_curr_ptr, d2_stored) = if !use_affine_kernel {
            let p = if self.wf_components.memory_modular {
                reuse_or_alloc_slot!(get_d2_ptr_slot, out_slot)
            } else {
                self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi)
            };
            (p, true)
        } else {
            (victim_ptr, false)
        };

        // Compute kernel using raw pointers for multi-borrow.
        // When all O2/E2 inputs are null, dispatch to the affine (3-component) kernel
        // which avoids loading/storing I2/D2 entirely (matching C's dispatcher).
        unsafe {
            let null_wf = &self.wf_components.wavefront_null as *const _;
            let wf_m_misms = if !ptr_m_misms.is_null() { ptr_m_misms as *const _ } else { null_wf };
            let wf_m_open1 = if !ptr_m_open1.is_null() { ptr_m_open1 as *const _ } else { null_wf };

            // set_limits on output wavefronts
            (*m_curr_ptr).set_limits(lo, hi);
            (*i1_curr_ptr).set_limits(lo, hi);
            (*d1_curr_ptr).set_limits(lo, hi);

            if use_affine_kernel {
                // Delegate to affine kernel (3 components) — no I2/D2 loads or stores
                let wf_i1_ext_ref = if !ptr_i1_ext.is_null() { ptr_i1_ext as *const _ } else { null_wf };
                let wf_d1_ext_ref = if !ptr_d1_ext.is_null() { ptr_d1_ext as *const _ } else { null_wf };

                if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
                    piggyback::compute_affine_idm(
                        plen,
                        tlen,
                        &*wf_m_misms,
                        &*wf_m_open1,
                        &*wf_i1_ext_ref,
                        &*wf_d1_ext_ref,
                        &mut *m_curr_ptr,
                        &mut *i1_curr_ptr,
                        &mut *d1_curr_ptr,
                        lo,
                        hi,
                    );
                    // I2/D2 have no inputs, so their occupancy (C) is 0: no offload.
                    let low = self.memory_mode == MemoryMode::Low;
                    piggyback::offload_affine(
                        &*wf_m_misms,
                        &*wf_m_open1,
                        &*wf_i1_ext_ref,
                        &*wf_d1_ext_ref,
                        &mut *m_curr_ptr,
                        &mut *i1_curr_ptr,
                        &mut *d1_curr_ptr,
                        lo,
                        hi,
                        low,
                        buf,
                    );
                } else {
                    affine::compute_affine_idm(
                        plen,
                        tlen,
                        &*wf_m_misms,
                        &*wf_m_open1,
                        &*wf_i1_ext_ref,
                        &*wf_d1_ext_ref,
                        &mut *m_curr_ptr,
                        &mut *i1_curr_ptr,
                        &mut *d1_curr_ptr,
                        lo,
                        hi,
                    );
                }
            } else {
                // Full affine2p kernel (5 components)
                let wf_m_open2 = if !ptr_m_open2.is_null() { ptr_m_open2 as *const _ } else { null_wf };
                let wf_i1_ext_ref = if !ptr_i1_ext.is_null()  { ptr_i1_ext  as *const _ } else { null_wf };
                let wf_i2_ext_ref = if !ptr_i2_ext.is_null()  { ptr_i2_ext  as *const _ } else { null_wf };
                let wf_d1_ext_ref = if !ptr_d1_ext.is_null()  { ptr_d1_ext  as *const _ } else { null_wf };
                let wf_d2_ext_ref = if !ptr_d2_ext.is_null()  { ptr_d2_ext  as *const _ } else { null_wf };

                if i2_stored { (*i2_curr_ptr).set_limits(lo, hi); }
                if d2_stored { (*d2_curr_ptr).set_limits(lo, hi); }

                if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
                    piggyback::compute_affine2p_idm(
                        plen,
                        tlen,
                        &*wf_m_misms,
                        &*wf_m_open1,
                        &*wf_m_open2,
                        &*wf_i1_ext_ref,
                        &*wf_i2_ext_ref,
                        &*wf_d1_ext_ref,
                        &*wf_d2_ext_ref,
                        &mut *m_curr_ptr,
                        &mut *i1_curr_ptr,
                        &mut *i2_curr_ptr,
                        &mut *d1_curr_ptr,
                        &mut *d2_curr_ptr,
                        lo,
                        hi,
                    );
                    let low = self.memory_mode == MemoryMode::Low;
                    piggyback::offload_affine2p(
                        &*wf_m_misms,
                        &*wf_m_open1,
                        &*wf_m_open2,
                        &*wf_i1_ext_ref,
                        &*wf_i2_ext_ref,
                        &*wf_d1_ext_ref,
                        &*wf_d2_ext_ref,
                        &mut *m_curr_ptr,
                        &mut *i1_curr_ptr,
                        &mut *i2_curr_ptr,
                        &mut *d1_curr_ptr,
                        &mut *d2_curr_ptr,
                        lo,
                        hi,
                        low,
                        buf,
                    );
                } else {
                    affine2p::compute_affine2p_idm(
                        plen,
                        tlen,
                        &*wf_m_misms,
                        &*wf_m_open1,
                        &*wf_m_open2,
                        &*wf_i1_ext_ref,
                        &*wf_i2_ext_ref,
                        &*wf_d1_ext_ref,
                        &*wf_d2_ext_ref,
                        &mut *m_curr_ptr,
                        &mut *i1_curr_ptr,
                        &mut *i2_curr_ptr,
                        &mut *d1_curr_ptr,
                        &mut *d2_curr_ptr,
                        lo,
                        hi,
                    );
                }
            }
        }

        // Store in components using pre-computed slot — NULL for skipped I2/D2
        // wavefronts, releasing the modular slot's previous ones (else they
        // leak, and the slab outgrows its reserve and moves the wavefronts).
        if self.wf_components.memory_modular {
            if !i2_stored {
                self.wavefront_slab
                    .free_ptr(self.wf_components.get_i2_ptr_slot(out_slot));
            }
            if !d2_stored {
                self.wavefront_slab
                    .free_ptr(self.wf_components.get_d2_ptr_slot(out_slot));
            }
        }
        self.wf_components.set_m_ptr_slot(out_slot, m_curr_ptr);
        self.wf_components.set_i1_ptr_slot(out_slot, i1_curr_ptr);
        self.wf_components.set_i2_ptr_slot(out_slot, if i2_stored { i2_curr_ptr } else { WF_PTR_NONE });
        self.wf_components.set_d1_ptr_slot(out_slot, d1_curr_ptr);
        self.wf_components.set_d2_ptr_slot(out_slot, if d2_stored { d2_curr_ptr } else { WF_PTR_NONE });

        // Trim all output wavefronts independently, as C's `process_ends` does.
        // I/D may keep valid diagonals outside M's trimmed range.
        if let Some(k) = self.endsfree_begin_k(score) {
            self.endsfree_init(unsafe { &mut *m_curr_ptr }, k);
        }
        unsafe {
            compute::trim_ends(plen, tlen, &mut *m_curr_ptr);
            compute::trim_ends(plen, tlen, &mut *i1_curr_ptr);
            compute::trim_ends(plen, tlen, &mut *d1_curr_ptr);
            if i2_stored {
                compute::trim_ends(plen, tlen, &mut *i2_curr_ptr);
            }
            if d2_stored {
                compute::trim_ends(plen, tlen, &mut *d2_curr_ptr);
            }
        }
    }

    // =========================================================================
    // BiWFA: Bidirectional Wavefront Alignment
    // =========================================================================

    /// Bidirectional wavefront alignment (BiWFA).
    ///
    /// Uses divide-and-conquer with O(s) memory via bidirectional
    /// wavefront expansion and breakpoint detection.
    pub fn align_biwfa(&mut self, pattern: &[u8], text: &[u8]) -> i32 {
        self.sequences.init_ascii(pattern, text, false);
        self.alignment_scope = AlignmentScope::ComputeAlignment;

        let plen = pattern.len() as i32;
        let tlen = text.len() as i32;

        self.cigar.resize(cigar_buffer_len(plen, tlen));
        self.status = AlignStatus::Ok;

        // Take sub-aligners out of self so we can borrow self.cigar/penalties separately.
        let (mut fwd, mut rev, mut base) = self.take_biwfa_subs();

        let result = Self::bialign_recursive(
            &mut self.cigar,
            &self.penalties,
            &self.heuristic,
            self.max_alignment_steps,
            &mut fwd,
            &mut rev,
            &mut base,
            pattern,
            text,
            0,
            plen,
            0,
            tlen,
            ComponentType::M,
            ComponentType::M,
            // As C, short sequences go straight to unidirectional WFA; the
            // length check applies only at the top level.
            if plen.max(tlen) <= BIALIGN_FALLBACK_MIN_LENGTH {
                0
            } else {
                i32::MAX
            },
            0,
        );

        // C's `wavefront_bialign` status mapping; on failure the CIGAR holds
        // whatever was appended before it.
        let score = match result {
            Ok(score) => {
                self.cigar.end_v = plen;
                self.cigar.end_h = tlen;
                self.status = AlignStatus::Completed;
                score
            }
            Err(status) => {
                // C reports the other failures as UNATTAINABLE.
                self.status = if status == BiStatus::MaxSteps {
                    AlignStatus::MaxStepsReached
                } else {
                    AlignStatus::EndUnreachable
                };
                0
            }
        };

        // Return sub-aligners to self for reuse on the next call.
        self.biwfa_fwd = Some(fwd);
        self.biwfa_rev = Some(rev);
        self.biwfa_base = Some(base);

        score
    }

    /// BiWFA score only (C's `wavefront_bialign_compute_score`): the score of
    /// the first breakpoint, or of the end if one direction reaches it first.
    /// Without heuristics this equals [`Self::align_end2end`]'s score; with
    /// them, the forward and reverse searches are cut off independently. On
    /// failure the status is `EndUnreachable` or `MaxStepsReached` and the
    /// returned score is meaningless.
    pub fn align_biwfa_score(&mut self, pattern: &[u8], text: &[u8]) -> i32 {
        self.sequences.init_ascii(pattern, text, false);
        self.cigar.clear();
        let (mut fwd, mut rev, base) = self.take_biwfa_subs();

        // Top-level band shift: the sub-problem is the whole alignment.
        fwd.heuristic = self.heuristic.clone();
        rev.heuristic = self.heuristic.clone();
        if let Some((min_k, max_k)) = heuristic_band(&self.heuristic) {
            let diag = text.len() as i32 - pattern.len() as i32;
            rev.heuristic = with_band(&self.heuristic, diag - max_k, diag - min_k);
        }

        let mut breakpoint = BiAlignBreakpoint::new();
        let status = Self::bialign_find_breakpoint(
            &self.penalties,
            self.max_alignment_steps,
            &mut fwd,
            &mut rev,
            pattern,
            text,
            ComponentType::M,
            ComponentType::M,
            &mut breakpoint,
        );
        let score = match status {
            BiStatus::Ok => breakpoint.score,
            BiStatus::EndReached(score) => score,
            _ => 0,
        };
        self.status = match status {
            BiStatus::Ok | BiStatus::EndReached(_) => AlignStatus::Completed,
            BiStatus::MaxSteps => AlignStatus::MaxStepsReached,
            _ => AlignStatus::EndUnreachable,
        };

        self.biwfa_fwd = Some(fwd);
        self.biwfa_rev = Some(rev);
        self.biwfa_base = Some(base);
        score
    }

    /// Take the BiWFA sub-aligners (forward, reverse, base) out of `self`,
    /// creating them on first use; callers put them back. fwd/rev use modular
    /// storage; base uses full CIGAR storage.
    #[allow(clippy::type_complexity)]
    fn take_biwfa_subs(&mut self) -> (Box<Self>, Box<Self>, Box<Self>) {
        if self.biwfa_fwd.is_none() {
            self.biwfa_fwd = Some(Box::new(Self::new_biwfa_sub(self.penalties.clone())));
            self.biwfa_rev = Some(Box::new(Self::new_biwfa_sub(self.penalties.clone())));
            let mut b = WavefrontAligner::<N>::new(self.penalties.clone());
            b.alignment_scope = AlignmentScope::ComputeAlignment;
            self.biwfa_base = Some(Box::new(b));
        }
        (
            self.biwfa_fwd.take().unwrap(),
            self.biwfa_rev.take().unwrap(),
            self.biwfa_base.take().unwrap(),
        )
    }

    /// End-to-end alignment using a custom match function (Lambda mode).
    pub fn align_lambda(
        &mut self,
        match_funct: Box<dyn Fn(i32, i32) -> bool>,
        pattern_length: i32,
        text_length: i32,
    ) -> i32 {
        self.sequences
            .init_lambda(match_funct, pattern_length, text_length, false);
        self.align_end2end_inner()
    }

    /// End-to-end alignment using packed 2-bit DNA sequences.
    pub fn align_packed2bits(
        &mut self,
        pattern: &[u8],
        pattern_length: usize,
        text: &[u8],
        text_length: usize,
    ) -> i32 {
        self.sequences
            .init_packed2bits(pattern, pattern_length, text, text_length, false);
        self.align_end2end_inner()
    }

    /// Create a sub-aligner for BiWFA with modular wavefront storage (score-only).
    fn new_biwfa_sub(penalties: WavefrontPenalties) -> Self {
        let wf_components = WavefrontComponents::<N>::new(10, 10, &penalties, true, false);
        Self {
            sequences: WavefrontSequences::new(),
            wf_components,
            wavefront_slab: WavefrontSlab::new(10, false, SlabMode::Reuse),
            penalties,
            alignment_scope: AlignmentScope::ComputeScore,
            alignment_span: AlignmentSpan::End2End,
            cigar: Cigar::new(0),
            alignment_end_pos: WavefrontPos::default(),
            status: AlignStatus::Ok,
            num_null_steps: 0,
            pattern_begin_free: 0,
            pattern_end_free: 0,
            text_begin_free: 0,
            text_end_free: 0,
            extension: false,
            compat: Compat::C,
            endsfree_best: None,
            heuristic: HeuristicStrategy::None,
            heuristic_state: HeuristicState::default(),
            memory_mode: MemoryMode::High,
            max_memory_compact: MAX_MEMORY_COMPACT,
            max_alignment_steps: i32::MAX,
            biwfa_fwd: None,
            biwfa_rev: None,
            biwfa_base: None,
        }
    }

    /// Recursive BiWFA divide-and-conquer alignment (C's
    /// `wavefront_bialign_alignment`). `heuristic` is the top-level aligner's;
    /// sub-aligners inherit it with the band shifted to each sub-problem
    /// (Marco-Sola-Lab/WFA2-lib#125). Returns the alignment score, or the status that
    /// stopped it.
    #[allow(clippy::too_many_arguments, clippy::only_used_in_recursion)]
    fn bialign_recursive(
        cigar: &mut Cigar,
        penalties: &WavefrontPenalties,
        heuristic: &HeuristicStrategy,
        max_steps: i32,
        fwd: &mut WavefrontAligner<N>,
        rev: &mut WavefrontAligner<N>,
        base: &mut WavefrontAligner<N>,
        pattern: &[u8],
        text: &[u8],
        pb: i32,
        pe: i32,
        tb: i32,
        te: i32,
        component_begin: ComponentType,
        component_end: ComponentType,
        score_remaining: i32,
        level: i32,
    ) -> Result<i32, BiStatus> {
        let plen = pe - pb;
        let tlen = te - tb;

        // Trivial cases: empty sequences
        if tlen == 0 {
            if plen > 0 {
                cigar.append_deletion(plen as usize);
            }
            return Ok(0);
        }
        if plen == 0 {
            if tlen > 0 {
                cigar.append_insertion(tlen as usize);
            }
            return Ok(0);
        }

        // Fallback to standard WFA for low-score subproblems
        let fallback = |cigar: &mut Cigar, base: &mut WavefrontAligner<N>| {
            Self::bialign_base(
                cigar,
                heuristic,
                base,
                pattern,
                text,
                pb,
                pe,
                tb,
                te,
                component_begin,
                component_end,
            )
        };
        if score_remaining <= BIALIGN_FALLBACK_MIN_SCORE {
            return fallback(cigar, base);
        }

        // Sub-aligners inherit the heuristic, with the band shifted to this
        // sub-problem (C's `wavefront_bialign_set_subsidiary_band`).
        fwd.heuristic = heuristic.clone();
        rev.heuristic = heuristic.clone();
        if let Some((min_k, max_k)) = heuristic_band(heuristic) {
            let shift = tb - pb;
            let (fmin, fmax) = (min_k - shift, max_k - shift);
            let diag = tlen - plen;
            fwd.heuristic = with_band(heuristic, fmin, fmax);
            rev.heuristic = with_band(heuristic, diag - fmax, diag - fmin);
        }

        // Find breakpoint via bidirectional expansion
        let p_sub = &pattern[pb as usize..pe as usize];
        let t_sub = &text[tb as usize..te as usize];
        let mut breakpoint = BiAlignBreakpoint::new();

        match Self::bialign_find_breakpoint(
            penalties,
            max_steps,
            fwd,
            rev,
            p_sub,
            t_sub,
            component_begin,
            component_end,
            &mut breakpoint,
        ) {
            BiStatus::Ok => {}
            // One direction reached the end first: fall back if cheap enough
            // (C's `wavefront_bialign_find_breakpoint_exception`).
            BiStatus::EndReached(score) if score <= BIALIGN_RECOVERY_MIN_SCORE => {
                return fallback(cigar, base);
            }
            BiStatus::EndReached(_) => return Err(BiStatus::Unreachable),
            status => return Err(status),
        }

        // Split at breakpoint (forward h,v coordinates)
        let bp_h = wavefront_h(breakpoint.k_forward, breakpoint.offset_forward);
        let bp_v = wavefront_v(breakpoint.k_forward, breakpoint.offset_forward);

        // Recurse on first half [pb..pb+bp_v, tb..tb+bp_h]
        Self::bialign_recursive(
            cigar,
            penalties,
            heuristic,
            max_steps,
            fwd,
            rev,
            base,
            pattern,
            text,
            pb,
            pb + bp_v,
            tb,
            tb + bp_h,
            component_begin,
            breakpoint.component,
            breakpoint.score_forward,
            level + 1,
        )?;

        // Recurse on second half [pb+bp_v..pe, tb+bp_h..te]
        Self::bialign_recursive(
            cigar,
            penalties,
            heuristic,
            max_steps,
            fwd,
            rev,
            base,
            pattern,
            text,
            pb + bp_v,
            pe,
            tb + bp_h,
            te,
            breakpoint.component,
            component_end,
            breakpoint.score_reverse,
            level + 1,
        )?;

        Ok(breakpoint.score)
    }

    /// Two-phase breakpoint detection for BiWFA (C's
    /// `wavefront_bialign_find_breakpoint`). The forward aligner ends in
    /// `component_end` and the reverse one in `component_begin`; either may
    /// finish first (end reached, or unreachable under a heuristic), which is
    /// reported instead of a breakpoint.
    #[allow(clippy::too_many_arguments)]
    fn bialign_find_breakpoint(
        penalties: &WavefrontPenalties,
        max_steps: i32,
        fwd: &mut WavefrontAligner<N>,
        rev: &mut WavefrontAligner<N>,
        pattern: &[u8],
        text: &[u8],
        component_begin: ComponentType,
        component_end: ComponentType,
        breakpoint: &mut BiAlignBreakpoint,
    ) -> BiStatus {
        let plen = pattern.len() as i32;
        let tlen = text.len() as i32;

        // Init sub-aligners with correct starting components
        // Forward starts from component_begin, reverse starts from component_end
        fwd.sequences.init_ascii(pattern, text, false);
        fwd.init_alignment_component(component_begin);
        rev.sequences.init_ascii(pattern, text, true);
        rev.init_alignment_component(component_end);
        // Each direction terminates at the other end's component.
        let (fwd_end, rev_end) = (component_end, component_begin);
        let quit = |a: &WavefrontAligner<N>| {
            if a.status == AlignStatus::EndReached {
                BiStatus::EndReached(a.alignment_end_pos.score)
            } else {
                BiStatus::Unreachable
            }
        };

        let has_affine = matches!(
            penalties.distance_metric,
            DistanceMetric::GapAffine | DistanceMetric::GapAffine2p
        );
        let has_affine2p = matches!(penalties.distance_metric, DistanceMetric::GapAffine2p);
        let gap_opening = match penalties.distance_metric {
            DistanceMetric::GapAffine => penalties.gap_opening1,
            DistanceMetric::GapAffine2p => penalties.gap_opening1.max(penalties.gap_opening2),
            _ => 0,
        };

        let max_antidiag = plen + tlen - 1;
        let mut score_fwd = 0i32;
        let mut score_rev = 0i32;
        let mut last_was_forward = false;

        // Extend at score 0
        let (done, mut fwd_max_ak) = fwd.extend_biwfa_step(0, fwd_end, true);
        if done {
            return quit(fwd);
        }
        let (done, mut rev_max_ak) = rev.extend_biwfa_step(0, rev_end, true);
        if done {
            return quit(rev);
        }

        // Phase 1: Alternate compute+extend until antidiagonals meet
        while fwd_max_ak + rev_max_ak < max_antidiag {
            // Forward
            score_fwd += 1;
            fwd.compute_step(score_fwd);
            let (done, ak) = fwd.extend_biwfa_step(score_fwd, fwd_end, true);
            fwd_max_ak = fwd_max_ak.max(ak);
            last_was_forward = true;
            if done {
                return quit(fwd);
            }
            if fwd_max_ak + rev_max_ak >= max_antidiag {
                break;
            }
            // Reverse
            score_rev += 1;
            rev.compute_step(score_rev);
            let (done, ak) = rev.extend_biwfa_step(score_rev, rev_end, true);
            rev_max_ak = rev_max_ak.max(ak);
            last_was_forward = false;
            if done {
                return quit(rev);
            }
            if score_fwd + score_rev >= max_steps {
                return BiStatus::MaxSteps;
            }
        }

        // Phase 2: Alternate compute+extend with overlap checks
        // Match C reference: overlap check first, then compute opposite direction
        let max_score_scope = fwd.wf_components.max_score_scope as i32;
        loop {
            if last_was_forward {
                // Check forward overlap
                let min_score_rev = if score_rev > max_score_scope - 1 {
                    score_rev - (max_score_scope - 1)
                } else {
                    0
                };
                if score_fwd + min_score_rev - gap_opening >= breakpoint.score {
                    break;
                }
                bialign::bialign_overlap(
                    &fwd.wf_components,
                    &fwd.wavefront_slab,
                    &rev.wf_components,
                    &rev.wavefront_slab,
                    score_fwd,
                    score_rev,
                    plen,
                    tlen,
                    penalties.gap_opening1,
                    penalties.gap_opening2,
                    has_affine,
                    has_affine2p,
                    true,
                    breakpoint,
                );
                // Compute+extend reverse
                score_rev += 1;
                rev.compute_step(score_rev);
                if rev.extend_biwfa_step(score_rev, rev_end, false).0 {
                    return quit(rev);
                }
            }

            // Check reverse overlap
            let min_score_fwd = if score_fwd > max_score_scope - 1 {
                score_fwd - (max_score_scope - 1)
            } else {
                0
            };
            if min_score_fwd + score_rev - gap_opening >= breakpoint.score {
                break;
            }
            bialign::bialign_overlap(
                &rev.wf_components,
                &rev.wavefront_slab,
                &fwd.wf_components,
                &fwd.wavefront_slab,
                score_rev,
                score_fwd,
                plen,
                tlen,
                penalties.gap_opening1,
                penalties.gap_opening2,
                has_affine,
                has_affine2p,
                false,
                breakpoint,
            );

            // Compute+extend forward
            score_fwd += 1;
            fwd.compute_step(score_fwd);
            if fwd.extend_biwfa_step(score_fwd, fwd_end, false).0 {
                return quit(fwd);
            }
            if score_fwd + score_rev >= max_steps {
                return BiStatus::MaxSteps;
            }

            last_was_forward = true;
        }

        BiStatus::Ok
    }

    /// Base case: unidirectional WFA between the given components (C's
    /// `wavefront_bialign_base`). The base aligner takes the top-level
    /// heuristic, band shifted to the sub-problem. Its CIGAR is appended
    /// whether or not it completed (empty, or max-trimmed after a Z-drop).
    #[allow(clippy::too_many_arguments)]
    fn bialign_base(
        cigar: &mut Cigar,
        heuristic: &HeuristicStrategy,
        base: &mut WavefrontAligner<N>,
        pattern: &[u8],
        text: &[u8],
        pb: i32,
        pe: i32,
        tb: i32,
        te: i32,
        component_begin: ComponentType,
        component_end: ComponentType,
    ) -> Result<i32, BiStatus> {
        let p_sub = &pattern[pb as usize..pe as usize];
        let t_sub = &text[tb as usize..te as usize];

        base.heuristic = match heuristic_band(heuristic) {
            Some((min_k, max_k)) => with_band(heuristic, min_k - (tb - pb), max_k - (tb - pb)),
            None => heuristic.clone(),
        };

        // Init sequences and alignment with the specified starting component
        base.sequences.init_ascii(p_sub, t_sub, false);
        base.init_alignment_component(component_begin);

        // C's `wavefront_unialign`, terminating on `component_end`.
        let mut score = 0;
        loop {
            if base.extend_biwfa_step(score, component_end, false).0 {
                break;
            }
            score += 1;
            if !base.wf_components.has_slot(score as usize) {
                base.status = AlignStatus::MaxStepsReached;
                break;
            }
            base.compute_step(score);
            if score >= base.max_alignment_steps {
                base.status = AlignStatus::MaxStepsReached;
                break;
            }
        }

        // Record how many flat slots were used (for efficient clear on next call).
        base.wf_components.prev_score_count = (score + 1) as usize;

        // C's `wavefront_unialign_terminate`: backtrace when the end was
        // reached, or from the best cell after a Z-drop (then max-trimmed).
        let pos = base.alignment_end_pos;
        let result = match base.status {
            AlignStatus::EndReached => {
                base.backtrace_components(score, pos.k, pos.offset, component_begin, component_end);
                Ok(score)
            }
            AlignStatus::EndUnreachable if pos.offset != OFFSET_NULL => {
                base.backtrace_components(score, pos.k, pos.offset, component_begin, component_end);
                base.maxtrim_cigar();
                Err(BiStatus::Unattainable)
            }
            _ => {
                base.cigar.clear();
                Err(BiStatus::Unattainable)
            }
        };
        cigar.append_forward(&base.cigar);
        result
    }

    /// Initialize alignment with a specific starting component (for BiWFA base case).
    fn init_alignment_component(&mut self, component: ComponentType) {
        let plen = self.sequences.pattern_length;
        let tlen = self.sequences.text_length;

        if self.alignment_scope == AlignmentScope::ComputeAlignment {
            self.cigar.resize(cigar_buffer_len(plen, tlen));
        }

        self.wf_components.resize(plen, tlen, &self.penalties);
        let modular = self.wf_components.memory_modular;
        self.wavefront_slab.set_mode(if modular { SlabMode::Reuse } else { SlabMode::Tight });
        if modular {
            self.wavefront_slab.clear();
        } else {
            self.wavefront_slab.reap_cigar();
        }

        // Pre-size slab to avoid Vec reallocations that would invalidate stored *mut Wavefront ptrs.
        {
            let num_types = match self.penalties.distance_metric {
                DistanceMetric::Edit | DistanceMetric::Indel => 1,
                DistanceMetric::GapLinear => 1,
                DistanceMetric::GapAffine => 3,
                DistanceMetric::GapAffine2p => 5,
            };
            let capacity = if self.wf_components.memory_modular {
                self.wf_components.max_score_scope * num_types + num_types
            } else {
                self.wf_components.num_wavefronts * num_types
            };
            self.wavefront_slab.reserve(capacity);
        }

        self.num_null_steps = 0;
        self.status = AlignStatus::Ok;
        self.alignment_end_pos = WavefrontPos::default();
        self.heuristic_state.reset(&self.heuristic);

        let max_score_scope = self.wf_components.max_score_scope as i32;
        let eff_lo = -(plen + max_score_scope + 2);
        let eff_hi = tlen + max_score_scope + 2;

        self.wf_components.historic_min_lo = eff_lo;
        self.wf_components.historic_max_hi = eff_hi;

        // Allocate initial wavefront at score 0
        let (alloc_lo, alloc_hi) = self.output_alloc_range(0, 0, 0);
        let wf_ptr = self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi);
        unsafe {
            (*wf_ptr).set_offset(0, 0);
            (*wf_ptr).set_limits(0, 0);
        }

        // Set the wavefront for the specified starting component
        match component {
            ComponentType::M => self.wf_components.set_m_ptr(0, wf_ptr),
            ComponentType::I1 => self.wf_components.set_i1_ptr(0, wf_ptr),
            ComponentType::D1 => self.wf_components.set_d1_ptr(0, wf_ptr),
            ComponentType::I2 => self.wf_components.set_i2_ptr(0, wf_ptr),
            ComponentType::D2 => self.wf_components.set_d2_ptr(0, wf_ptr),
        }
    }

    /// Check if the specified wavefront component has reached the end diagonal.
    fn check_component_termination(
        base: &mut WavefrontAligner<N>,
        score: i32,
        component: ComponentType,
    ) -> bool {
        let plen = base.sequences.pattern_length;
        let tlen = base.sequences.text_length;
        let k_end = tlen - plen;

        let ptr = match component {
            ComponentType::M => base.wf_components.get_m_ptr(score as usize),
            ComponentType::I1 => base.wf_components.get_i1_ptr(score as usize),
            ComponentType::D1 => base.wf_components.get_d1_ptr(score as usize),
            ComponentType::I2 => base.wf_components.get_i2_ptr(score as usize),
            ComponentType::D2 => base.wf_components.get_d2_ptr(score as usize),
        };

        if ptr.is_null() {
            return false;
        }

        let wf = unsafe { &*ptr };
        if k_end >= wf.lo && k_end <= wf.hi && wf.get_offset(k_end) >= tlen {
            base.alignment_end_pos = WavefrontPos {
                score,
                k: k_end,
                offset: tlen,
            };
            base.status = AlignStatus::EndReached;
            true
        } else {
            false
        }
    }

    /// Clamp a computed diagonal range to the heuristic band, as C's
    /// `wavefront_compute_limits_heuristic` (Marco-Sola-Lab/WFA2-lib#125) does for
    /// banded-static and banded-adaptive in the linear/affine/affine2p steps.
    #[inline(always)]
    fn clamp_to_band(&self, lo: i32, hi: i32) -> (i32, i32) {
        match self.heuristic {
            HeuristicStrategy::BandedStatic { min_k, max_k }
            | HeuristicStrategy::BandedAdaptive { min_k, max_k, .. } => {
                (lo.max(min_k), hi.min(max_k))
            }
            _ => (lo, hi),
        }
    }

    /// C's `wavefront_compute_limits_output`, with the ends-free limits of
    /// `wavefront_compute_allocate_output`: the diagonal range to allocate an
    /// output wavefront for `lo..=hi` at `score` over. Pads by
    /// `max_score_scope + 1` and merges into the alignment's running range,
    /// so inputs read up to `max_score_scope` scores later stay covered.
    /// Also resizes the null/victim wavefronts to it.
    #[inline(always)]
    fn output_alloc_range(&mut self, score: i32, lo: i32, hi: i32) -> (i32, i32) {
        let (mut lo, mut hi) = (lo, hi);
        if self.alignment_span == AlignmentSpan::EndsFree
            && score >= 0
            && let Some(k) = self.endsfree_begin_k(score)
        {
            lo = lo.min(-k);
            hi = hi.max(k);
        }
        let pad = self.wf_components.max_score_scope as i32 + 1;
        let c = &mut self.wf_components;
        if lo - pad < c.alloc_min_lo || hi + pad > c.alloc_max_hi {
            c.alloc_min_lo = c.alloc_min_lo.min(lo - pad);
            c.alloc_max_hi = c.alloc_max_hi.max(hi + pad);
            c.resize_null_victim(c.alloc_min_lo, c.alloc_max_hi);
        }
        (c.alloc_min_lo, c.alloc_max_hi)
    }

    /// An output wavefront over `alloc_lo..=alloc_hi`: `old` (the modular
    /// slot's previous wavefront) re-initialized if it is large enough, else
    /// a slab allocation (`old` is returned to the slab).
    #[inline(always)]
    fn reuse_or_allocate(
        &mut self,
        old: *mut crate::wavefront::Wavefront,
        alloc_lo: i32,
        alloc_hi: i32,
    ) -> *mut crate::wavefront::Wavefront {
        if !old.is_null() {
            let fits = unsafe { (*old).wf_elements_allocated }
                >= crate::offset::wavefront_length(alloc_lo, alloc_hi);
            if fits {
                unsafe {
                    (*old).init(alloc_lo, alloc_hi);
                    (*old).status = WavefrontStatus::Busy;
                }
                return old;
            }
            self.wavefront_slab.free_ptr(old);
        }
        self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi)
    }

    /// C's `wavefront_compute_endsfree_required`: with match < 0, the
    /// begin-free cells (h,0) / (0,v) with h = v = d are reached at score
    /// d * -match. Returns d if `score` is such a score within the free range.
    /// (With match == 0 they are all seeded at score 0 instead.)
    fn endsfree_begin_k(&self, score: i32) -> Option<i32> {
        let m = self.penalties.match_;
        if m == 0 || self.alignment_span != AlignmentSpan::EndsFree {
            return None;
        }
        if self.text_begin_free == 0 && self.pattern_begin_free == 0 {
            return None;
        }
        if score % (-m) != 0 {
            return None;
        }
        let k = score / (-m);
        (self.text_begin_free >= k || self.pattern_begin_free >= k).then_some(k)
    }

    /// C's `wavefront_compute_endsfree_init`: add the begin-free cells on
    /// diagonals `k` (text) and `-k` (pattern) to a computed M-wavefront,
    /// before its ends are trimmed.
    fn endsfree_init(&mut self, wf: &mut crate::wavefront::Wavefront, k: i32) {
        let (lo, hi) = (wf.lo, wf.hi);
        if self.text_begin_free >= k {
            if hi >= k {
                if wf.get_offset(k) <= k {
                    self.endsfree_init_offset(wf, k, 0, k); // (h, v) = (k, 0)
                }
            } else {
                for kk in hi + 1..k {
                    wf.set_null_cell(kk);
                }
                self.endsfree_init_offset(wf, k, 0, k);
                wf.hi = k;
            }
        }
        if self.pattern_begin_free >= k {
            if lo <= -k {
                if wf.get_offset(-k) <= 0 {
                    self.endsfree_init_offset(wf, -k, k, 0); // (h, v) = (0, k)
                }
            } else {
                self.endsfree_init_offset(wf, -k, k, 0);
                for kk in -k + 1..lo {
                    wf.set_null_cell(kk);
                }
                wf.lo = -k;
            }
        }
    }

    /// C's `wavefront_compute_endsfree_init_offset`: set a begin-free cell
    /// (v, h) on diagonal `k`, starting a backtrace chain there if piggyback.
    fn endsfree_init_offset(
        &mut self,
        wf: &mut crate::wavefront::Wavefront,
        k: i32,
        v: i32,
        h: i32,
    ) {
        wf.set_offset(k, h);
        if let Some(buf) = self.wf_components.bt_buffer.as_mut() {
            set_bt_init(wf, buf, k, v, h);
            wf.mark_written(k, k);
        }
    }

    /// C's `wavefront_compute_endsfree_allocate_null`: when every input is
    /// null on a begin-free score, M holds just the begin-free cells.
    fn endsfree_null_m(&mut self, k: i32) -> *mut crate::wavefront::Wavefront {
        let text = self.text_begin_free >= k;
        let pattern = self.pattern_begin_free >= k;
        let (lo, hi) = match (text, pattern) {
            (true, true) => (-k, k),
            (true, false) => (k, k),
            (false, true) => (-k, -k),
            (false, false) => (0, 0),
        };
        let (alloc_lo, alloc_hi) = self.output_alloc_range(-1, lo, hi);
        let wf_ptr = self.wavefront_slab.allocate_ptr(alloc_lo, alloc_hi);
        let wf = unsafe { &mut *wf_ptr };
        for kk in lo..=hi {
            wf.set_null_cell(kk);
        }
        if text {
            self.endsfree_init_offset(wf, k, 0, k);
        }
        if pattern {
            self.endsfree_init_offset(wf, -k, k, 0);
        }
        wf.set_limits(lo, hi);
        wf_ptr
    }

    /// Store a null M-wavefront at `score` (edit/indel/linear), releasing any
    /// modular-slot wavefront it replaces.
    #[inline(always)]
    fn set_null_output_m(&mut self, score: usize) {
        let slot = self.wf_components.slot_of(score);
        if self.wf_components.memory_modular {
            self.wavefront_slab
                .free_ptr(self.wf_components.get_m_ptr_slot(slot));
        }
        self.wf_components.set_m_ptr_slot(slot, WF_PTR_NONE);
    }

    /// Free old wavefronts at the modular slot for affine (M, I1, D1).
    #[inline(always)]
    fn free_output_wavefronts_affine(&mut self, slot: usize) {
        self.wavefront_slab.free_ptr(self.wf_components.get_m_ptr_slot(slot));
        self.wavefront_slab.free_ptr(self.wf_components.get_i1_ptr_slot(slot));
        self.wavefront_slab.free_ptr(self.wf_components.get_d1_ptr_slot(slot));
        self.wf_components.set_m_ptr_slot(slot, WF_PTR_NONE);
        self.wf_components.set_i1_ptr_slot(slot, WF_PTR_NONE);
        self.wf_components.set_d1_ptr_slot(slot, WF_PTR_NONE);
    }

    #[inline(always)]
    fn free_output_wavefronts_affine2p(&mut self, slot: usize) {
        self.wavefront_slab.free_ptr(self.wf_components.get_m_ptr_slot(slot));
        self.wavefront_slab.free_ptr(self.wf_components.get_i1_ptr_slot(slot));
        self.wavefront_slab.free_ptr(self.wf_components.get_d1_ptr_slot(slot));
        self.wavefront_slab.free_ptr(self.wf_components.get_i2_ptr_slot(slot));
        self.wavefront_slab.free_ptr(self.wf_components.get_d2_ptr_slot(slot));
        self.wf_components.set_m_ptr_slot(slot, WF_PTR_NONE);
        self.wf_components.set_i1_ptr_slot(slot, WF_PTR_NONE);
        self.wf_components.set_d1_ptr_slot(slot, WF_PTR_NONE);
        self.wf_components.set_i2_ptr_slot(slot, WF_PTR_NONE);
        self.wf_components.set_d2_ptr_slot(slot, WF_PTR_NONE);
    }

    /// Dispatch compute step based on distance metric.
    fn compute_step(&mut self, score: i32) {
        match self.penalties.distance_metric {
            DistanceMetric::Edit | DistanceMetric::Indel => self.compute_edit_step(score),
            DistanceMetric::GapLinear => self.compute_linear_step(score),
            DistanceMetric::GapAffine => self.compute_affine_step(score),
            DistanceMetric::GapAffine2p => self.compute_affine2p_step(score),
        }
    }

    /// One BiWFA extend step (C's `wavefront_extend_end2end_max`, or
    /// `wavefront_extend_end2end` when `max` is false): extend M, then report
    /// done if the alignment is unreachable (null steps), `component_end` is
    /// reached, or the heuristic drops it. Returns `(done, max_antidiagonal)`;
    /// the antidiagonal is 0 unless `max` and not done.
    fn extend_biwfa_step(
        &mut self,
        score: i32,
        component_end: ComponentType,
        max: bool,
    ) -> (bool, i32) {
        let m_ptr = self.wf_components.get_m_ptr(score as usize);
        if m_ptr.is_null() {
            if self.num_null_steps > self.wf_components.max_score_scope as i32 {
                self.status = AlignStatus::EndUnreachable;
                return (true, 0);
            }
            return (false, 0);
        }
        let wf = unsafe { &mut *m_ptr };
        let lambda = self.sequences.mode == SequenceMode::Lambda;
        let max_ak = match (max, lambda) {
            (true, true) => extend::extend_matches_custom_end2end_max(&self.sequences, wf),
            (true, false) => extend::extend_matches_packed_end2end_max(&self.sequences, wf),
            (false, true) => {
                extend::extend_matches_custom_end2end(&self.sequences, wf);
                0
            }
            (false, false) => {
                extend::extend_matches_packed_end2end(&self.sequences, wf);
                0
            }
        };
        if Self::check_component_termination(self, score, component_end) {
            return (true, 0);
        }
        if !matches!(self.heuristic, HeuristicStrategy::None)
            && self.apply_heuristic_cutoff(score, m_ptr)
        {
            self.status = AlignStatus::EndUnreachable;
            return (true, 0);
        }
        (false, max_ak)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::penalties::WavefrontPenalties;

    #[test]
    fn test_identical_sequences() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_one_mismatch() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 1);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_one_insertion() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"ACT", b"ACGT");
        assert_eq!(score, 1);
    }

    #[test]
    fn test_one_deletion() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"ACGT", b"ACT");
        assert_eq!(score, 1);
    }

    #[test]
    fn test_multiple_edits() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        // kitten → sitting: 3 edits (k→s, e→i, +g)
        let score = aligner.align_end2end(b"kitten", b"sitting");
        assert_eq!(score, 3);
    }

    #[test]
    fn test_empty_vs_nonempty() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"", b"ACGT");
        assert_eq!(score, 4);
    }

    #[test]
    fn test_nonempty_vs_empty() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"ACGT", b"");
        assert_eq!(score, 4);
    }

    #[test]
    fn test_completely_different() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"AAAA", b"TTTT");
        assert_eq!(score, 4);
    }

    #[test]
    fn test_longer_sequences() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_end2end(b"ACGTACGTACGTACGT", b"ACGTACGTACGTACGT");
        assert_eq!(score, 0);
    }

    #[test]
    fn test_indel_distance() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_indel());
        // Indel distance: only insertions and deletions, no substitutions
        // "ACGT" vs "ACTT": need to delete G, insert T → indel dist = 2
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 2);
    }

    #[test]
    fn test_indel_identical() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_indel());
        let score = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score, 0);
    }

    #[test]
    fn test_indel_insertion() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_indel());
        let score = aligner.align_end2end(b"ACT", b"ACGT");
        assert_eq!(score, 1);
    }

    #[test]
    fn test_reuse_aligner() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score1 = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score1, 1);

        let score2 = aligner.align_end2end(b"AAA", b"BBB");
        assert_eq!(score2, 3);
    }

    #[test]
    fn test_max_steps_limit() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.set_max_alignment_steps(2);
        let score = aligner.align_end2end(b"AAAA", b"TTTT");
        assert_eq!(aligner.status(), AlignStatus::MaxStepsReached);
        // Score should be 2 (the limit), not the actual edit distance 4
        assert!(score <= 4);
    }

    /// Pattern/text whose end diagonal (k = -30) lies outside a [-5, 5] band.
    fn band_excludes_target() -> (Vec<u8>, Vec<u8>) {
        let pattern = b"ACGTTGCAAGCT".repeat(10);
        let text = pattern[..90].to_vec();
        (pattern, text)
    }

    #[test]
    fn test_banded_unreachable_does_not_overflow_components() {
        // The band keeps the end diagonal out of reach, so the score grows past
        // `num_wavefronts`; this previously wrote out of bounds.
        let (pattern, text) = band_excludes_target();
        let band = HeuristicStrategy::BandedStatic {
            min_k: -5,
            max_k: 5,
        };

        let mut edit = EditAligner::new(WavefrontPenalties::new_edit());
        edit.alignment_scope = AlignmentScope::ComputeAlignment;
        edit.set_heuristic(band.clone());
        edit.align_end2end(&pattern, &text);
        assert_ne!(edit.status(), AlignStatus::Completed);

        let mut affine2p = Affine2pAligner::new(WavefrontPenalties::new_affine2p(
            crate::penalties::Affine2pPenalties {
                match_: 0,
                mismatch: 4,
                gap_opening1: 6,
                gap_extension1: 2,
                gap_opening2: 24,
                gap_extension2: 1,
            },
        ));
        affine2p.alignment_scope = AlignmentScope::ComputeAlignment;
        affine2p.set_heuristic(band);
        affine2p.align_end2end(&pattern, &text);
        assert_ne!(affine2p.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_banded_linear_init_ends_empty_range() {
        // Band cutoffs leave lo below the requested init range; the init_ends
        // fill must then be a no-op rather than a wrapped slice length.
        let (pattern, text) = band_excludes_target();
        let mut linear = EditAligner::new(WavefrontPenalties::new_linear(
            crate::penalties::LinearPenalties {
                match_: 0,
                mismatch: 4,
                indel: 2,
            },
        ));
        linear.alignment_scope = AlignmentScope::ComputeAlignment;
        linear.set_heuristic(HeuristicStrategy::BandedStatic {
            min_k: -5,
            max_k: 5,
        });
        linear.align_end2end(&pattern, &text);
        assert_ne!(linear.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_banded_unreachable_edit_linear_end_unreachable() {
        // Once the band and trimming empty the wavefront, edit/indel/linear must
        // emit null wavefronts (C's null-input check) so extend reports
        // END_UNREACHABLE. Previously they kept computing from the empty
        // wavefront, widening lo/hi each step past the allocated range and,
        // in score-only mode, never terminating. The step cap turns a
        // regression into a failure instead of a hang.
        let (pattern, text) = band_excludes_target();
        let penalties = [
            WavefrontPenalties::new_edit(),
            WavefrontPenalties::new_indel(),
            WavefrontPenalties::new_linear(crate::penalties::LinearPenalties {
                match_: 0,
                mismatch: 4,
                indel: 2,
            }),
        ];
        for p in penalties {
            for scope in [
                AlignmentScope::ComputeAlignment,
                AlignmentScope::ComputeScore,
            ] {
                let mut aligner = EditAligner::new(p.clone());
                aligner.alignment_scope = scope;
                aligner.set_max_alignment_steps(10_000);
                aligner.set_heuristic(HeuristicStrategy::BandedStatic {
                    min_k: -5,
                    max_k: 5,
                });
                aligner.align_end2end(&pattern, &text);
                assert_eq!(
                    aligner.status(),
                    AlignStatus::EndUnreachable,
                    "{:?} {:?}",
                    p.distance_metric,
                    scope
                );

                // The aligner must still produce correct results afterwards.
                aligner.set_heuristic(HeuristicStrategy::None);
                aligner.alignment_scope = AlignmentScope::ComputeAlignment;
                aligner.align_end2end(b"ACGTACGT", b"ACGAACGT");
                assert_eq!(aligner.status(), AlignStatus::Completed);
            }
        }
    }

    fn affine_penalties() -> WavefrontPenalties {
        WavefrontPenalties::new_affine(crate::penalties::AffinePenalties {
            match_: 0,
            mismatch: 4,
            gap_opening: 6,
            gap_extension: 2,
        })
    }

    fn affine2p_penalties() -> WavefrontPenalties {
        WavefrontPenalties::new_affine2p(crate::penalties::Affine2pPenalties {
            match_: 0,
            mismatch: 4,
            gap_opening1: 6,
            gap_extension1: 2,
            gap_opening2: 24,
            gap_extension2: 1,
        })
    }

    #[test]
    fn test_banded_unreachable_affine_end_unreachable() {
        // Affine/affine2p must also emit null wavefronts once every input is
        // trimmed empty, so unreachable banded alignments end as END_UNREACHABLE
        // (score-only mode previously never terminated). The step cap turns a
        // regression into a failure instead of a hang.
        let (pattern, text) = band_excludes_target();
        let band = HeuristicStrategy::BandedStatic {
            min_k: -5,
            max_k: 5,
        };
        for scope in [
            AlignmentScope::ComputeAlignment,
            AlignmentScope::ComputeScore,
        ] {
            let mut affine = AffineAligner::new(affine_penalties());
            affine.alignment_scope = scope;
            affine.set_max_alignment_steps(10_000);
            affine.set_heuristic(band.clone());
            let score = affine.align_end2end(&pattern, &text);
            assert_eq!(affine.status(), AlignStatus::EndUnreachable, "{scope:?}");
            // Give-up score reported by C WFA2-lib (16969b5b) score-only mode.
            assert_eq!(score, 58, "{scope:?}");

            let mut affine2p = Affine2pAligner::new(affine2p_penalties());
            affine2p.alignment_scope = scope;
            affine2p.set_max_alignment_steps(10_000);
            affine2p.set_heuristic(band.clone());
            affine2p.align_end2end(&pattern, &text);
            assert_eq!(affine2p.status(), AlignStatus::EndUnreachable, "{scope:?}");
        }
    }

    #[test]
    fn test_banded_end_outside_band_is_unreachable() {
        // End diagonal k = 12 lies outside the [-10, 10] band. C WFA2-lib
        // 16969b5b still completed this pair (I/D wavefronts carried it past
        // the band: 210 affine / 204 affine2p); since Marco-Sola-Lab/WFA2-lib#125
        // clamps the computed range to the band, C (3282f83) reports it
        // unreachable, and so must we.
        let pattern = b"GCAAGACGATTAGACTTGATAAATGCAATCATCGACTGTAAACAGACGAGCTGTCCTGGGACACTACTGACCGGATCGACTACCGCGT";
        let text = b"ACAAGACACGCTCAGTGACTTGCGTCACCATCGAGCATCCGTCTTAGCAACAGGGCTGTCTGGGGCTGCCATACGACCCGCTGACCGCACGTACGGGCCT";
        let band = HeuristicStrategy::BandedStatic {
            min_k: -10,
            max_k: 10,
        };

        let mut affine = AffineAligner::new(affine_penalties());
        affine.alignment_scope = AlignmentScope::ComputeAlignment;
        affine.set_heuristic(band.clone());
        affine.align_end2end(pattern, text);
        assert_eq!(affine.status(), AlignStatus::EndUnreachable);
        assert!(affine.cigar().operations_slice().is_empty());

        let mut affine2p = Affine2pAligner::new(affine2p_penalties());
        affine2p.alignment_scope = AlignmentScope::ComputeAlignment;
        affine2p.set_heuristic(band);
        affine2p.align_end2end(pattern, text);
        assert_eq!(affine2p.status(), AlignStatus::EndUnreachable);
    }

    /// Pairs whose best ends-free alignment uses begin gaps (P0 = 5, T0 = 3).
    const BEGIN_FREE_PAIRS: [(&[u8], &[u8]); 3] = [
        (b"TTTTTACGTACGTACGT", b"ACGTACGTACGT"),
        (b"ACGTACGTACGT", b"GGGACGTACGTACGT"),
        (b"ACGTTCGTACGA", b"CCACGTACGTACGA"),
    ];

    #[test]
    fn test_endsfree_begin_free_match_zero() {
        // With match == 0 the begin-free cells are seeded at score 0. Expected
        // CIGARs and scores are from C WFA2-lib (16969b5b and 3282f83 agree);
        // previously `*_begin_free` was ignored and the result was global.
        let expected = [("5D12M", 0), ("3I12M", 0), ("2I4M1X7M", 4)];
        let mut affine = AffineAligner::new(affine_penalties());
        affine.alignment_scope = AlignmentScope::ComputeAlignment;
        for ((pattern, text), (cigar, score)) in BEGIN_FREE_PAIRS.iter().zip(expected) {
            affine.set_alignment_free_ends(5, 0, 3, 0);
            assert_eq!(affine.align_endsfree(pattern, text), score);
            assert_eq!(affine.cigar().to_string_rle(true), cigar);
            assert!(affine.cigar().check_alignment(pattern, text).is_ok());
        }

        let mut edit = EditAligner::new(WavefrontPenalties::new_edit());
        for ((pattern, text), score) in BEGIN_FREE_PAIRS.iter().zip([0, 0, 1]) {
            edit.set_alignment_free_ends(5, 0, 3, 0);
            assert_eq!(edit.align_endsfree(pattern, text), score);
        }
    }

    #[test]
    fn test_empty_sequences_cigar() {
        // HPCBio/wfa2lib-rs#41: a fresh aligner's CIGAR buffer was empty for
        // empty pattern and text, and the high-memory backtraces underflowed
        // `begin_offset` (a panic in debug builds).
        let affine = WavefrontPenalties::new_affine(crate::penalties::AffinePenalties {
            match_: 0,
            mismatch: 4,
            gap_opening: 6,
            gap_extension: 2,
        });
        for pen in [WavefrontPenalties::new_edit(), affine] {
            for mode in [MemoryMode::High, MemoryMode::Med, MemoryMode::Low] {
                for endsfree in [false, true] {
                    let mut a = WavefrontAligner::<3>::new(pen.clone());
                    a.set_memory_mode(mode);
                    a.set_alignment_scope(AlignmentScope::ComputeAlignment);
                    let score = if endsfree {
                        a.set_alignment_free_ends(0, 0, 0, 0);
                        a.align_endsfree(b"", b"")
                    } else {
                        a.align_end2end(b"", b"")
                    };
                    assert_eq!(score, 0, "{mode:?} endsfree={endsfree}");
                    assert_eq!(a.status, AlignStatus::Completed);
                    assert!(a.cigar().is_null());
                }
            }
            let mut a = WavefrontAligner::<3>::new(pen.clone());
            assert_eq!(a.align_biwfa(b"", b""), 0);
            assert!(a.cigar().is_null());
        }
    }

    #[test]
    fn test_endsfree_begin_free_seed_backtrace_med_low() {
        // match == 0: the alignment ends on a begin-free seed's diagonal at
        // score 0. Med/low read the seed's backtrace as unwritten and
        // started the CIGAR at (0, 0) (`1M12I`, invalid); C gives `2I1M10I`.
        let (pattern, text) = (b"G".as_slice(), b"TCGGGACATCGAG".as_slice());
        for mode in [MemoryMode::High, MemoryMode::Med, MemoryMode::Low] {
            let mut a = EditAligner::new(WavefrontPenalties::new_linear(
                crate::penalties::LinearPenalties {
                    match_: 0,
                    mismatch: 3,
                    indel: 3,
                },
            ));
            a.set_memory_mode(mode);
            a.set_alignment_scope(AlignmentScope::ComputeAlignment);
            a.set_alignment_free_ends(0, 0, 4, 11);
            a.align_endsfree(pattern, text);
            assert_eq!(a.cigar().to_string_rle(true), "2I1M10I", "{mode:?}");
        }
    }

    #[test]
    fn test_affine2p_modular_slots_release_i2_d2() {
        // Score-only affine2p: steps that skip I2/D2 (affine kernel) must
        // release the modular slot's old I2/D2, else the slab outgrows its
        // reserve and reallocates under the components' pointers (debug
        // builds assert; release builds could segfault).
        let penalties = WavefrontPenalties::new_affine2p(crate::penalties::Affine2pPenalties {
            match_: -3,
            mismatch: 6,
            gap_opening1: 1,
            gap_extension1: 5,
            gap_opening2: 9,
            gap_extension2: 1,
        });
        let mut a = Affine2pAligner::new(penalties);
        a.set_alignment_free_ends(0, 33, 0, 35);
        a.align_endsfree(
            b"TGGATTGAGATTCCCCGGTATGCTGACAGAAGC",
            b"ACATTCCGACATAGAGCCATATGATAAAAGTTAGT",
        );
        assert_eq!(a.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_extension_ignores_begin_free() {
        // C's extension preset forces P0 = T0 = 0, even if begin gaps were set.
        let mut a = AffineAligner::new(affine_penalties());
        a.alignment_scope = AlignmentScope::ComputeAlignment;
        a.set_alignment_free_ends(5, 0, 3, 0);
        a.set_alignment_extension();
        a.align_endsfree(b"TTTTTACGTACGTACGT", b"ACGTACGTACGT");
        assert!(!a.cigar().to_string_rle(true).starts_with("5D"));
    }

    // --- Gap-linear tests ---

    #[test]
    fn test_linear_identical() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_linear(
            crate::penalties::LinearPenalties {
                match_: 0,
                mismatch: 4,
                indel: 2,
            },
        ));
        let score = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_linear_one_mismatch() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_linear(
            crate::penalties::LinearPenalties {
                match_: 0,
                mismatch: 4,
                indel: 2,
            },
        ));
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        // Mismatch cost = 4, or ins+del = 2+2 = 4. Either way score = 4.
        assert_eq!(score, 4);
    }

    #[test]
    fn test_linear_one_insertion() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_linear(
            crate::penalties::LinearPenalties {
                match_: 0,
                mismatch: 4,
                indel: 2,
            },
        ));
        let score = aligner.align_end2end(b"ACT", b"ACGT");
        assert_eq!(score, 2);
    }

    #[test]
    fn test_linear_with_cigar() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_linear(
            crate::penalties::LinearPenalties {
                match_: 0,
                mismatch: 4,
                indel: 2,
            },
        ));
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACT", b"ACGT");
        assert_eq!(score, 2);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"ACT", b"ACGT").unwrap();
    }

    // --- Gap-affine tests ---

    #[test]
    fn test_affine_identical() {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        let score = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_affine_one_mismatch() {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 4); // One mismatch = X=4
    }

    #[test]
    fn test_affine_one_insertion() {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        let score = aligner.align_end2end(b"ACT", b"ACGT");
        assert_eq!(score, 8); // One gap = O+E = 6+2 = 8
    }

    #[test]
    fn test_affine_with_cigar() {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 4);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"ACGT", b"ACTT").unwrap();
    }

    #[test]
    fn test_affine_gap_extend() {
        // Two consecutive insertions should cost O + 2*E = 6 + 4 = 10
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"AT", b"ACGT");
        assert_eq!(score, 10); // O + 2*E = 6 + 4 = 10
        let cigar = aligner.cigar();
        cigar.check_alignment(b"AT", b"ACGT").unwrap();
    }

    #[test]
    fn test_affine_reuse() {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        let score1 = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score1, 0);
        let score2 = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score2, 4);
    }

    // --- Gap-affine 2-piece tests ---

    fn make_affine2p_aligner() -> Affine2pAligner {
        Affine2pAligner::new(WavefrontPenalties::new_affine2p(
            crate::penalties::Affine2pPenalties {
                match_: 0,
                mismatch: 4,
                gap_opening1: 6,
                gap_extension1: 2,
                gap_opening2: 24,
                gap_extension2: 1,
            },
        ))
    }

    #[test]
    fn test_affine2p_identical() {
        let mut aligner = make_affine2p_aligner();
        let score = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_affine2p_one_mismatch() {
        let mut aligner = make_affine2p_aligner();
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 4); // One mismatch = X=4
    }

    #[test]
    fn test_affine2p_one_insertion() {
        let mut aligner = make_affine2p_aligner();
        let score = aligner.align_end2end(b"ACT", b"ACGT");
        assert_eq!(score, 8); // One gap = O1+E1 = 6+2 = 8
    }

    #[test]
    fn test_affine2p_with_cigar() {
        let mut aligner = make_affine2p_aligner();
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score, 4);
        let cigar = aligner.cigar();
        cigar.check_alignment(b"ACGT", b"ACTT").unwrap();
    }

    #[test]
    fn test_affine2p_reuse() {
        let mut aligner = make_affine2p_aligner();
        let score1 = aligner.align_end2end(b"ACGT", b"ACGT");
        assert_eq!(score1, 0);
        let score2 = aligner.align_end2end(b"ACGT", b"ACTT");
        assert_eq!(score2, 4);
    }

    // --- Ends-free tests ---

    #[test]
    fn test_endsfree_text_suffix_free() {
        // Pattern "ACGT" (4) vs Text "ACGTXXXX" (8)
        // With text_end_free=4: pattern fully consumed, remaining text (4) is free
        // Edit distance: score should be 0 (perfect match of pattern against text prefix)
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.set_alignment_free_ends(0, 0, 0, 4);
        let score = aligner.align_endsfree(b"ACGT", b"ACGTXXXX");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_endsfree_pattern_suffix_free() {
        // Pattern "ACGTXXXX" (8) vs Text "ACGT" (4)
        // With pattern_end_free=4: text fully consumed, remaining pattern (4) is free
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.set_alignment_free_ends(0, 4, 0, 0);
        let score = aligner.align_endsfree(b"ACGTXXXX", b"ACGT");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    #[test]
    fn test_endsfree_with_cigar() {
        // Ends-free with CIGAR computation
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        aligner.set_alignment_free_ends(0, 0, 0, 4);
        let score = aligner.align_endsfree(b"ACGT", b"ACGTXXXX");
        assert_eq!(score, 0);
        // CIGAR should cover pattern fully but not all of text
        let cigar = aligner.cigar();
        assert!(!cigar.is_null());
    }

    #[test]
    fn test_endsfree_with_mismatch() {
        // Pattern "ACGT" vs Text "ACTTXXXX"
        // With text_end_free=4: one mismatch at position 2
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        aligner.set_alignment_free_ends(0, 0, 0, 4);
        let score = aligner.align_endsfree(b"ACGT", b"ACTTXXXX");
        assert_eq!(score, 1);
    }

    #[test]
    fn test_endsfree_affine() {
        // Affine ends-free: pattern fully consumed with text suffix free
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        aligner.set_alignment_free_ends(0, 0, 0, 4);
        let score = aligner.align_endsfree(b"ACGT", b"ACGTXXXX");
        assert_eq!(score, 0);
    }

    // --- Extension alignment tests ---

    #[test]
    fn test_extension_basic() {
        // Extension mode: find maximal-scoring prefix
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        aligner.set_alignment_extension();
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        // "ACGTXXXX" vs "ACGTTTTT": first 4 match, then mismatches
        // Maxtrim should keep the matching prefix
        let _score = aligner.align_endsfree(b"ACGTXXXX", b"ACGTTTTT");
        let cigar = aligner.cigar();
        assert!(!cigar.is_null());
        // The maximal-scoring prefix should be the 4 matches
        assert_eq!(cigar.end_v, 4);
        assert_eq!(cigar.end_h, 4);
        assert_eq!(aligner.status(), AlignStatus::Partial);
    }

    #[test]
    fn test_extension_identical() {
        // Extension with identical sequences: no trim needed
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        aligner.set_alignment_extension();
        aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let score = aligner.align_endsfree(b"ACGT", b"ACGT");
        assert_eq!(score, 0);
        assert_eq!(aligner.status(), AlignStatus::Completed);
    }

    // --- BiWFA tests ---

    #[test]
    fn test_biwfa_base_case_matches_standard() {
        // This tests the base case fallback specifically (short sequences)
        let penalties = WavefrontPenalties::new_affine(crate::penalties::AffinePenalties {
            match_: 0,
            mismatch: 4,
            gap_opening: 6,
            gap_extension: 2,
        });
        let pattern = b"GCTGTTTGCCAGAATATGGTTTGACTCGAGAGGCTAGTCTACATGGCCGGTGTCACTTAGAAACAAGCGCCGTCCCTAGATACACCTGCCCACCACGGCTCGT";
        let text = b"GCCTGGTTCAAAACTACTGTGACTGCCAGAGCCTCTCAGGCTTGAGGGCAGACTTACAAAAGCGCGGTTCGAACTATCGCGACAACACCACGGTAACGTC";

        // Standard alignment
        let mut std_aligner = AffineAligner::new(penalties.clone());
        std_aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let std_score = std_aligner.align_end2end(pattern, text);
        let std_cigar = std_aligner.cigar().to_string_rle(true);
        std_aligner.cigar().check_alignment(pattern, text).unwrap();

        // BiWFA
        let mut biwfa_aligner = AffineAligner::new(penalties);
        let biwfa_score = biwfa_aligner.align_biwfa(pattern, text);
        let biwfa_cigar = biwfa_aligner.cigar().to_string_rle(true);

        assert_eq!(biwfa_score, std_score, "Scores differ");

        if let Err(e) = biwfa_aligner.cigar().check_alignment(pattern, text) {
            panic!(
                "BiWFA CIGAR invalid: {}\nStd CIGAR: {}\nBiWFA CIGAR: {}\nStd begin/end: {}/{}, BiWFA begin/end: {}/{}",
                e,
                std_cigar,
                biwfa_cigar,
                std_aligner.cigar().begin_offset,
                std_aligner.cigar().end_offset,
                biwfa_aligner.cigar().begin_offset,
                biwfa_aligner.cigar().end_offset,
            );
        }
    }

    #[test]
    fn test_biwfa_identical_sequences() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let score = aligner.align_biwfa(b"ACGTACGT", b"ACGTACGT");
        assert_eq!(score, 0);
        aligner
            .cigar()
            .check_alignment(b"ACGTACGT", b"ACGTACGT")
            .unwrap();
    }

    #[test]
    fn test_biwfa_edit_simple() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let pattern = b"ACGTACGT";
        let text = b"ACTTACGT";
        let biwfa_score = aligner.align_biwfa(pattern, text);
        aligner.cigar().check_alignment(pattern, text).unwrap();

        // Compare with standard WFA
        let mut std_aligner = EditAligner::new(WavefrontPenalties::new_edit());
        std_aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let std_score = std_aligner.align_end2end(pattern, text);
        assert_eq!(biwfa_score, std_score);
    }

    #[test]
    fn test_biwfa_edit_with_indels() {
        let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());
        let pattern = b"ACGTACGTACGT";
        let text = b"ACGACGTAGT";
        let biwfa_score = aligner.align_biwfa(pattern, text);
        aligner.cigar().check_alignment(pattern, text).unwrap();

        let mut std_aligner = EditAligner::new(WavefrontPenalties::new_edit());
        std_aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let std_score = std_aligner.align_end2end(pattern, text);
        assert_eq!(biwfa_score, std_score);
    }

    #[test]
    fn test_biwfa_affine_simple() {
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        let pattern = b"ACGTACGTACGTACGT";
        let text = b"ACTTACGTACGTACTT";
        let biwfa_score = aligner.align_biwfa(pattern, text);
        aligner.cigar().check_alignment(pattern, text).unwrap();

        let mut std_aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        std_aligner.alignment_scope = AlignmentScope::ComputeAlignment;
        let std_score = std_aligner.align_end2end(pattern, text);
        assert_eq!(biwfa_score, std_score);
    }

    /// Deterministic 300 bp pair with a few edits (long enough for BiWFA to
    /// search for a breakpoint rather than fall back to unidirectional WFA).
    fn biwfa_test_pair() -> (Vec<u8>, Vec<u8>) {
        let mut x = 12345u32;
        let pattern: Vec<u8> = (0..300)
            .map(|_| {
                x = x.wrapping_mul(1103515245).wrapping_add(12345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect();
        let mut text = pattern.clone();
        text[50] = if text[50] == b'A' { b'C' } else { b'A' };
        text.remove(120);
        text.insert(200, b'G');
        (pattern, text)
    }

    #[test]
    fn test_piggyback_compaction_keeps_alignment() {
        // Unrelated 3 kb sequences: affine score ~6900, past the first
        // compaction probe (score 6000).
        let mut x = 7u32;
        let mut seq = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(1103515245).wrapping_add(12345);
                    b"ACGT"[(x >> 16) as usize % 4]
                })
                .collect()
        };
        let (pattern, text) = (seq(3000), seq(3000));
        let align = |mode: MemoryMode, compact_at: u64| {
            let mut a = AffineAligner::new(WavefrontPenalties::new_affine(
                crate::penalties::AffinePenalties {
                    match_: 0,
                    mismatch: 4,
                    gap_opening: 6,
                    gap_extension: 2,
                },
            ));
            a.set_alignment_scope(AlignmentScope::ComputeAlignment);
            a.set_memory_mode(mode);
            a.set_max_memory_compact(compact_at);
            let score = a.align_end2end(&pattern, &text);
            assert!(
                score > 6000,
                "score {score} must reach the compaction probe"
            );
            let compactions = a.wf_components.bt_buffer.as_ref().unwrap().num_compactions;
            (a.cigar().to_string_rle(true), compactions)
        };
        for mode in [MemoryMode::Med, MemoryMode::Low] {
            let (plain, n) = align(mode, u64::MAX);
            assert_eq!(n, 0);
            let (compacted, n) = align(mode, 0);
            assert_eq!(compacted, plain, "{mode:?}");
            // Med keeps counting compactions; low resets after each one.
            assert_eq!(n, if mode == MemoryMode::Med { 1 } else { 0 });
        }
    }

    #[test]
    fn test_biwfa_score_matches_end2end() {
        let (pattern, text) = biwfa_test_pair();
        let penalties = WavefrontPenalties::new_affine(crate::penalties::AffinePenalties {
            match_: 0,
            mismatch: 4,
            gap_opening: 6,
            gap_extension: 2,
        });
        let mut aligner = AffineAligner::new(penalties.clone());
        let score = aligner.align_biwfa_score(&pattern, &text);
        assert_eq!(aligner.status(), AlignStatus::Completed);
        let mut std_aligner = AffineAligner::new(penalties);
        assert_eq!(score, std_aligner.align_end2end(&pattern, &text));
    }

    #[test]
    fn test_biwfa_band_excluding_end_is_unreachable() {
        // End diagonal k = -30 is outside the [-5, 5] band, which the
        // sub-aligners now inherit (as C's Marco-Sola-Lab/WFA2-lib#125); previously
        // BiWFA ignored the heuristic and completed.
        let (pattern, _) = biwfa_test_pair();
        let text = &pattern[..270];
        let mut aligner = AffineAligner::new(WavefrontPenalties::new_affine(
            crate::penalties::AffinePenalties {
                match_: 0,
                mismatch: 4,
                gap_opening: 6,
                gap_extension: 2,
            },
        ));
        aligner.set_heuristic(HeuristicStrategy::BandedStatic {
            min_k: -5,
            max_k: 5,
        });
        aligner.align_biwfa(&pattern, text);
        assert_eq!(aligner.status(), AlignStatus::EndUnreachable);
        assert!(aligner.cigar().operations_slice().is_empty());
        aligner.align_biwfa_score(&pattern, text);
        assert_eq!(aligner.status(), AlignStatus::EndUnreachable);
        // A band containing the end diagonal still aligns.
        aligner.set_heuristic(HeuristicStrategy::BandedStatic {
            min_k: -40,
            max_k: 5,
        });
        aligner.align_biwfa(&pattern, text);
        assert_eq!(aligner.status(), AlignStatus::Completed);
        aligner.cigar().check_alignment(&pattern, text).unwrap();
    }
}
