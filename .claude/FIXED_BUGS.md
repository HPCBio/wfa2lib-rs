# Fixed Bugs

Non-trivial bugs found and fixed in the Rust port, as asked for in
`PROJECT.md`. Each entry records the symptom, the root cause, the C reference
behaviour, which commit introduced the bug, and the fix.

C reference: WFA2-lib `16969b5b` (2026-02-28), the commit the port was made
from (see "Port baseline" in `PROGRESS.md`), plus current C `main` where noted.
Issue and PR numbers refer to the HPCBio/wfa2lib-rs fork unless marked
COMBINE-lab.

| # | Bug | Introduced | Status |
|---|---|---|---|
| 1 | Out-of-bounds write past score-indexed component arrays | d235d45 | Fixed: `d3618de` (HPCBio/wfa2lib-rs#11, merged); upstream COMBINE-lab/wfa2lib-rs#6 (open) |
| 2 | `init_ends_{lower,higher}` negative length wraps to a huge slice | 4d99296 | Fixed: `e494181` (HPCBio/wfa2lib-rs#11, merged); upstream COMBINE-lab/wfa2lib-rs#6 (open) |
| 3 | Edit/indel/linear kernels overrun the wavefront allocation | original port | Fixed: `06c593c` (HPCBio/wfa2lib-rs#12, merged `12ba05b`) |
| 4 | Unreachable banded alignments never terminate; affine banded output differs from C | 2c112b2, 35bfcde + original port | Fixed: `320611d`, `ddae4d8`, `29a9c89` (HPCBio/wfa2lib-rs#14, merged `59db5a4`) |
| 5 | CLI prints `0` and an empty CIGAR field when no alignment is found (C: `-1  -`) | original port (CLI) | Fixed: `5e6ce55` (HPCBio/wfa2lib-rs#15, merged) |
| 6 | Ends-free begin gaps ignored (output equals global) | original port | Fixed: `28cb857`, `6a06d8f` (HPCBio/wfa2lib-rs#17, merged) |
| 7 | X-drop/Z-drop had no effect; backtrace and null-step counting differed from C | original port | Fix on branch `fix/xdrop-zdrop` |

All four were found by differential testing against C (HPCBio/wfa2lib-rs#10): the
Rust and C `align_benchmark` binaries run over the same configuration matrix,
which includes banded heuristics that the stock `wfa.utest.sh` doesn't cover.

---

## 1. Out-of-bounds write past the score-indexed component arrays

- **Symptom:** segfault (release) or UB-check abort (debug) with
  `banded-static` whenever the band excludes the end diagonal, in full-memory
  (CIGAR) mode. Reproduce: `>` 120 bp / `<` its first 90 bp, `-a edit-wfa
  --wfa-heuristic=banded-static --wfa-heuristic-parameters=-5,5`.
- **Root cause:** an unreachable alignment keeps raising the score past
  `WavefrontComponents::num_wavefronts`. `set_*_ptr` writes through
  `get_unchecked_mut`. Running past `num_wavefronts` also exceeds the slab
  reservation (`num_wavefronts * num_types`) that keeps the stored
  `*mut Wavefront` pointers valid.
- **Introduced:** d235d45 removed the original grow-on-demand bounds checks, on
  the assumption (stated in a `// SAFETY:` comment) that `compute_dimensions()`
  pre-sizes the arrays to the maximum reachable score. That only holds when the
  alignment terminates.
- **C:** `wavefront_compute_allocate_output` checks
  `score_mod >= num_wavefronts` ("Maximum allocated wavefronts reached") and
  calls `exit(1)`.
- **Fix (`d3618de`):** `WavefrontComponents::has_slot()` is checked once per
  score step in the end-to-end, ends-free and BiWFA-base loops, stopping with
  `AlignStatus::MaxStepsReached`. The hot accessors stay unchecked.
- **Test:** `test_banded_unreachable_does_not_overflow_components`.

## 2. `Wavefront::init_ends_{lower,higher}` wrapped slice length

- **Symptom:** UB-check abort at `wavefront.rs` (`slice::from_raw_parts_mut`)
  with gap-linear + banded-static; silent memory corruption in release.
- **Root cause:** when a cutoff leaves `lo < min_lo < wf_elements_init_min`,
  `(min_init - start) as usize` is negative and wraps.
  `init_ends_higher` has the mirror-image case.
- **Introduced:** the original `for k in start..min_init` loop matched C.
  cbe976e changed it to a checked slice range (which would panic instead), and
  4d99296 to `from_raw_parts_mut` with an unchecked cast.
- **C:** `wavefront_compute_init_ends_wf_{lower,higher}` use `for` loops, which
  do nothing when the range is empty.
- **Fix (`e494181`):** fill only when the range is non-empty.
- **Tests:** `test_init_ends_{lower,higher}_empty_range`,
  `test_banded_linear_init_ends_empty_range`.

## 3. Edit/indel/linear kernels overrun the wavefront allocation

- **Symptom:** heap corruption with banded-static ±10.
  - gap-linear: on `wfa.utest.seq` pair 251 alone (`lo = -1265`).
  - edit: depends on order, appearing on pair 261 only after pair 251.
  - score-only mode never terminated.
- **Root cause:** once trimming empties M (`null = true`, `lo > hi`), the edit
  and linear steps kept computing from it, because they only checked for a
  missing *pointer*. Each step widened the range by one diagonal. The band
  cutoff skips `lo > hi` wavefronts, so `lo/hi` drifted past the historic range
  the wavefronts are allocated with.
- **C:** linear/affine/affine2p check `->null` on every input. If all are null
  they store a NULL output and count a null step, and extend reports
  `END_UNREACHABLE` after `max_score_scope` such steps. (C's edit step lacks
  the check and never terminates: a C bug, still to be reported upstream.)
- **Fix (`06c593c`, HPCBio/wfa2lib-rs#12, merged):**
  - null output when every input is null;
  - edit/linear `lo/hi` clamped to the historic range minus the kernel's k±1
    margin, as affine already does;
  - null outputs release the modular slot.
- **Test:** `test_banded_unreachable_edit_linear_end_unreachable`.

## 4. Banded alignments: no termination and divergent results (affine/affine2p)

- **Symptoms:**
  - affine/affine2p score-only mode hung on unreachable banded alignments;
  - CIGAR mode gave up on pairs that C completes (utest pairs 233 and 293);
  - give-up scores differed from C.
- **Root causes** (three; each one alone left differences from C):
  1. `apply_heuristic_cutoff` equated I/D to M on every call. C's
     `wavefront_heuristic_cufoff` skips an empty or unchanged M, and marks M
     null if the cutoff emptied it. Equating unconditionally empties I/D
     whenever M is empty.
  2. 2c112b2 dropped `trim_ends` on I/D wavefronts ("matching C's
     process_ends"), and 35bfcde replaced it with clamping I/D to M's range.
     C's `wavefront_compute_process_ends` trims every output independently, so
     I/D keep valid diagonals outside M and can be recognised as null.
  3. Inputs that were present but trimmed empty were read as live, and the
     range union skipped missing inputs. C's `wavefront_compute_get_*wavefront`
     substitutes `wavefront_null` for missing *or* `->null` inputs, and
     `wavefront_compute_limits_input` includes its bounds (lo = 1, hi = -1) in
     the union.
- **How it was found:** per-score traces of M/I/D `lo/hi/null` and the null-step
  count, from Rust and from an instrumented copy of C, on utest pair 202. The
  traces first diverge at s = 980 (M = [-11,-5] in Rust vs [-11,0] in C).
- **Fix (HPCBio/wfa2lib-rs#14, merged):** `320611d` (cutoff equate),
  `ddae4d8` (I/D trim), `29a9c89` (null inputs).
- **Result:** banded-static output for affine, affine2p and gap-linear matches
  C `16969b5b` exactly on `wfa.utest.seq` and a 300 × 2 kb / 8% set (bands
  ±10 and ±50, alignment and score-only modes), including the give-up score for
  unreachable pairs.
- **Tests:** `test_banded_unreachable_affine_end_unreachable`,
  `test_banded_affine_equates_id_only_when_m_narrowed`.

## 5. CLI output when no alignment is found

- **Symptom:** in alignment (CIGAR) mode, `align_benchmark` printed
  `0<TAB>` (score 0, empty CIGAR field) for pairs with no alignment, e.g.
  END_UNREACHABLE under a band. C prints `-1<TAB>-`. Because of this, banded
  output never diffed clean against C, even where the alignments matched.
- **Root cause:** the CLI always re-scored the CIGAR (an empty CIGAR scores 0)
  and printed its RLE string (empty).
- **C:** `benchmark_print_output` uses score `-1` when the CIGAR is empty
  (`begin_offset >= end_offset`) and prints `-` for an empty CIGAR.
- **Fix (`5e6ce55`, HPCBio/wfa2lib-rs#15, merged):** the same rule in
  `compute_output_score` and the CIGAR column. Banded-static output is now
  byte-identical to C `16969b5b` in high, med and score-only modes (30/30
  cells).
- **Tests:** `tests/cli_output.rs` runs the binary: an empty CIGAR prints
  `-1  -`, and the completed and score-only formats are unchanged.

## 6. Ends-free begin gaps ignored

- **Symptom:** with `set_alignment_free_ends(P0, Pf, T0, Tf)` and P0 or T0 > 0,
  the result was identical to the global alignment. Visible once the CLI
  honoured `--wfa-span` (HPCBio/wfa2lib-rs#16): every begin-gap cell differed from C,
  for every metric.
- **Root cause:** `pattern_begin_free` / `text_begin_free` were stored but
  never read. C has two mechanisms, and neither was ported:
  - `match == 0`: `wavefront_aligner_init_wf_m` seeds M0 over `[-P0, T0]`,
    putting cell (h,0) on k = h and (0,v) on k = -v.
  - `match != 0`: `wavefront_compute_endsfree_{required,limits,init,allocate_null}`
    add cell (d,0) / (0,d) at score d * -match, before the ends are trimmed.
- **Fix (HPCBio/wfa2lib-rs#17, merged):** `28cb857` (match == 0) and
  `6a06d8f` (match != 0). Also C's extension preset (P0 = T0 = 0) and smarco/WFA2-lib#126's
  restriction of the seeding to the ends-free span.
- **Parity details reproduced:**
  - C sets only `lo/hi` on M0; the initialized range stays `[0, 0]`.
  - C's edit step nulls the boundary cells directly instead of calling
    `init_ends`.
  - Together these mean the edit/indel backtrace can't step onto M0's
    begin-free diagonals. It falls back to `min(v,h)` unchecked `M`s, so
    **C can emit an invalid CIGAR for edit/indel with begin gaps** (its own
    `--check=correct` flags it; e.g. `ACGTTCGTACGA` vs `CCACGTACGTACGA`,
    T0 = 3 → `2I12M`).
  - Scores are correct. Gap metrics are unaffected, because their
    `init_ends` widens the initialized range.
  - It's a C regression from `380eb31b` (backtrace bounds checks moved
    from `lo/hi` to the initialized range). C at `380eb31b^` gives valid
    CIGARs.
  - Rust reproduces it for output identity, documented on
    `set_alignment_free_ends`. **Known issue HPCBio/wfa2lib-rs#18.**
  - Reported upstream as smarco/WFA2-lib#142, with a reproducer and a
    verified one-line fix: set `wf_elements_init_min/max = lo/hi` on M0 in
    `wavefront_aligner_init_wf_m`.
- **Result:** ends-free (begin and/or end gaps) and extension are
  byte-identical to C for all metrics, match == 0 and match < 0, in
  alignment, score-only and med modes.
- **Tests:** `test_endsfree_begin_free_match_zero`,
  `test_extension_ignores_begin_free`, and
  `tests/cli_output.rs::ends_free_begin_gaps_match_c`.

## 7. X-drop / Z-drop had no effect; backtrace and null-step counting differed

- **Symptom:** with xdrop/zdrop, end-to-end output equalled the unbanded
  alignment. After porting the heuristics, three more differences from C
  surfaced under pruning.
- **Root causes:**
  1. *Heuristics.* The SW-score formula differed (C:
     `(swg_match*(v+h) - s)/2`, `swg_match = -match` or -1 when match == 0).
     There was no `steps_between_cutoffs`. The first call pruned instead of
     only recording. X-drop updated the best before pruning. Z-drop never set
     the end cell, and Z-dropped alignments weren't backtraced as partial
     alignments.
  2. *Backtrace.* Rust's affine state M compared the mismatch with the
     stored I1/D1 cells at score s; C reads every predecessor from its source
     wavefront. These are equivalent until pruning narrows the computed
     range.
  3. *Null steps.* Rust counted them in extend; C counts in compute, including
     ends-free begin-gap "null" M-wavefronts (so Rust gave up much later with
     ends-free + match < 0).
- **C details reproduced:**
  - Z-drop backtraces from the best cell starting at the *current* score.
  - The CIGAR buffer is `2*(plen+tlen)`.
  - The end position starts unset (`k = INT_MAX`, `offset = OFFSET_NULL`).
  - Score-only output for unreachable pairs is evaluated there, with wrapping
    (HPCBio/wfa2lib-rs#23).
  - A negative match run moves the CIGAR write position forward
    (HPCBio/wfa2lib-rs#21).
- **CLI:** xdrop/zdrop require P1,P2 and are rejected with edit/indel, as in C.
- **Fix (branch `fix/xdrop-zdrop`):** `a3a4c4f` (backtrace port), `763ed8a`
  (null-step counting), `f4eb349` (heuristics, partial alignments, CLI).
- **Result:** 450/576 X-drop/Z-drop cells byte-identical to C 3282f83 (up
  from 0); all 126 others are `med` (C piggyback; HPCBio/wfa2lib-rs#8,
  HPCBio/wfa2lib-rs#22). Everything else is unchanged.
- **Tests:** `tests/cli_output.rs::{xdrop_zdrop_match_c,
  drops_reject_what_c_rejects}`, `heuristic::tests::{test_compute_sw_scores_c_formula,
  test_drop_state_reset_uses_steps}`.


## 8. BiWFA ignored heuristics; banded-adaptive and banded limits differed from C

- **Symptom:** with `--wfa-memory=ultralow`, every heuristic was ignored: the
  BiWFA sub-aligners and the base-case aligner ran unpruned. Separately,
  banded-adaptive applied the static band every N steps, unlike C, and
  unidirectional banded output followed the port baseline rather than C HEAD.
- **Root causes:**
  1. *Banded-adaptive.* C samples the wavefront at 4 points and shifts the
     `[lo, hi]` window towards the furthest-reaching diagonals; Rust clamped
     to the static band.
  2. *smarco/WFA2-lib#125 (C 6565731 / 7ea485e).* C HEAD clamps computed
     limits to `[min_k, max_k]` in `wavefront_compute_limits_heuristic`
     (linear, affine, affine2p), and shifts the band into each BiWFA
     sub-problem. C also creates components for 1000 × 1000 up front, which
     decides whether a banded alignment runs out of score slots.
  3. *BiWFA status handling.* C's forward and reverse aligners inherit the
     heuristic (band shifted per sub-problem). Their extend step includes the
     null-step check, termination and the cutoff. When one direction reaches
     the end first, C falls back to unidirectional WFA if the end score is at
     most `WF_BIALIGN_RECOVERY_MIN_SCORE` (500), and otherwise reports
     unreachable. The base case takes the heuristic too (band shifted). It
     appends its CIGAR whether or not it completed (empty, or max-trimmed
     after a Z-drop), and a failure makes the whole alignment UNATTAINABLE.
  4. *BiWFA score-only.* C's `wavefront_bialign_compute_score` only searches
     for the breakpoint. On failure the score stays at `cigar_clear`'s
     `INT32_MIN`. Rust ran unidirectional WFA, which gives the same score
     without heuristics but not with them.
- **Fix (branch `feat/biwfa-heuristics-banded`):** `19eac8b` (banded-adaptive),
  `3b16283` (#125 clamp, C-sized components), then the BiWFA commit (status
  plumbing via `BiStatus`, heuristic inheritance, base case as C's
  `wavefront_unialign`, `align_biwfa_score`). BiWFA failures report
  `AlignStatus::EndUnreachable`; C's UNATTAINABLE has no `AlignStatus`
  equivalent, and adding one would change the public enum.
- **CLI note:** C's driver uses −1 as "unset" for heuristic parameters, so it
  rejects a band edge of −1 (e.g. `-1,1`). Rust accepts it.
- **Result:** BiWFA × {none, banded-static, banded-adaptive, xdrop, zdrop,
  wfa-adaptive} × {edit, linear, affine, affine match < 0, affine2p} ×
  {CIGAR, score-only} × {utest, long}: 168/168 byte-identical to C 3282f83.
  Before the fix, every heuristic run diverged.
- **Tests:** `tests/cli_output.rs::biwfa_heuristics_match_c`,
  `aligner::tests::{test_biwfa_band_excluding_end_is_unreachable,
  test_biwfa_score_matches_end2end}`.

## 9. BiWFA tie-breaks with match < 0 (HPCBio/wfa2lib-rs#7)

- **Symptom:** 3 stock `wfa.utest.sh` BiWFA cases with match < 0 gave a
  different co-optimal CIGAR than C (same score).
- **Root cause:** Rust fell back to unidirectional WFA for any sub-problem of
  at most 100 bp, at every recursion level. In C,
  `WF_BIALIGN_FALLBACK_MIN_LENGTH` applies only at the top level; deeper
  levels fall back only by score (`score_remaining <= 250`). Splitting at
  different points can pick different co-optimal paths.
- **Fix:** check the length only at the top level, as C does.
- **Result:** stock `wfa.utest.sh` 60/60 identical (was 57/60, the other 3
  matching by score only). BiWFA match < 0 is now byte-identical to C on utest
  and long, with and without heuristics.
- **Test:** `tests/cli_output.rs::biwfa_match_reward_tie_breaks_match_c`
  (`wfa.utest.seq` pair 298).

## 10. `--wfa-memory=med|low` behaved like `high` (no piggyback backtrace)

- **Symptom:** `med` and `low` gave `high`'s results and memory use. C's
  `med`/`low` record the alignment alongside the wavefronts ("piggyback"
  backtrace) and keep only the wavefronts still needed. That can pick a
  different co-optimal path, and C's heuristic results differ between the
  modes (all X-drop/Z-drop `med` cells of the matrix differed).
- **C design ported** (`wavefront_compute_*_piggyback`,
  `wavefront_backtrace_offload.c`, `wavefront_backtrace_buffer.c`,
  `wavefront_pcigar.c`):
  - Each cell carries a 16-operation packed CIGAR (`pcigar`) and the index
    of the backtrace block it continues (`bt_prev`). The kernels copy them
    from the winning predecessor with C's tie order and push the operation.
    Matches aren't stored; unpacking re-extends them.
  - When a component's occupancy nears 16 operations, live cells' pcigars
    move to the backtrace buffer (`med`: at 9+ operations; `low`: when full,
    or 15+ for affine). Edit/indel offload every 16 scores.
  - Initial cells (score 0, ends-free begin cells) start chains at init
    blocks holding their (v, h).
  - Unpacking reproduces C exactly, including: begin I/D before the first
    block; padding to the sequence ends with matches then I/D; the dropped X
    after I/D (gap close); and a negative match count moving the position
    back.
  - After a Z-drop, C reads the end cell from the *current* score's
    wavefront, not the drop score's (see below).
- **Rust specifics:**
  - Kernels are split into per-component passes over slices (noalias), so
    they vectorize as C's `PRAGMA_LOOP_VECTORIZE` loops do. A fused
    raw-pointer loop ran at 40+ ns per cell from false store/load aliasing.
  - Backtrace-buffer segments commit memory as blocks are stored (C:
    `calloc`). The first version zero-filled 64 MB per aligner.
- **#22 root cause (C bug, not reproducible):** after a Z-drop, C backtraces
  from the end diagonal in the current score's M-wavefront. That cell is
  often dead (NULL offset) or outside the computed range, so its chain comes
  from memory C never initialized: leftovers of earlier steps or earlier
  alignments, or zeroes in fresh memory. C's result for a pair then depends
  on what it aligned before. Rust reads never-written cells, and the NULL
  padding the kernels read, as zeroes. Each result is deterministic and
  matches C with fresh memory (pair aligned alone) in all but a few cases.
- **Result:**
  - `med`/`low` × {edit, indel, linear, affine, affine match < 0, affine
    `0,1,2,1`, affine2p} × {utest, long}: 28/28 byte-identical to C 3282f83.
  - X-drop/Z-drop matrix: 556/576 (was 450). The other 20 are Z-drop `med`
    with affine match < 0 (C's uninitialized reads above).
  - 50k × 1 kb: byte-identical to C, and faster than C (affine `med` 8.3 s
    vs 11.0 s). Max RSS 3.6 MB (`high`: 11.3 MB).
- **Tests:** `tests/cli_output.rs::memory_med_low_match_c` (fixtures in
  `tests/data/`).
