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
  `6a06d8f` (match != 0). Also C's extension preset (P0 = T0 = 0) and Marco-Sola-Lab/WFA2-lib#126's
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
  - Rust reproduced it for output identity (HPCBio/wfa2lib-rs#18).
  - Reported upstream as Marco-Sola-Lab/WFA2-lib#142, with a reproducer and a
    verified one-line fix: set `wf_elements_init_min/max = lo/hi` on M0 in
    `wavefront_aligner_init_wf_m`. Merged as Marco-Sola-Lab/WFA2-lib#144
    (6176f94) and ported (see 14b).
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
  2. *Marco-Sola-Lab/WFA2-lib#125 (C 6565731 / 7ea485e).* C HEAD clamps computed
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

## 11. Backtrace buffer never compacted (med/low)

- **Symptom:** with piggyback backtrace, the backtrace buffer only grows. C
  compacts it during long alignments, so `med`/`low` stay bounded.
- **C design ported** (`wavefront_unialign_reached_limits`,
  `wavefront_components_compact_bt_buffer`,
  `wf_backtrace_buffer_compact_marked`):
  - Every 6000 scores, if the buffer exceeds `max_memory_compact` (256 MB),
    C marks the blocks reachable from live cells of the wavefronts still in
    scope. It moves those blocks down after the already compacted ones,
    translating their `prev` indices by exclusive rank, and translates the
    cells' indices.
  - After `max_partial_compacts` compactions (4 for `med`, 1 for `low`),
    the next one starts over from block 0.
  - C quirks kept: the compacted-block count is set one short (the last
    index); only live cells are translated; and the "telescopic" limit
    update assigns a local, so the threshold stays 256 MB.
- **Rust:** `set_max_memory_compact` exposes the threshold (fixed in C) so
  compaction can be tested on small inputs. It never changes results.
- **Test:** `aligner::tests::test_piggyback_compaction_keeps_alignment`
  (forced compaction at score 6000 gives the same CIGAR, `med` and `low`).

## 12. Med/low: begin-free seed at score 0 backtraced from (0, 0) (match == 0)

- **Symptom:** ends-free with begin gaps, match == 0, `med`/`low`: if the
  alignment ended on a begin-free seed's diagonal at score 0, the CIGAR
  started at (0, 0) instead of at the seed, and was invalid. Example: gap-linear `0,3,3`,
  `P=G`, `T=TCGGGACATCGAG`, ends-free `0,0,4,11`: Rust `med` gave `1M12I`,
  C and Rust `high` give `2I1M10I`.
- **Cause:** `init_alignment` set the seeds' backtrace blocks but not the
  wavefront's written range, so `bt_or_zero` read them as unwritten (block
  0, the (0, 0) start). C reads the cells directly.
- **Fix:** mark the seeded range `[-P0, T0]` written. (Unrelated to
  HPCBio/wfa2lib-rs#18, which is about the *initialized* range in C's
  high-memory backtrace.)
- **Test:** `aligner::tests::test_endsfree_begin_free_seed_backtrace_med_low`.
  Found by a DP-oracle sweep (HPCBio/wfa2lib-rs#32 work).

## 13. Affine2p modular slots leaked I2/D2 wavefronts (memory safety)

- **Symptom:** score-only (or med/low) affine2p could segfault. Debug builds
  hit `wavefront slab Vec would reallocate`. Found with ends-free and
  match < 0, but the bug is in the step itself.
- **Cause:** when no I2/D2 input is live, `compute_affine2p_step` uses the
  affine kernel and stores NULL for I2/D2. In modular mode it overwrote the
  slot's previous I2/D2 wavefronts without returning them to the slab. The
  slab then outgrew the capacity reserved in `init_alignment`, its Vec
  reallocated, and every wavefront pointer held by the components dangled.
- **Fix:** free the slot's old I2/D2 before storing NULL (as
  `reuse_or_allocate` does for stored wavefronts).
- **Test:** `aligner::tests::test_affine2p_modular_slots_release_i2_d2`.

## 14. Corrected mode (`Compat`, HPCBio/wfa2lib-rs#32)

`Compat::Corrected` (CLI `--no-c-compat`) is a per-aligner opt-in that fixes
known C bugs; the default `Compat::C` stays output-identical to C. It covers:

### 14a. Ends-free with match < 0 (HPCBio/wfa2lib-rs#37, HPCBio/wfa2lib-rs#38, Marco-Sola-Lab/WFA2-lib#102)

- **Symptom (C and C-compatible Rust):** with a match reward, ends-free
  alignment can return a suboptimal result: a trailing-end miss of a point
  or two (Marco-Sola-Lab/WFA2-lib#102), a shift by one repeat unit (#37), or, with
  free ends spanning whole sequences, the all-gap alignment (#38).
- **Root cause:** with `match < 0`, WFA scores an alignment as
  `2 * penalty + M * (aligned bases)` (M = -match). For global alignment
  `aligned bases = |P| + |T|` is fixed, so this ranks alignments correctly.
  For ends-free it isn't fixed. C already charges begin-free gaps M per base
  (the begin-free cell (0, d) or (d, 0) is seeded at score `d * M`). But it
  ends at the first cell to reach a free end, charging nothing for the
  trailing bases left over. An alignment that skips more trailing bases
  scores too low, and the first one found wins.
- **Fix (opt-in, `Compat::Corrected`, `--no-c-compat`):** a cell reaching
  a free end at score `s` with `left` bases to skip is a complete alignment
  of total `s + M * left`, which is `2 * penalty + M * (|P| + |T|)` for
  every alignment. The extend step extends all diagonals, keeps the best
  candidate (fewest left, then lowest k), and stops once the score reaches
  its total. No later cell can do better. The backtrace starts from the
  candidate's own score. Med/low keep its packed-CIGAR chain through buffer
  compaction (it is marked and translated as a pinned chain). The returned
  score is the total, so `wf_score_to_sw_score(|P|, |T|, score)` is the
  real score. With match == 0 the output is unchanged (the totals equal the
  scores).
- **Cost:** it runs on until the optimal total, not just the first boundary
  hit. Pruning keeps that cheap: from cell (v, h) any completion totals at
  least `s + M * |k - k_end|` (`k_end = |T| - |P|`; the bases left differ by
  `k - k_end` all along a diagonal, and each of them is a gap or free-end
  base, at least M). Once there is a candidate, the cells at each score are
  trimmed to the diagonals that could still beat it (ties kept) and the
  dropped cells nulled, so no later step reads them. The result is
  unchanged. Only with no heuristic or `BandedStatic`, since the other
  heuristics read the wavefronts. The extend kernel checks C's cheap
  termination test first and counts the bases left only at an end. With
  both, corrected mode is +2–2.5% instructions over C mode when the first
  end found is near-optimal; when it isn't, the extra search is what fixes
  the result. Measured in `PERFORMANCE_LOG.md`.
- **Tests:** `tests/corrected_endsfree.rs`: the #37/#38 reproducers, and
  3000 random cases (linear/affine/affine2p, small alphabets and repeats,
  whole-length, extension-like and random free ends). Each is checked
  against a Gotoh DP oracle for score-only and for the CIGAR in
  high/med/low. Also: match == 0 identical to C mode; a pending candidate
  survives forced compaction in med/low; #38's pair B (16S, 1.4 kb).
  240k further random cases (up to 40 bp) passed during development.
  `random_long_against_oracle`: 300 longer pairs (to ~200 bp, with flanks)
  where pruning trims most cells; a band one diagonal narrower fails it.

### 14c. Corrected mode: null steps ended the search with a candidate pending

- **Symptom:** corrected mode with `BandedStatic` could return a worse end
  than the optimum within the band (e.g. SW 1 instead of 2) when the end
  diagonal lies outside the band.
- **Root cause:** the band (and, with pruning, the pruning band) can empty
  every wavefront well before the candidate's total. Ending on null steps,
  as C does, stopped the search there, but begin-free seeds still to come
  (one per `M` score, outside the band) can beat the candidate.
- **Fix:** with a candidate, null steps don't end the search; it ends when
  the score reaches the candidate's total, as before. Without a candidate,
  unchanged. Found while checking the pruning (14a) for identical output.
- **Test:** `banded_static_pending_seeds` (`tests/data/corrected_banded_seeds.seq`,
  linear `-2,4,3`, whole-length free ends, band `-10,200`) vs the DP oracle,
  high/med/low.

### 14b. Edit/indel ends-free with begin gaps: invalid CIGAR (HPCBio/wfa2lib-rs#18)

- **C bug:** a regression from `380eb31b` (reported as Marco-Sola-Lab/WFA2-lib#142).
  For match == 0, M0's begin-free seeds get `lo/hi` but the initialized range
  stays `[0, 0]`, and edit/indel never widen it. So the backtrace can't step
  onto the seeds and emits unchecked `M`s.
- **Fix (corrected mode):** set M0's initialized range to `[-P0, T0]`, the
  fix proposed upstream. Gap metrics are unaffected (identical in both modes).
- **Upstream (2026-10-17):** merged as Marco-Sola-Lab/WFA2-lib#144 (C 6176f94).
  `C_COMMIT` bumped and the fix applied in both modes; no longer a
  corrected-mode difference.
- **Tests:** `issue_18_edit_begin_gaps` (now checks C mode), `edit_indel_random_against_oracle`
  (3000 cases vs the DP oracle, all memory modes).

### 14d. X-drop with steps > 1: invalid CIGARs (HPCBio/wfa2lib-rs#21)

- **C bug:** with `steps_between_cutoffs > 1` only every n-th wavefront is
  pruned, and a cutoff can trim it to empty (`lo > hi`, `null`). Compute
  reads an empty input as the null wavefront, but its offsets stay within
  the initialized range, and the high-memory backtrace bounds its reads by
  that range only. So it can take a stale source that compute never used.
  If that source's offset exceeds the current cell's, the match run is
  negative, and C moves the CIGAR write position forward over earlier
  operations (e.g. `wfa.utest.seq` pair 201, affine `50,10`: at score 506,
  k = 6, offset 547, the mismatch source M(502) was trimmed empty but still
  holds 548). C's `med` piggyback backtrace records the real sources, so it
  is valid.
- **Fix (corrected mode):** when a cutoff trims M (or I/D, by equating) to
  empty, also empty its initialized range (`Wavefront::clear_init_range`),
  so the backtrace reads it as null, as compute does. Scores are unchanged
  in all modes.
- **Result:** X-drop `50,10` / `100,5` × affine / linear / affine2p ×
  {`wfa.utest.seq`, `long.seq`}: every non-empty CIGAR is valid (C mode:
  e.g. affine `100,5` on `long.seq` 0/300 valid). Score-only output
  identical in both modes.
- **Not a bug:** with pruning, the WFA path can open a gap right after
  another (I-open, then I-open with no match between), when the extension
  source was pruned. The CIGAR merges them into one gap, so it rescores
  below the WFA score (never above): e.g. pair 298, affine `100,5`, 242 vs
  248. C's valid CIGARs show the same (affine2p `50,10`: 5/305 pairs).
- **Tests:** `tests/corrected_xdrop.rs` (all three gap metrics). First pair
  201 with default penalties; since 14f changes X-drop's scores at
  `match == 0`, pair 305 (`tests/data/issue21_xdrop.seq`) with match-reward
  penalties and X-drop `20,2`.

### 14e. Z-drop ends backtraced at the wrong score (HPCBio/wfa2lib-rs#22)

- **C bug:** on a Z-drop, `wavefront_unialign_terminate` backtraces the best
  cell's diagonal and offset (found at `max_wf_score`) starting at the
  *current* score, where that cell needn't exist.
  - High memory: the backtrace finds no source, stops, and pads with
    unchecked `M`s plus leading indels. With default penalties that is
    usually max-trimmed to empty; with `match < 0` the CIGARs are invalid
    (e.g. `2I58M` on a 2 kb, 8%-error pair; up to 300/300 invalid on
    `long.seq`, affine `-3,5,3,2`, Z-drop `100,5`).
  - Med/low: the cell at the current score is often dead, and its
    `bt_pcigar/bt_prev` were never initialized (leftovers of earlier
    alignments in C; zeroes in Rust), so C's high and med results differ and
    C's med depends on what it aligned before.
- **Fix (corrected mode):** backtrace the best cell at `max_wf_score`
  (`backtrace_zdrop_end`): high memory from the stored wavefronts; med/low
  from the cell's `(pcigar, prev)` kept when it became the best
  (`zdrop_best_bt`), pinned through buffer compaction (which now pins a
  slice: the ends-free candidate and the Z-drop cell).
- **Result:** Z-drop `20,1` / `50,1` / `50,10` / `100,5` × default and
  `-3,5,3,2` × global / extension / ends-free × {`wfa.utest.seq`,
  `long.seq`}: every partial CIGAR valid, and high = med except co-optimal
  ties that score the same (an empty CIGAR vs `1M`/`2M`, score 0). C mode
  is unchanged.
- **Note:** with `match == 0`, C's Z-drop score makes the start the best
  cell, so this fix alone leaves most such partial alignments empty; 14f
  fixes that.
- **BiWFA:** a base sub-alignment that Z-drops is backtraced the same way
  in C (current score), and also starts in the sub-problem's end component,
  which can be I/D; the best cell is an M cell. Corrected mode starts in M
  at the best cell's score. BiWFA Z-drop `20,1` / `50,10` / `100,5` ×
  default and `-3,5,3,2` × {`wfa.utest.seq`, `long.seq`}: no invalid
  CIGARs (C: up to 64), and with default penalties far fewer empty results
  (`long.seq` `50,10`: 0, C: 300).
- **Tests:** `tests/corrected_zdrop.rs` (`long.seq` pair 1,
  `tests/data/issue22_zdrop.seq`; a generated pair whose Z-drop comes after
  two compaction probes, which fails without the pin; BiWFA on `long.seq`
  pair 21, `tests/data/issue22_zdrop_biwfa.seq`).

### 14f. X-drop/Z-drop with match == 0: sign-flipped SW score (HPCBio/wfa2lib-rs#55)

- **C bug:** a regression from 51fb1aa (2023-03-19). X-drop and Z-drop score
  cells with `(swg_match·(v+h) − wf_score)/2`. Before 51fb1aa, `match == 0`
  used `(v+h) − wf_score` ("use match-score=1"); since then it is
  `swg_match = -1`, so making progress lowers the score.
  - X-drop prunes the diagonals that progress most: with default penalties
    it abandons nearly every divergent pair (`long.seq`: 0–81 of 300
    optimal, typically 0).
  - Z-drop's best cell is the start, so its partial alignment is empty or
    the leading matches (hidden in C by 14e).
- **Fix (corrected mode):** `heuristic::swg_match` returns +1 for
  `match == 0` (C: -1), passed to the heuristics in place of the match
  penalty. The shared formula keeps the `/2`, so +1 has half the scale of
  the pre-51fb1aa special case.
- **Result:** X-drop with default penalties is optimal on 300/300 `long.seq`
  pairs for affine/linear/affine2p at `50,10` and `100,5` (293–300 at
  `50,1`); Z-drop gives real partial alignments (e.g. `50,10` on
  `long.seq`: the full optimal alignment), high = med.
- **Tests:** `test_compute_sw_scores_c_formula` (both stand-ins).

### 14g. Corrected mode didn't reach BiWFA's sub-aligners (Rust bug)

- **Symptom:** under BiWFA, corrected mode gave C's output: the #21 fix
  (14d) and the #55 fix (14f) never applied. #54's "BiWFA identical to C"
  in corrected mode was this, not parity. E.g. BiWFA X-drop `50,10`,
  default penalties, `long.seq`: 107 invalid CIGARs and 2/300 optimal
  scores in both modes.
- **Root cause:** `take_biwfa_subs` hands out the forward, reverse and base
  sub-aligners, which run the heuristics and the base backtrace, but they
  were created with `Compat::C` and never given the parent's mode (only its
  heuristic).
- **Fix:** `take_biwfa_subs` sets their `compat` to the parent's on every
  call (`set_compat` can change it between alignments).
- **Result:** BiWFA X-drop `50,10` / `100,5` / `20,2`, default and
  `-3,5,3,2`, `wfa.utest.seq` and `long.seq`: no invalid CIGARs in
  corrected mode (C mode: up to 213/300); with default penalties 168–303
  optimal (C mode: 0–240).
- **Tests:** `tests/corrected_xdrop.rs` checks BiWFA too.

### 14h. Score-only output without an alignment (HPCBio/wfa2lib-rs#23)

- **C bug:** in score-only mode C always prints a score, even when there is
  no alignment (in CIGAR mode it prints `-1  -`):
  - unreachable (banded, X-drop) with a match reward: the classic score is
    evaluated at the unset end position (`k = INT_MAX`,
    `offset = OFFSET_NULL`); `v = offset - k` is a signed overflow and wraps
    to `v + h == 1`, giving `(swg_match - score) / 2`. E.g. `wfa.utest.seq`
    pair 262, affine `-3,5,3,2`, X-drop `20,1`: `-113`;
  - unreachable with `match == 0`: the negated score where it gave up;
  - a failed BiWFA score search: `cigar_clear`'s `INT32_MIN`.
  (The example in #23, pair 101, is a reachable pair: off by one.)
- **Rust:** the library returns the WFA score with an `EndUnreachable` /
  `Partial` status (its score is documented as meaningless then); only the
  CLI converts it, reproducing C. So the fix is in `align_benchmark`.
- **Fix (corrected mode, CLI):** print `-1` in score-only mode when the
  alignment ended unreachable with no Z-drop cell, or BiWFA's score search
  failed. Z-drop results keep their score (a real cell).
- **Result:** banded-static and X-drop no-alignment pairs print `-1`,
  unidirectional and BiWFA; only score-only output changes. On the sampled
  drop/banded configurations, every `-1` matches a CIGAR-mode `-1  -`,
  except as noted below.
- **Left as is:**
  - BiWFA score-only has no partial result: when its search Z-drops it
    prints `-1` (C: `INT32_MIN`), though CIGAR mode returns the partial
    alignment (e.g. utest affine `-3,5,3,2`, Z-drop `100,5`: 7/305 pairs).
  - BiWFA's score search can succeed where its CIGAR search fails in a
    sub-problem (independent cutoffs).
  - A Z-drop whose CIGAR max-trims to empty keeps its score.
  - With a match reward, `-1` can also be a real score (e.g. utest pair
    281, 10 × 10 bp, extension), as in C's CIGAR mode.
- **Tests:** `unreachable_score_only_corrected_prints_minus_one`
  (`tests/cli_output.rs`, `tests/data/issue23_unreachable.seq`).

## 15. Empty pattern and text: CIGAR offset underflow (debug panic, HPCBio/wfa2lib-rs#41)

- **Symptom:** a fresh aligner given an empty pattern and an empty text,
  with a CIGAR in high memory, panicked in debug builds ("attempt to
  subtract with overflow", `backtrace.rs`). Release builds wrapped and
  happened to return the right result (score 0, empty CIGAR). CI runs the
  tests in release, so it never showed.
- **Root cause:** the CIGAR buffer is `2 * (|P| + |T|)` ops, as in C, and the
  backtraces write backwards, leaving `begin_offset` one below the first
  op. C's offsets are `int`, so -1 is harmless there. In Rust they are
  `usize`, and `max_ops - 1` underflows when the buffer is empty (it only
  grows, so only a fresh aligner hit this).
- **Fix:** `cigar_buffer_len` sizes the buffer with one extra slot, at all
  three sites, which also covers a completely full buffer.
- **Test:** `aligner::tests::test_empty_sequences_cigar` (edit and affine,
  high/med/low, end-to-end and ends-free, BiWFA); fails without the fix.
  The whole suite now also passes in debug.


## 16. End-to-end alignment kept earlier ends-free settings (Rust bug)

- **Symptom:** on an aligner reused after `set_alignment_free_ends` with
  begin-free limits, `align_end2end` (and `align_lambda`,
  `align_packed2bits`) could return a suboptimal score and a CIGAR starting
  with a gap. Found by dada2-rs (HPCBio/dada2-rs#289), which caches one
  aligner per thread: real PacBio pairs scored 3 points below NW's global
  optimum after an ends-free pass. Both compat modes; the CLI and the
  differential harness never mix the two on one aligner.
- **Root cause:** `set_alignment_free_ends` sets the span to ends-free and
  nothing reset it, so end-to-end alignment still seeded the begin-free
  cells (at score `k · |match|`, or 0 with `match == 0`): a leading gap was
  taken as a free skip. Since #60 the dominance pruning was enabled too.
  C has no such state: its span is the aligner's, and `wavefront_align`
  follows it.
- **Fix:** `align_end2end_inner` runs with an end-to-end span and restores
  the previous one, so later `align_endsfree` calls keep their settings.
  Also adds C's missing setter, `set_alignment_end_to_end`, which resets the
  span (and the free-end limits, so `align_endsfree` is end-to-end too).
- **Test:** `aligner::tests::test_end2end_ignores_free_ends` (linear,
  match −5 and 0, both compat modes; fails without the fix). CLI output is
  unchanged.
