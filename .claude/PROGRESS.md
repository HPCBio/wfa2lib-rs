# Progress

Status of the Rust port against the goals in `PROJECT.md`: output-identical,
performance-competitive and feature-complete relative to C WFA2-lib. See
`PLAN.md` for phases and next steps, `FIXED_BUGS.md` for bugs, and
`PERFORMANCE_LOG.md` for benchmarks.

> These files were started in the HPCBio/wfa2lib-rs fork. Commit hashes,
> statuses and issue/PR numbers refer to that fork unless qualified otherwise.

## Port baseline

Nothing in the repository records which C commit the port came from. Fix
fingerprints in the code place it on C `main` between `ac9f2ab5` (2026-02-08,
static-band backtrace fixes present) and **`16969b5b` (2026-02-28)**, the head
when the port began (2026-03-04). Algorithmically, `16969b5b` is the same as
the v2.3.6 tag. Changes in C after that date are tracked under "Upstream C
sync" below.

## Status (as of 2026-10-12)

| Area | Status |
|---|---|
| Phases 1–11 (types → CLI) | Done (initial port, 2026-03-04) |
| Phase 12 (performance) | Done: NEON/AVX2 kernels, arenas, slab reuse (2026-03-05..09) |
| Phase 13 (polish) | Partial: tracking files added 2026-10-12; plotting, rustdoc pass, examples and a full `cargo fmt` pass outstanding |

### Output parity with C (differential testing, HPCBio/wfa2lib-rs#10)

Checked with `tests/differential/run.sh` (1262 cases against the C commit in
`tests/differential/C_COMMIT`; accepted differences in `known_diffs.tsv`).

| Configuration | vs C |
|---|---|
| indel, edit, gap-linear, gap-affine, gap-affine2p; custom penalties; wfa-adaptive | Byte-identical, all memory modes (high / score-only / med / BiWFA / BiWFA score-only) |
| Stock `wfa.utest.sh` | 60/60 identical (BiWFA match < 0 tie-breaks fixed, HPCBio/wfa2lib-rs#7) |
| banded-static (affine, affine2p, linear) | Byte-identical to C 3282f83 (with Marco-Sola-Lab/WFA2-lib#125's band clamp) in high / med / score-only modes, including unreachable pairs |
| banded-static (edit, indel) | Rust stops as unreachable; C never terminates (C bug, HPCBio/wfa2lib-rs#24) |
| banded-adaptive | C's re-centred band ported (HPCBio/wfa2lib-rs#6): byte-identical to C 3282f83 in high / score-only / med / low (48/48) |
| xdrop / zdrop (affine, linear, affine2p; global / extension / ends-free) | Byte-identical to C 3282f83 in high / score-only / med (556/576 matrix cells). The other 20 are Z-drop `med` with match < 0, where C reads uninitialized memory (HPCBio/wfa2lib-rs#22). Includes C bugs kept for parity: HPCBio/wfa2lib-rs#21, HPCBio/wfa2lib-rs#23 |
| BiWFA + any heuristic | Byte-identical to C 3282f83 (168/168: banded-static/adaptive, xdrop, zdrop, wfa-adaptive; CIGAR and score-only; HPCBio/wfa2lib-rs#6) |
| `--wfa-memory=med\|low` (piggyback backtrace) | Byte-identical to C 3282f83 for all metrics (28/28 matrix cells, 50k × 1 kb), with heuristics as above (HPCBio/wfa2lib-rs#8) |
| extension; ends-free with end gaps only | Byte-identical to C (all metrics incl. match < 0, alignment and score-only) |
| ends-free with begin gaps (P0/T0 > 0) | Byte-identical to C, match == 0 and match < 0, alignment / score-only / med (HPCBio/wfa2lib-rs#9). Edit/indel CIGARs valid since C 6176f94 (Marco-Sola-Lab/WFA2-lib#144, HPCBio/wfa2lib-rs#18) |

### Upstream C sync

C changes after the port baseline that affect behaviour (HPCBio/wfa2lib-rs#10):

- [x] 6565731 / 7ea485e (Marco-Sola-Lab/WFA2-lib#125): band inherited into BiWFA sub-aligners, and
  compute limits clamped to `[min_k, max_k]`. This also changes unidirectional
  banded output (HPCBio/wfa2lib-rs#6, `3b16283` and the BiWFA commit).
- [x] 90582eb (Marco-Sola-Lab/WFA2-lib#126): ends-free begin seeding only for ends-free span. Implemented this way in HPCBio/wfa2lib-rs#9; the baseline's end-to-end widening was never ported
  (HPCBio/wfa2lib-rs#9).
- [x] 11b225e / 6176f94 (Marco-Sola-Lab/WFA2-lib#144): M0's initialized range
  covers the ends-free begin seeds, so edit/indel begin-gap CIGARs are valid
  (HPCBio/wfa2lib-rs#18). Was corrected-mode only; now unconditional.
- [x] cda19a0 (Marco-Sola-Lab/WFA2-lib#140): the slab reaps only the wavefronts
  handed out since the last reap. No output change. Ported in
  `WavefrontSlab::reap_cigar` (busy list, recorded only while the slab is
  large): after long alignments, short ones in high memory were up to 3.5×
  slower than needed (PERFORMANCE_LOG 2026-10-17).
- n/a 061509a SAM CIGAR table (no SAM output in Rust); build, packaging and
  C++-binding commits.

### Known issues kept for C parity

- HPCBio/wfa2lib-rs#21: X-drop with `steps_between_cutoffs > 1` returns invalid
  CIGARs (C's high-memory backtrace reads wavefronts trimmed to empty; C's
  `med` is valid). **Fixed in corrected mode** (`Compat::Corrected`;
  FIXED_BUGS 14d).
- HPCBio/wfa2lib-rs#22: Z-drop results differ between C's high-memory and `med`
  modes. Root cause: in `med`/`low`, C backtraces from the end diagonal at
  the current score, often a dead cell whose chain comes from uninitialized
  memory, so C's output depends on earlier alignments. Rust reads such
  cells as zeroes (deterministic; matches C on fresh memory). C's high mode
  backtraces at the same wrong score (empty or invalid CIGARs). **Fixed in
  corrected mode** (FIXED_BUGS 14e).
- HPCBio/wfa2lib-rs#55: X-drop/Z-drop with `match == 0` use a sign-flipped
  Smith-Waterman score (C regression from 51fb1aa): X-drop abandons most
  divergent pairs, Z-drop's best cell is the start. **Fixed in corrected
  mode** (FIXED_BUGS 14f).
- HPCBio/wfa2lib-rs#23: score-only score for unreachable alignments (match < 0)
  is evaluated at C's unset end position (signed overflow); reproduced.
  **Fixed in corrected mode** in the CLI, which prints -1 for every
  score-only result without an alignment (FIXED_BUGS 14h).
- HPCBio/wfa2lib-rs#24: edit/indel + banded-static never terminates in C;
  Rust deliberately terminates (divergence).
- HPCBio/wfa2lib-rs#34: C's ends-free/extension results depend on the CPU.
  Built with `-march=native` on AVX2/AVX-512, C's ends-free extend checks
  termination in a different diagonal order from its scalar kernel, so it
  can pick a different one of several equally optimal alignments. The
  ends-free score is the same and optimal (checked against an exact DP);
  `align_benchmark`'s printed score differs only because it rescores the
  CIGAR globally, charging free end gaps. Rust follows the scalar order on
  every CPU; the harness builds C without `-march=native`.
- Marco-Sola-Lab/WFA2-lib#102, HPCBio/wfa2lib-rs#37, HPCBio/wfa2lib-rs#38: ends-free
  with match < 0 can be suboptimal (trailing-end miss, repeat-unit shift,
  all-gap alignment with whole-length free ends). This showed up in
  HPCBio/dada2-rs#51. Note that dada2-rs's "free leading end-gap not
  credited" failure mode was the Rust HPCBio/wfa2lib-rs#9 bug, fixed by
  HPCBio/wfa2lib-rs#17; honoring begin gaps as C does then exposed #38.
  **Fixed in corrected mode** (`Compat::Corrected`, `--no-c-compat`;
  HPCBio/wfa2lib-rs#32), optimal against a DP oracle (FIXED_BUGS 14a). The
  default stays C-identical.

### Upstream contributions

All upstream work (COMBINE-lab/wfa2lib-rs PRs grouped by dependency, C
reports, inbound C ports) is tracked in HPCBio/wfa2lib-rs#19.

- COMBINE-lab/wfa2lib-rs#6: bugs 1–2 (bounds check, `init_ends`). Open.
- Queued: bugs 3–4 (HPCBio/wfa2lib-rs#12 + HPCBio/wfa2lib-rs#14), then the
  CLI changes and ends-free begin gaps.
- Marco-Sola-Lab/WFA2-lib#142: edit/indel begin-gap invalid CIGAR. Fixed upstream
  by Marco-Sola-Lab/WFA2-lib#144 (merged 2026-10-16, 6176f94); ported.

## Log

- **2026-03-04..09:** initial port (phases 1–12) by Rob Patro with Claude Code.
- **2026-06-18:** `cli` feature gate, `set_alignment_scope` setter, license
  metadata (COMBINE-lab/wfa2lib-rs#2, COMBINE-lab/wfa2lib-rs#3, COMBINE-lab/wfa2lib-rs#4).
- **2026-10-12:**
  - Identified the port baseline (above).
  - Built the differential test harness (Rust vs C HEAD vs C baseline, 28
    configs × 5 modes × 2 datasets) and filed HPCBio/wfa2lib-rs#1–HPCBio/wfa2lib-rs#10.
  - Fixed bugs 1–2 (HPCBio/wfa2lib-rs#11; COMBINE-lab/wfa2lib-rs#6), bug 3 (HPCBio/wfa2lib-rs#12) and
    bug 4 (HPCBio/wfa2lib-rs#14), all merged into the fork's `main`.
  - Added these tracking files.
- **2026-10-13:** HPCBio/wfa2lib-rs#6 (banded-adaptive as in C,
  Marco-Sola-Lab/WFA2-lib#125 band clamp and BiWFA band shift, BiWFA heuristics and
  statuses, BiWFA score-only) and HPCBio/wfa2lib-rs#7 (BiWFA length fallback at
  the top level only). BiWFA × every heuristic is now byte-identical to C
  3282f83, and the stock `wfa.utest.sh` is 60/60.
- **2026-10-13:** HPCBio/wfa2lib-rs#8 `--wfa-memory=med|low`: ported C's
  piggyback backtrace (kernels, offloading, backtrace buffer, unpacking) and
  buffer compaction. Root cause of HPCBio/wfa2lib-rs#22 found.
- **2026-10-13:** committed the differential harness (`tests/differential/`):
  C built at a pinned commit, 1262 cases over metrics, spans, heuristics and
  memory modes, with known differences listed (HPCBio/wfa2lib-rs#10). CI
  proposed in HPCBio/wfa2lib-rs#29.
- **2026-10-13:** CI (HPCBio/wfa2lib-rs#29): GitHub Actions runs clippy,
  the tests, a library-only build, `wfa_utest.sh` and the differential
  harness on Linux (x86-64) and macOS (arm64), on every PR and push to
  `main`, so Rust is checked to be deterministic across platforms; Linux
  also runs with `target-cpu=x86-64-v3` (AVX2 enabled at compile time). CI
  overrides `.cargo/config.toml`'s `target-cpu=native`, since cached builds
  from one runner's CPU can crash on another. Its
  first run on x86 found HPCBio/wfa2lib-rs#34 (C's AVX2 ends-free kernel)
  and x86-only clippy errors.
- **2026-10-13:** HPCBio/wfa2lib-rs#31 high-memory allocation: wavefronts
  allocated over C's growing range and a "tight" slab in CIGAR mode, as C.
  Long pairs now match C's memory and time (utest affine 927 → 275 MB,
  affine2p 3.5 → 1.8 GB; ends-free affine2p match < 0 18.8 s → 3.4 s vs
  C 4.2 s); 1 kb benchmark within ±3%. See PERFORMANCE_LOG.
- **2026-10-14:** two Rust bugs found by a DP-oracle sweep: the `med`/`low`
  begin-free seed backtrace (match == 0, invalid CIGAR; FIXED_BUGS 12) and
  the affine2p modular-slot leak, which could segfault (FIXED_BUGS 13).
- **2026-10-14:** HPCBio/wfa2lib-rs#32 per-aligner `Compat` setting
  (default `C`; `Corrected` fixes known C bugs, CLI `--no-c-compat`). Covers
  ends-free with match < 0 (HPCBio/wfa2lib-rs#37, HPCBio/wfa2lib-rs#38,
  Marco-Sola-Lab/WFA2-lib#102), now optimal, and edit/indel begin-gap CIGARs
  (HPCBio/wfa2lib-rs#18). Checked against a DP oracle.
- **2026-10-14:** corrected-mode speed (HPCBio/wfa2lib-rs#31): its extend
  kernel was the overhead (+26–28% instructions over C mode on `long.seq`,
  now +2–2.5%). Also exact pruning: with an end candidate, wavefronts are
  trimmed to the diagonals that could still beat it (output unchanged).
  Null steps no longer end a corrected search while a candidate is pending
  (FIXED_BUGS 14c, with `BandedStatic`).
  Also HPCBio/wfa2lib-rs#41: empty pattern + text no longer panics in debug
  (FIXED_BUGS 15). CI now also runs the tests in a debug build (Linux and
  macOS), so overflows like it fail CI.
- **2026-10-15:** Rust-vs-C sweep (HPCBio/wfa2lib-rs#31): 301 configs, all
  output identical to C; Rust faster or at parity almost everywhere. Gaps
  (med/low kernels, indel med/low, high-memory CIGAR RSS, BiWFA sequence
  copies) are now a checklist in HPCBio/wfa2lib-rs#31. `bench/sweep.py`
  makes the sweep reproducible on other machines (C also built with
  `-march=native`).
- **2026-10-15:** med/low piggyback kernels vectorized and fused
  (HPCBio/wfa2lib-rs#31): none vectorized before (LLVM turned selects
  between loads into gathers). med/low now 0.52–0.83× its previous time and
  faster than C on ARM; output unchanged.
- **2026-10-15:** x86 sweeps (`bench/sweep.py` on the HPCBio cluster): Rust
  up to 2.16× slower than C on long/divergent edit/indel/linear, traced with
  `perf` to AVX2 gather stalls in the extend kernel; the gather kernels are
  removed (HPCBio/wfa2lib-rs#31).
- **2026-10-15:** BiWFA overlap scan vectorized (HPCBio/wfa2lib-rs#31): the
  breakpoint checks were 73% of affine2p BiWFA time and scalar except the
  NEON M scan. affine2p ultralow is now 0.38–0.52× C on ARM (was ~1.0×);
  output unchanged.
- **2026-10-17:** x86 prefetches (HPCBio/wfa2lib-rs#31): on the Xeon
  (Skylake-SP), affine2p high-memory CIGAR waited on ownership requests for
  its fresh output lines (store buffer full). The input prefetches are removed
  and the affine2p kernel issues `PREFETCHW` on its outputs. Xeon: 5–35%
  faster on affine2p high memory, now 0.83–0.89× C; EPYC slightly faster;
  output unchanged.

- **2026-10-17:** C `C_COMMIT` 3282f83 → 6176f94. Marco-Sola-Lab/WFA2-lib#144
  (our fix for Marco-Sola-Lab/WFA2-lib#142) merged upstream, so the edit/indel
  begin-gap fix (HPCBio/wfa2lib-rs#18) now applies in default C mode too.
- **2026-10-17:** HPCBio/wfa2lib-rs#21 fixed in corrected mode: a wavefront
  trimmed empty by a heuristic cutoff is null for the backtrace too, as for
  compute. X-drop with steps > 1 now gives valid CIGARs; scores unchanged.
- **2026-10-17:** HPCBio/wfa2lib-rs#22 fixed in corrected mode: a Z-drop end
  is backtraced from the best cell at the score it was found at (its
  backtrace kept and pinned for med/low). All memory modes agree and every
  partial CIGAR is valid; C's high mode had invalid ones with match < 0.
- **2026-10-17:** HPCBio/wfa2lib-rs#55 (found while fixing #22): with
  `match == 0`, C's X-drop/Z-drop score has its sign flipped since 51fb1aa.
  Corrected mode uses +1: X-drop with default penalties goes from ~0/300 to
  300/300 optimal on `long.seq`.
- **2026-10-17:** corrected mode now reaches BiWFA's sub-aligners (a Rust bug
  from #54, FIXED_BUGS 14g): the #21 and #55 fixes apply under BiWFA too.
- **2026-10-17:** HPCBio/wfa2lib-rs#23 fixed in corrected mode (CLI): score-only
  output without an alignment is -1, as in CIGAR mode, instead of C's
  overflowed, give-up or INT32_MIN scores.
- **2026-10-17:** ported C's busy-list slab reap (Marco-Sola-Lab/WFA2-lib#140):
  mixed long/short workloads in high memory 0.30× the time (C's benchmark:
  1231 → 365 ms, C 421 ms); a dada2-like interleaved workload 735 → 209 ms.
  Output identical.
- **2026-10-18:** begin-free ends-free cost (HPCBio/dada2-rs#51) investigated:
  inherent step count, not prunable exactly; two per-step overhead cuts
  (short NULL fills, banded cutoff fast path) give 9–16% on dada2-like pairs.
