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

| Configuration | vs C |
|---|---|
| indel, edit, gap-linear, gap-affine, gap-affine2p; custom penalties; wfa-adaptive | Byte-identical, all memory modes (high / score-only / med / BiWFA / BiWFA score-only) |
| Stock `wfa.utest.sh` | 60/60 identical (BiWFA match < 0 tie-breaks fixed, HPCBio/wfa2lib-rs#7) |
| banded-static (affine, affine2p, linear) | Byte-identical to C 3282f83 (with smarco/WFA2-lib#125's band clamp) in high / med / score-only modes, including unreachable pairs |
| banded-static (edit, indel) | Rust stops as unreachable; C never terminates (C bug, HPCBio/wfa2lib-rs#24) |
| banded-adaptive | C's re-centred band ported (HPCBio/wfa2lib-rs#6): 34/36 matrix cells byte-identical to C 3282f83; the other 2 are `med` (C piggyback, HPCBio/wfa2lib-rs#8) |
| xdrop / zdrop (affine, linear, affine2p; global / extension / ends-free) | Byte-identical to C 3282f83 in high and score-only modes (450/576 matrix cells; the other 126 are `med`, where C uses its piggyback backtrace, HPCBio/wfa2lib-rs#8 / HPCBio/wfa2lib-rs#22). Includes C bugs kept for parity: HPCBio/wfa2lib-rs#21, HPCBio/wfa2lib-rs#23 |
| BiWFA + any heuristic | Byte-identical to C 3282f83 (168/168: banded-static/adaptive, xdrop, zdrop, wfa-adaptive; CIGAR and score-only; HPCBio/wfa2lib-rs#6) |
| extension; ends-free with end gaps only | Byte-identical to C (all metrics incl. match < 0, alignment and score-only) |
| ends-free with begin gaps (P0/T0 > 0) | Byte-identical to C, match == 0 and match < 0, alignment / score-only / med (HPCBio/wfa2lib-rs#9). Note: for edit/indel this includes a C bug that can emit an invalid CIGAR (known issue HPCBio/wfa2lib-rs#18, from C `380eb31b`) |

### Upstream C sync

C changes after the port baseline that affect behaviour (HPCBio/wfa2lib-rs#10):

- [x] 6565731 / 7ea485e (smarco/WFA2-lib#125): band inherited into BiWFA sub-aligners, and
  compute limits clamped to `[min_k, max_k]`. This also changes unidirectional
  banded output (HPCBio/wfa2lib-rs#6, `3b16283` and the BiWFA commit).
- [x] 90582eb (smarco/WFA2-lib#126): ends-free begin seeding only for ends-free span. Implemented this way in HPCBio/wfa2lib-rs#9; the baseline's end-to-end widening was never ported
  (HPCBio/wfa2lib-rs#9).
- n/a 061509a SAM CIGAR table (no SAM output in Rust); build, packaging and
  C++-binding commits.

### Known issues kept for C parity

- HPCBio/wfa2lib-rs#18: edit/indel ends-free with begin gaps can return an invalid CIGAR
  (C regression from `380eb31b`; reported as smarco/WFA2-lib#142 with a
  proposed fix).
- HPCBio/wfa2lib-rs#21: X-drop with `steps_between_cutoffs > 1` returns invalid
  CIGARs (C high-memory backtrace over pruned diagonals; C's `med` is valid).
- HPCBio/wfa2lib-rs#22: Z-drop results differ between C's high-memory and `med`
  modes (Rust matches high; root cause open).
- HPCBio/wfa2lib-rs#23: score-only score for unreachable alignments (match < 0)
  is evaluated at C's unset end position (signed overflow); reproduced.
- HPCBio/wfa2lib-rs#24: edit/indel + banded-static never terminates in C;
  Rust deliberately terminates (divergence).
- smarco/WFA2-lib#102: ends-free with match < 0 can be suboptimal, depending on argument
  order. This showed up in HPCBio/dada2-rs#51. Note that dada2-rs's
  "free leading end-gap not credited" failure mode was the Rust HPCBio/wfa2lib-rs#9 bug,
  fixed by HPCBio/wfa2lib-rs#17.

### Upstream contributions

All upstream work (COMBINE-lab/wfa2lib-rs PRs grouped by dependency, C
reports, inbound C ports) is tracked in HPCBio/wfa2lib-rs#19.

- COMBINE-lab/wfa2lib-rs#6: bugs 1–2 (bounds check, `init_ends`). Open.
- Queued: bugs 3–4 (HPCBio/wfa2lib-rs#12 + HPCBio/wfa2lib-rs#14), then the
  CLI changes and ends-free begin gaps.
- smarco/WFA2-lib#142: edit/indel begin-gap invalid CIGAR. Reported.

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
  smarco/WFA2-lib#125 band clamp and BiWFA band shift, BiWFA heuristics and
  statuses, BiWFA score-only) and HPCBio/wfa2lib-rs#7 (BiWFA length fallback at
  the top level only). BiWFA × every heuristic is now byte-identical to C
  3282f83, and the stock `wfa.utest.sh` is 60/60.
