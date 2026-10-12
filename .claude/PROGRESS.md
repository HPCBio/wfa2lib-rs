# Progress

Status of the Rust port against the goals in `PROJECT.md`: output-identical,
performance-competitive and feature-complete relative to C WFA2-lib. See
`PLAN.md` for phases and next steps, `FIXED_BUGS.md` for bugs, and
`PERFORMANCE_LOG.md` for benchmarks.

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

### Output parity with C (differential testing, fork issue #10)

| Configuration | vs C |
|---|---|
| indel, edit, gap-linear, gap-affine, gap-affine2p; custom penalties; wfa-adaptive | Byte-identical, all memory modes (high / score-only / med / BiWFA / BiWFA score-only) |
| Stock `wfa.utest.sh` | 57/60 identical; 3 are BiWFA CIGAR tie-breaks with match < 0 (#7) |
| banded-static (affine, affine2p, linear) | Byte-identical to C `16969b5b` in high / med / score-only modes, including unreachable pairs |
| banded-static (edit) | Rust stops as unreachable; C never terminates (C bug) |
| banded-adaptive | Differs: Rust applies the static clamp every N steps instead of C's re-centred band (#6) |
| xdrop / zdrop | Not applied in end-to-end mode (#5) |
| BiWFA + any heuristic | Heuristics ignored in BiWFA sub-aligners (#6) |
| ends-free / extension | Not reachable from the CLI (`--wfa-span` unused, #8); `match == 0` begin-free seeding missing (#9) |

### Upstream C sync

C changes after the port baseline that affect behaviour (fork issue #10):

- [ ] 6565731 / 7ea485e (C #125): band inherited into BiWFA sub-aligners, and
  compute limits clamped to `[min_k, max_k]`. This also changes unidirectional
  banded output (#6).
- [ ] 90582eb (C #126): don't apply ends-free begin seeding in end-to-end mode
  (#9).
- n/a 061509a SAM CIGAR table (no SAM output in Rust); build, packaging and
  C++-binding commits.

### Upstream Rust (COMBINE-lab/wfa2lib-rs) PRs

- COMBINE-lab #6: bugs 1–2 (bounds check, `init_ends`). Open.
- Planned: bugs 3–4 combined into one PR (fork #12 + #14).

## Log

- **2026-03-04..09:** initial port (phases 1–12) by Rob Patro with Claude Code.
- **2026-06-18:** `cli` feature gate, `set_alignment_scope` setter, license
  metadata (COMBINE-lab #2–#4).
- **2026-10-12:**
  - Identified the port baseline (above).
  - Built the differential test harness (Rust vs C HEAD vs C baseline, 28
    configs × 5 modes × 2 datasets) and filed fork issues #1–#10.
  - Fixed bugs 1–2 (fork PR #11; COMBINE-lab #6), bug 3 (fork PR #12) and
    bug 4 (fork PR #14), all merged into the fork's `main`.
  - Added these tracking files.
