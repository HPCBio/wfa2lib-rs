# Differential tests against C WFA2-lib

The Rust port aims to be output-identical to C WFA2-lib. This harness runs
C's and Rust's `align_benchmark` over a matrix of configurations and
compares their outputs byte for byte. It covers much more than the stock
`wfa.utest.sh`, which exercises no heuristics, spans or memory modes.

```bash
tests/differential/run.sh              # all suites (~1300 cases)
tests/differential/run.sh drops biwfa  # selected suites
JOBS=4 TIMEOUT=60 tests/differential/run.sh core
```

`run.sh` builds Rust (`cargo build --release`) and C (`build_c.sh`), runs
every case, and writes `target/differential/results.tsv` (suite, dataset,
arguments, result). It exits non-zero if any case differs from C and isn't
listed in `known_diffs.tsv`. It also reports listed differences that now
match.

## Files

| File | Purpose |
|---|---|
| `C_COMMIT` | The pinned C commit (the parity target). |
| `build_c.sh [COMMIT]` | Clones and builds C `align_benchmark` into `target/differential/c-<commit>/`. Source: `$WFA2_LIB_REPO`, else `./WFA2-lib`, else GitHub. It strips `-march=native`, so C uses its scalar kernels everywhere (its AVX2 ends-free kernel gives CPU-dependent results: HPCBio/wfa2lib-rs#34), and on macOS `-fopenmp` (Apple clang has no OpenMP). On macOS it also makes C's benchmark timer use `clock_gettime`, as on Linux: C's Mach clock calls cost about 1.3 µs per alignment inside the timed interval, most of the measured time for short pairs. Cached builds are redone when the build recipe changes (`RECIPE`). |
| `cases.sh` | Generates the cases: `suite <TAB> dataset <TAB> arguments`. |
| `run.sh` | Runs and compares the cases (in parallel). |
| `known_diffs.tsv` | Accepted differences, each with the issue explaining it. |
| `data/long.seq.gz` | 300 × 2 kb pairs, 8% error plus long indels (C `generate_dataset -n 300 -l 2000 -e 0.08 --indels 3,40`). Committed because `generate_dataset` seeds from the clock. |

The other dataset is C's `tests/wfa.utest.seq`, taken from the C checkout.

## Suites

| Suite | Covers |
|---|---|
| `core` | edit, indel, gap-linear, gap-affine (incl. match < 0), gap-affine2p, wfa-adaptive × high / score-only / med / low / BiWFA / BiWFA score-only |
| `span` | extension and ends-free (begin and/or end gaps, fractional and absolute), incl. match < 0 |
| `banded` | banded-static (gap metrics; edit/indel excluded, C never terminates: HPCBio/wfa2lib-rs#24) |
| `adaptive` | banded-adaptive, incl. med/low |
| `drops` | X-drop / Z-drop × global / extension / ends-free × high / score-only / med |
| `biwfa` | BiWFA with every heuristic |

Results are `identical`, `same-error` (both exit with the same error, e.g.
options C rejects), or `differ` (with counts of differing scores and
CIGARs).

## Changing the C commit

To move to a newer C commit:

1. Run `build_c.sh <commit>` and `C_BIN=... run.sh`.
2. Port or document the differences.
3. Update `C_COMMIT` and `known_diffs.tsv`.
4. Record the sync point in `.claude/PROGRESS.md` ("Upstream C sync").
