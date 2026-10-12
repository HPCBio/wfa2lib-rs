# Performance Log

Benchmarks tied to git commits, as asked for in `PROJECT.md`. Check here before
merging a change that touches the compute, extend or allocation paths, to make
sure it doesn't introduce a regression.

**Method:** `bench/run.sh` (C vs Rust on the same `.seq` file; minimum of two
runs, `Time.Alignment` for Rust and `Time.Benchmark` for C). Dataset:
`bench/generate_data.sh`, which produces 50,000 pairs of 1000 bp at 5% error.
Ratio = Rust / C; below 1.0 means Rust is faster. Also record:
- the machine;
- the C commit;
- max RSS (`/usr/bin/time -l` on macOS, `-v` on Linux) when memory may change.

Run on an otherwise idle machine. Concurrent jobs, or the machine sleeping,
invalidate the numbers.

> `generate_dataset` output may vary between machines and C versions, so ratios
> are only comparable within one entry. Compare commits by re-running both on
> the same data.

---

## Historical (initial port, 2026-03)

Copied from the README and `AVX2_MIGRATION.md`; raw logs weren't kept.

### 2026-03-07 — `2e7ab34` / `3ac6766` (README numbers)

aarch64 (Apple M-series, NEON):

| Benchmark | Rust / C |
|---|---|
| Edit (score) | 0.56x |
| Gap-linear (score) | 0.68x |
| Gap-affine (score) | 0.85x |
| Gap-affine (CIGAR) | 0.90x |
| Gap-affine 2-piece (score) | 0.88x |
| Gap-affine 2-piece (CIGAR) | 0.86x |
| BiWFA (ultralow) | 0.86x |

x86_64 (Intel Broadwell, AVX2):

| Benchmark | Rust | C | Ratio |
|---|---|---|---|
| Edit (score) | 405 ms | 446 ms | 0.91x |
| Gap-linear (score) | 682 ms | 1090 ms | 0.63x |
| Gap-affine (score) | 3980 ms | 3970 ms | 1.00x |
| Gap-affine (CIGAR) | 5440 ms | 5380 ms | 1.01x |
| Gap-affine 2-piece (score) | 16800 ms | 19710 ms | 0.85x |
| Gap-affine 2-piece (CIGAR) | 20620 ms | 22070 ms | 0.93x |
| BiWFA (ultralow) | 7133 ms | 6660 ms | 1.07x |

### 2026-03-06 — `bc838bb` (`AVX2_MIGRATION.md`, NEON)

Edit 0.59x, Linear 0.69x, Affine score 0.94x, Affine CIGAR 0.95x, BiWFA 0.92x,
Affine2p score 1.00x, Affine2p CIGAR 1.03x.

---

## 2026-10 — C-parity fixes (Phase 14)

Machine: Apple M2 Pro (aarch64, NEON), macOS 26.6.2. C: WFA2-lib `3282f83`,
built with `make` minus OpenMP (Apple clang has no `-fopenmp`; strip it from
`tools/align_benchmark/Makefile`).

_Pending: a clean run (idle machine) of `bench/run.sh` for `6d284c8` (main
after #1/#2), `06c593c` (#3) and `29a9c89` (#4), plus max RSS._

Preliminary, informal measurements (not `bench/run.sh`):
- `06c593c` (#3) vs `6d284c8`, 20k × 1 kb / 5%, median of 5: edit and linear
  unchanged (edit 113 vs 115 ms; linear 155 vs 175 ms, within run-to-run noise).
- `29a9c89` (#4) vs `06c593c`: affine/affine2p about 3–8% slower in every mode.
  **Measured with another benchmark running concurrently, so not reliable.**
  To be re-measured. The expected cost comes from restoring the I/D
  `trim_ends` that 2c112b2 removed for speed.
