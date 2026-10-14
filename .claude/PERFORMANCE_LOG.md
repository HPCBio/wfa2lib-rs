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

Data: `bench/generate_data.sh` settings (50,000 × 1000 bp, 5% error),
regenerated with C `3282f83`'s `generate_dataset`. Timing as in `bench/run.sh`
(minimum of 2 runs). Max RSS from `/usr/bin/time -l` (the larger of the 2 runs).
The machine was otherwise idle apart from an endpoint-backup agent (~0.5 core).
Within each mode, C and all three commits ran back to back.

### 2026-10-12 — `6d284c8` (main, HPCBio/wfa2lib-rs#1 and HPCBio/wfa2lib-rs#2), `06c593c` (HPCBio/wfa2lib-rs#3), `29a9c89` (HPCBio/wfa2lib-rs#4)

| Mode | C (ms) | `6d284c8` | `06c593c` | `29a9c89` | 29a9c89 vs 06c593c |
|---|---|---|---|---|---|
| Edit score | 403 | 223 (0.55x) | 223 (0.55x) | 228 (0.57x) | +2.0% (noise) |
| Linear score | 555 | 369 (0.66x) | 328 (0.59x) | 323 (0.58x) | -1.3% |
| Affine score | 2370 | 1881 (0.79x) | 1887 (0.80x) | 1940 (0.82x) | +2.8% |
| Affine CIGAR | 2610 | 1938 (0.74x) | 1973 (0.76x) | 2074 (0.79x) | +5.2% |
| BiWFA (ultralow) | 4110 | 3106 (0.76x) | 3079 (0.75x) | 3181 (0.77x) | +3.3% |
| Affine2p score | 8430 | 6975 (0.83x) | 6919 (0.82x) | 7580 (0.90x) | +9.6% |
| Affine2p CIGAR | 8910 | 7433 (0.83x) | 7462 (0.84x) | 7899 (0.89x) | +5.9% |

- **`06c593c` (HPCBio/wfa2lib-rs#3):** no regression.
- **`29a9c89` (HPCBio/wfa2lib-rs#4):** costs 3–10% on affine/affine2p. Rust is still faster than
  C in every mode.

**Where `29a9c89`'s cost comes from** (minimum of 3 runs, the two most affected modes):

| Commit | Affine2p score | Affine CIGAR |
|---|---|---|
| `06c593c` (before `29a9c89`) | 7016 | 1969 |
| `320611d` cutoff equates I/D only when M narrowed | 6993 (±0%) | 1973 (±0%) |
| `ddae4d8` I/D `trim_ends` restored | 7388 (+5.3%) | 2022 (+2.7%) |
| `29a9c89` null-input handling | 7419 (+0.4%) | 2041 (+1.0%) |

Nearly all of it comes from restoring C's independent I/D trimming, which
2c112b2 removed for speed (C pays the same cost). That removal changed results,
and output-identical comes first in the project goals, so the trimming stays.
**Optimization target:** make I/D trimming cheaper without changing which
diagonals survive.

**Memory** (since fixed by streaming CLI I/O, `0f1eed5`; see below): Rust max RSS is about 200 MB in every mode, at all three commits,
vs 1.6–9.5 MB for C. That points to the `align_benchmark` CLI rather than the
library: it reads the whole input file (about 100 MB here) into memory and
buffers all output lines, while C streams. It fails PLAN's "memory within 1.5x"
target as measured. Tracked as a CLI gap (HPCBio/wfa2lib-rs#8). Library-level
memory should be measured separately, without CLI buffering.

### 2026-10-12 — `6a06d8f` (HPCBio/wfa2lib-rs#9 ends-free begin gaps) vs `eee854b` (main)

End-to-end, same machine and data, Rust only (minimum of 3 runs, ms):

| Mode | `eee854b` | `6a06d8f` |
|---|---|---|
| Edit score | 220 | 201 |
| Linear score | 321 | 324 |
| Affine score | 1939 | 1937 |
| Affine CIGAR | 2074 | 2070 |
| Affine2p score | 7420 | 7416 |

No regression from the per-step begin-free check. Edit is about 8% faster,
plausibly because the edit step now writes C's four boundary NULLs directly
instead of calling `init_ends`.

### 2026-10-13 — `fix/xdrop-zdrop` vs `3488e33` (main)

Same machine and data, Rust only, minimum of 3 runs (ms). This covers the
ported affine backtrace, null-step counting in compute, and the
`2*(plen+tlen)` CIGAR buffer:

| Mode | `3488e33` | branch | Δ |
|---|---|---|---|
| Edit score | 202 | 204 | +0.9% |
| Linear score | 325 | 328 | +1.0% |
| Affine score | 1907 | 1905 | -0.1% |
| Affine CIGAR | 2031 | 2046 | +0.7% |
| BiWFA (ultralow) | 3224 | 3256 | +1.0% |
| Affine2p score | 7532 | 7522 | -0.1% |
| Affine2p CIGAR | 7794 | 7952 | +2.0% |

No measurable regression (within ~2%). The largest change, affine2p CIGAR,
is consistent with the ported backtrace doing all of C's source lookups.

### 2026-10-13 — `0f1eed5` streaming CLI I/O: memory

Max RSS (`/usr/bin/time -l`), 50,000 × 1 kb, same machine. "main" is
`55231da`.

| Mode | C | main | `0f1eed5` | `0f1eed5`, system allocator |
|---|---|---|---|---|
| Edit score | 1.6 MB | 199.2 MB | 2.5 MB | 2.0 MB |
| Affine score | 1.7 MB | 199.2 MB | 2.6 MB | 2.1 MB |
| Affine CIGAR | 4.3 MB | 199.2 MB | 10.7 MB | 10.0 MB |
| BiWFA (ultralow) | 3.2 MB | 199.2 MB | 6.2 MB | 5.7 MB |

- The ~200 MB was the CLI buffering the whole input and all output. Streaming
  removes it; runtime is unchanged (affine score 1944 → 1932 ms).
- The library's own memory is now visible: about 1.2–2.5× C. mimalloc accounts
  for only ~0.5 MB.
- **Next memory target:** in CIGAR (non-modular) mode every wavefront is
  allocated over the full historic diagonal range, whereas C allocates each at
  its `lo..hi` (plus margins). That likely explains most of the affine CIGAR
  gap (10.7 vs 4.3 MB) against PLAN's "memory within 1.5x" target.
  Addressed by `perf/high-memory` (HPCBio/wfa2lib-rs#31; see below).


### 2026-10-13 — `feat/biwfa-heuristics-banded` vs `6031aa7` (main)

Same machine and data, Rust only, minimum of 3 runs (ms). Covers
banded-adaptive, smarco/WFA2-lib#125's band clamp with 1000 × 1000 initial
components, BiWFA heuristic/status handling, and the top-level-only BiWFA length
fallback (HPCBio/wfa2lib-rs#6, HPCBio/wfa2lib-rs#7):

| Mode | `6031aa7` | branch | Δ |
|---|---|---|---|
| Edit score | 200 | 199 | -0.5% |
| Affine score | 1919 | 1930 | +0.6% |
| Affine CIGAR | 2071 | 2082 | +0.5% |
| BiWFA (ultralow) | 3313 | 3357 | +1.3% |
| Affine2p score | 7636 | 7732 | +1.3% |

No measurable regression (within ~1.5%). Recursing to deeper BiWFA levels
instead of falling back at 100 bp costs about 1% on BiWFA.

### 2026-10-13 — `feat/piggyback-backtrace` vs `b1d1101` (main)

Same machine and data, Rust only, minimum of 3 runs (ms). Max RSS from
`/usr/bin/time -l`.

| Mode | `b1d1101` | branch | Δ |
|---|---|---|---|
| Edit score | 199 | 204 | +2.5% |
| Affine score | 1905 | 1940 | +1.8% |
| Affine CIGAR | 2063 | 2083 | +1.0% |
| BiWFA (ultralow) | 3346 | 3481 | +4.0% |
| Affine2p score | 7554 | 7662 | +1.4% |
| Affine2p CIGAR | 7970 | 8067 | +1.2% |
| Affine `med` | 2067 (= high) | 8562 | new |
| Affine `low` | 2079 (= high) | 10342 | new |

- `med`/`low` now do what C does, and are faster than C (affine `med` 8.3 s
  vs C 11.0 s, `low` 10.0 vs 12.1, edit `med` 0.53 vs 1.18). Affine `med`
  max RSS is 3.6 MB (`high` 11.3, C `med` 2.1).
- The other modes are 1–2.5% slower (4% BiWFA, measured before trimming
  per-step checks). The cost is per score step: piggyback branches and
  written-range tracking. Interleaved re-runs after making those checks
  cheaper: affine score +1.8%, CIGAR +1.0%.
- **Kernel note:** a fused raw-pointer piggyback loop (C's structure) ran
  at 40+ ns per cell on M2 from false store/load aliasing between
  identically aligned arrays. Per-component passes over slices
  (`compute/piggyback.rs`) vectorize: ~4 ns per cell.

### 2026-10-13 — `perf/high-memory` vs `22d97e8` (main): C's wavefront allocation

HPCBio/wfa2lib-rs#31. Same machine (M2 Pro), idle. Time is the minimum of 3
interleaved runs; max RSS from `/usr/bin/time -l`.

**Long pairs** (`wfa.utest.seq`, incl. 10 kb pairs; `long.seq`, 300 × 2 kb),
wall time:

| Configuration | main | branch | C | main RSS | branch RSS | C RSS |
|---|---|---|---|---|---|---|
| utest affine2p `-1,4,6,2,24,1`, ends-free `0.1,0,0.1,0` | 18.76 s | 3.43 s | 4.19 s | 4838 MB | 3665 MB | 3574 MB |
| utest affine2p | 2.54 s | 1.67 s | 1.69 s | 3529 MB | 1789 MB | 1753 MB |
| utest affine | 0.41 s | 0.35 s | 0.35 s | 927 MB | 275 MB | 272 MB |
| utest edit | 0.04 s | 0.03 s | 0.03 s | 56 MB | 15 MB | 22 MB |
| utest gap-linear | 0.05 s | 0.04 s | 0.04 s | 85 MB | 27 MB | 25 MB |
| long affine | 0.12 s | 0.11 s | 0.11 s | 42 MB | 17 MB | 11 MB |
| long affine2p | 0.60 s | 0.53 s | 0.53 s | 139 MB | 62 MB | 49 MB |

**Standard benchmark** (50,000 × 1 kb, 5% error), ms:

| Mode | main | branch | Δ | C |
|---|---|---|---|---|
| Edit score | 204 | 209 | +2.6% | 401 |
| Edit CIGAR | 251 | 258 | +2.5% | 750 |
| Linear score | 339 | 338 | −0.5% | 570 |
| Linear CIGAR | 390 | 383 | −1.6% | 962 |
| Affine score | 1948 | 1944 | −0.2% | 2300 |
| Affine CIGAR | 2098 | 2079 | −0.9% | 2660 |
| Affine `med` | 8475 | 8196 | −3.3% | 10540 |
| BiWFA (ultralow) | 3327 | 3320 | −0.2% | 4180 |
| Affine2p score | 7543 | 7663 | +1.6% | 8570 |
| Affine2p CIGAR | 7892 | 8051 | +2.0% | 8930 |

Max RSS on this set: affine CIGAR 10 → 6 MB (C 4), affine2p CIGAR 25 → 13 MB
(C 9). Run-to-run noise (same binary twice) is about 1%.

What changed (all as C):
- **Allocation range.** Every wavefront was allocated over the full diagonal
  range `[-(plen + scope + 2), tlen + scope + 2]` from score 0. C allocates
  over a running range that starts at `[0, 0]` and grows with the wavefronts
  (`lo - (scope + 1)` .. `hi + (scope + 1)`, merged), which still covers the
  cells kernels read up to `max_score_scope` scores later. The computed
  `lo`/`hi` and their clamps are unchanged.
- **Slab mode.** C uses a "tight" slab when non-modular (CIGAR): wavefronts up
  to the initial length (1000) are standard and reused across alignments;
  larger ones get their exact length and are freed after the alignment. Rust
  used "reuse" everywhere, sizing every wavefront at 1.5× the largest seen
  (47,141 wavefronts × 24,880 cells on utest affine2p). Oversize wavefronts
  now come from their own arena, reset after the alignment; standard ones
  stay on the free list.
- **Slab safety.** New wavefronts reuse deallocated slots, so the slab Vec
  can't outgrow its reserved capacity (a reallocation would invalidate the
  pointers components hold); `debug_assert`ed, as are cell accesses outside
  a wavefront's allocation. The arena's slow path now looks for an existing
  chunk that fits instead of appending (it had reached 80,763 chunks).

Notes:
- First version: CIGAR modes on 1 kb were 4–6% slower. Resetting the whole
  slab after any alignment with an oversize wavefront (on 1 kb affine2p,
  the cheap second gap makes some wavefronts wider than 1000) rebuilt every
  wavefront; the separate oversize arena fixed it (affine2p CIGAR +5.5% →
  +2.0%). A per-alignment scan for non-standard lengths cost ~1%; a flag
  replaced it.
- The remaining ±2–3% on short pairs isn't in any one function (profiles
  match main); padding the allocation range shifted affine2p by ~1%, so it's
  likely data placement. Not pursued.
- Parity: differential harness with debug assertions on, 1242/1262
  identical or same-error, only the 20 known HPCBio/wfa2lib-rs#22 cases differ.


## 2026-10-14: `Compat::Corrected` ends-free (HPCBio/wfa2lib-rs#32, #37, #38)

Branch `fix/corrected-endsfree-match-reward` (with HPCBio/wfa2lib-rs#39) vs `main` `5cdcea3`, Apple
Silicon (scalar kernels), `long.seq` (300 × 2 kb), best of 3 runs, CIGAR /
score-only, ms. "C" is the default mode, "corr." is `--no-c-compat`.

| Config | main | C | corr. | Max RSS |
|---|---|---|---|---|
| affine `-5,4,0,8`, ends-free 0.1 ×4 | 1042 / 704 | 1020 / 701 | 1391 / 1019 | 82 / 3 MB |
| affine `-5,4,0,8`, ends-free 1 ×4 (whole length) | 1649 / 1214 | 1652 / 1225 | 2301 / 1795 | 138 / 4 MB |
| affine `-5,4,0,8`, extension | 140 / 119 | 135 / 121 | 182 / 167 | 39 / 2 MB |
| affine `0,4,6,2`, ends-free 0.1 ×4 | 239 / 223 | 239 / 225 | 236 / 231 | 17 / 2 MB |
| affine2p `-1,6,2,4,12,1`, ends-free 0.1 ×4 | 1395 / 989 | 1392 / 970 | 1775 / 1323 | 123 / 5 MB |
| affine `-5,4,0,8`, global | 118 / 105 | 119 / 107 | 116 / 106 | 39 / 2 MB |

- C mode: no change from `main` (within noise). The affine2p slot-leak fix
  (FIXED_BUGS 13) and the med/low seed fix (FIXED_BUGS 12, both in
  HPCBio/wfa2lib-rs#39) don't show.
- Corrected mode costs +27–47% on ends-free/extension with match < 0, and
  nothing elsewhere (match == 0 and global take the C path). It runs on
  past the first boundary hit until the score reaches the best candidate's
  total, `s + M * left`. Memory is unchanged.
- Possible later optimization (exact, not a heuristic): from cell (v, h)
  at score s, any completion totals at least `s + M * |(|P| - v) - (|T| - h)|`.
  Matches cost 0, and every base of the length difference is either a gap
  (> M) or a free-end base (M). Once a candidate exists, diagonals whose
  bound reaches its total can be trimmed, narrowing the wavefronts.
