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
banded-adaptive, Marco-Sola-Lab/WFA2-lib#125's band clamp with 1000 × 1000 initial
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


## 2026-10-14: corrected-mode ends-free speed (HPCBio/wfa2lib-rs#31, #32)

Branch `perf/corrected-endsfree-pruning` vs `main` `e1ced46`, Apple Silicon,
`long.seq` (300 × 2 kb). The machine was loaded (an endpoint-backup service)
for part of the work, so instructions retired (`/usr/bin/time -l`, load-
independent) are the main measure; wall clock is best of 5, interleaved, at
load ~3. CIGAR / score-only. "C" is the default mode, "corr." `--no-c-compat`.

| Config | main C | main corr. | branch C | branch corr. |
|---|---|---|---|---|
| affine `-5,4,0,8`, ends-free 0.1 ×4 (G instr.) | 14.54 / 14.25 | 18.30 / 18.01 | 14.54 / 14.25 | 14.87 / 14.58 |
| same (ms) | 1092 / 715 | 1451 / 1054 | 1103 / 717 | 1181 / 796 |
| affine `-5,4,0,8`, ends-free 1 ×4 (G instr.) | 25.42 / 24.94 | 32.35 / 31.87 | 25.42 / 24.94 | 26.02 / 25.54 |
| same (ms) | 1713 / 1230 | 2344 / 1856 | 1710 / 1237 | 1854 / 1382 |
| affine `-5,4,0,8`, extension (G instr.) | 2.54 / 2.45 | 3.10 / 3.00 | 2.54 / 2.45 | 2.61 / 2.52 |
| same (ms) | 145 / 119 | 198 / 169 | 145 / 120 | 159 / 130 |
| affine2p `-1,6,2,4,12,1`, ends-free 0.1 ×4 (G instr.) | 17.79 / 17.33 | 21.78 / 21.33 | 17.81 / 17.34 | 18.15 / 17.70 |
| same (ms) | 1444 / 1008 | 1814 / 1372 | 1451 / 1008 | 1533 / 1087 |
| affine ends-free 0.1 ×4, banded-static ±200 (G instr.) | 11.63 / 11.33 | 14.46 / 14.17 | 11.62 / 11.32 | 11.88 / 11.59 |
| same (ms) | 904 / 560 | 1160 / 816 | 880 / 561 | 956 / 613 |
| affine `-5,4,0,8`, global (ms) | 128 / 105 | 127 / 106 | 128 / 105 | 129 / 105 |

- Corrected-mode overhead over C mode: +26–28% instructions (+26–51% time)
  → +2–2.5% instructions (+7–12% time). C mode unchanged.
- Cause: not the extra steps. Corrected mode ran only 0.8% more steps than
  C on `long.seq`, but its extend took 74% more samples (samply). The
  `_best` extend kernel looped `for k in lo..=hi` and called
  `termination::endsfree_left` (two `Option`s, merged) on every diagonal.
  Now it walks the offsets as the C-mode kernel does and checks
  `termination_endsfree` first, counting the bases left only at an end.
- Pruning (FIXED_BUGS 14a): exact, but worth only 0.3–2% here, since the
  first candidate is usually near-optimal on these pairs. On pairs with
  random flanks and whole-length free ends (300 × 1–2 kb + up to 300 bp
  flanks), corrected mode is still +15–55% instructions over C mode
  (+160% on 500 short pairs with `-1,4,6,2`). There C stops at a wrong,
  #38-style first hit, so the extra search is what makes the result
  optimal; pruning trims 0.3–2% of it.
- Output: pruned vs unpruned byte-identical on 360 configs × 506 generated
  pairs and 240 configs × `long.seq` (linear/affine/affine2p, five spans,
  high/med/low, CIGAR and score-only, no heuristic and two banded-static
  bands), once both include the null-step fix (FIXED_BUGS 14c).


## 2026-10-15: vectorized, fused piggyback kernels (med/low, HPCBio/wfa2lib-rs#31)

Branch `perf/piggyback-fused` vs `main` `99943f4`, Apple M2 Pro, best of 3
(load under 3), CIGAR, ms. Data from the #31 sweep: `m1k_5` (20k × 1 kb,
5%), `m1k_20` (5k × 1 kb, 20%). C is `3282f83` (scalar on ARM, so C and
C `-march=native` time the same here).

| Case | C | main | fused | fused/main | fused/C | instr. fused/main |
|---|---|---|---|---|---|---|
| m1k_5 edit med | 286 | 259 | 185 | 0.71 | 0.65 | 0.76 |
| m1k_5 edit low | 282 | 211 | 175 | 0.83 | 0.62 | 0.76 |
| m1k_5 indel med | 200 | 229 | 167 | 0.73 | 0.83 | 0.77 |
| m1k_5 indel low | 202 | 226 | 168 | 0.74 | 0.83 | 0.77 |
| m1k_5 linear med | 392 | 351 | 237 | 0.67 | 0.60 | 0.72 |
| m1k_5 linear low | 435 | 407 | 287 | 0.70 | 0.66 | 0.74 |
| m1k_5 affine med | 4160 | 3325 | 2241 | 0.67 | 0.54 | 0.57 |
| m1k_5 affine low | 4610 | 4032 | 2958 | 0.73 | 0.64 | 0.63 |
| m1k_5 affine2p med | 27180 | 20665 | 11721 | 0.57 | 0.43 | 0.49 |
| m1k_5 affine2p low | 29470 | 24332 | 15345 | 0.63 | 0.52 | 0.57 |
| m1k_20 edit med | 761 | 466 | 243 | 0.52 | 0.32 | 0.58 |
| m1k_20 edit low | 749 | 464 | 240 | 0.52 | 0.32 | 0.58 |
| m1k_20 indel med | 421 | 502 | 311 | 0.62 | 0.74 | 0.61 |
| m1k_20 indel low | 386 | 498 | 319 | 0.64 | 0.83 | 0.61 |
| m1k_20 linear med | 1090 | 891 | 478 | 0.54 | 0.44 | 0.56 |
| m1k_20 linear low | 1200 | 1079 | 647 | 0.60 | 0.54 | 0.61 |
| m1k_20 affine med | 9800 | 7584 | 4836 | 0.64 | 0.49 | 0.51 |
| m1k_20 affine low | 10880 | 9252 | 6474 | 0.70 | 0.60 | 0.59 |
| m1k_20 affine2p med | 73800 | 53301 | 30158 | 0.57 | 0.41 | 0.46 |
| m1k_20 affine2p low | 80400 | 64311 | 40923 | 0.64 | 0.51 | 0.54 |

- Cause: none of the piggyback kernels vectorized in the linked (LTO)
  binary, including the per-component passes split out for that purpose.
  LLVM turned `if c { a[i] } else { b[i] }` (a choice between two loads) into
  a load from a selected address, a gather, which stops vectorization.
  The crate's own `--emit asm` is misleading with `lto = true`: it is the
  pre-link pipeline, which skips loop vectorization. Check the linked
  binary (`objdump -d`).
- Fix: choose between loaded values with a mask blend
  (`(a & m) | (b & !m)`); fuse affine/affine2p into one loop per kernel,
  as C (M from this cell's I/D in registers); each loop in an
  `#[inline(never)]` function with every array a separate slice parameter
  (noalias). All five kernels now vectorize (NEON here).
- med/low is now 0.52–0.83× `main` and 0.32–0.83× C (before: indel med/low
  1.12–1.34× C). Output identical to C and `main` in every case; the
  differential harness is clean (1242/1262, the 20 known #22 cases).
- On x86 (EPYC 7713, quick sweep before this change) med/low was 1.3–1.9×
  slower than C built with `-march=native` (GCC vectorizes C's fused loops
  with AVX2). To re-check with `bench/sweep.py`.


## 2026-10-15: x86 extend stalls on AVX2 gathers (HPCBio/wfa2lib-rs#31)

`bench/sweep.py` on `hpcbio` (Xeon Gold 6150, Skylake-SP; `main` `aa091cf`)
had Rust up to 2.16× slower than generic (scalar) C on long or divergent
edit/indel/linear alignments (`l10k_5`, `m1k_20`), while executing about 0.6×
the instructions. The EPYC 7713 quick run showed the same pattern (indel
1.22–1.25×). `perf` on `l10k_5 indel score`:

| | C (scalar) | Rust |
|---|---|---|
| cycles | 0.70 G | 1.50 G |
| instructions | 2.74 G | 1.56 G |
| stall cycles (total) | 0.013 G | 0.487 G |
| stalls on L1 miss | 0.5 M | 4.6 M |
| 4K-alias events | 2.5 M | 5.3 M |

- 88% of Rust's time was in `extend_end2end_avx2`, 65% of that on its two
  `vpgatherdd` (gather 4 bytes of pattern and text for 8 diagonals). The
  stalls are gather latency, not memory (L1/L2/LLC misses, store buffer and
  4K aliasing are all small).
- Rust used the kernel only for wavefronts of 256+ diagonals, so the 1 kb/5%
  benchmark it was tuned on rarely ran it. C's AVX2 build uses the same
  gathers (slower than its scalar build on the EPYC too); on Skylake-SP
  C-native uses its AVX-512 kernel instead.
- Fix: removed the gather kernels (end-to-end and the BiWFA max variant);
  x86 now uses the blockwise 64-bit scalar extend, as on ARM. Output is
  unchanged (same offsets). To confirm with `bench/sweep.py` on x86.


- Confirmed on the EPYC 7713 (`bench/sweep.py`, `8ccfb64`, 3 runs per case,
  vs the pre-fix quick run): `l10k_5` indel score 1.69× → 0.80× C,
  `m1k_20` indel CIGAR 1.25× → 0.87×; edit unchanged (it rarely hit the
  kernel). With #45, 1k/10k med/low went from a median of 0.94× to 0.45×
  generic C.

## 2026-10-15: BiWFA overlap scan vectorized (HPCBio/wfa2lib-rs#31)

The full EPYC sweep after #45/#48 left affine2p BiWFA (ultralow) as the main
case slower than C: 1.08–1.15× generic C, ~2.2× C-native, at 1.2–1.4× C's
instructions. A `sample` profile on the M2 (`m1k_5` affine2p ultralow) put
73% of the time in `bialign_overlap`, against 11% each for compute and
extend.

- At every step BiWFA checks the new wavefronts against the last
  `max_score_scope` scores of the other direction: one M scan plus four
  I1/D1/I2/D2 scans for affine2p (scope ~25 with the default penalties).
  Rust had NEON only for the M scan; the indel scans and every scan on x86
  were scalar. C has AVX2/AVX-512 versions of both, used when built with
  `-march=native`: the C-native gap.
- Fix: one `scan_overlap` for all five scans. Chunks of 8 diagonals are
  tested at once (NEON: 2 × 4 lanes; AVX2: one `vpermd`-reversed load,
  inlined with `target-cpu=native`, runtime-detected otherwise); only a chunk
  with a candidate is rechecked one diagonal at a time. The visiting order is
  unchanged, so the breakpoint is the same as the scalar scan (and C's).
- M2 (`data-quick`, Rust vs `main` vs C 3282f83, ultralow, CIGAR output):

| Case | `main` | this | C | this / C |
|---|---|---|---|---|
| m1k_5 affine2p | 1612 ms | 732 ms | 1.58 s | 0.46 |
| m1k_5 affine2p m-1 | 7118 ms | 2799 ms | 7.45 s | 0.38 |
| l10k_5 affine2p | 1325 ms | 702 ms | 1.36 s | 0.52 |
| l10k_5 affine2p m-1 | 3576 ms | 1474 ms | 3.63 s | 0.41 |
| m1k_5 affine | 132 ms | 122 ms | 149 ms | 0.82 |
| m1k_5 edit / indel / linear | | unchanged | | 0.83 / 0.85 / 0.81 |

- Output identical to C and `main` in all cases above, and in x86 builds
  (AVX2 and generic, run under Rosetta) on subsets. Differential harness:
  1242/1262 identical, the 20 known #22 cases. New unit tests compare the
  scan with a naive scan (tails, nulls, rejected indel hits) and the SIMD
  chunk tests with scalar; both catch a deliberately broken lane pairing
  and a scan that skips the rest of a chunk.
- To confirm on x86 with `bench/sweep.py` (expect affine2p ultralow near or
  below C-native).

## 2026-10-17: x86 prefetches: claim output lines instead (HPCBio/wfa2lib-rs#31)

After #49, the Xeon Gold 6150 (Skylake-SP) was the only machine where
gap-affine 2-piece high-memory CIGAR wasn't faster than C (1.01×; banded 1.28×,
against 0.75–0.77× on the EPYC 7713). Rust's CIGAR-minus-score overhead there
was +32% against +8% for C, and +8–10% for Rust on the M2 and EPYC.

- `perf` (`l2_rqsts.rfo_miss`, `resource_stalls.sb`, sampled per instruction):
  the store buffer was full for 4.67 G cycles against C's 0.38 G, about 60%
  of them on the AVX2 affine2p kernel's five output stores. Those stores
  missed L2 for ownership (RFO) 158 M times against about 88 M in C's kernel.
  In high memory the outputs are fresh lines.
- Ruled out by A/B or counters: per-component victim outputs (no change),
  split stores (12% vs C's 10%), slab cache coloring (helps only banded on
  the Xeon, costs up to 4% on the EPYC), and the input prefetches (not the
  cause, but they didn't help either).
- Fix: drop all x86 input prefetches (the hardware prefetcher covers these
  sequential streams) and issue `PREFETCHW` on the five affine2p outputs,
  32 diagonals (two lines) ahead.

A/B at `45f419e`, best of 3, time relative to the unchanged build, with
output identical in all cases:

| Case | Xeon: 32 | Xeon: 64 | Xeon: 128 | EPYC: 32 |
|---|---|---|---|---|
| affine2p CIGAR (m1k_5) | **0.88** | 0.98 | 1.08 | 0.96 |
| affine2p banded CIGAR (m1k_20) | **0.65** | 0.91 | 1.05 | 1.03 |
| affine2p CIGAR (l10k_5) | **0.89** | 0.92 | 0.96 | 1.00 |
| affine2p score | **0.93** | 0.95 | 0.97 | 0.94 |
| affine2p BiWFA | **0.95** | 0.98 | 1.03 | 0.96 |
| gap-affine CIGAR / score | 0.98 / 0.98 | | | 1.01 / 0.97 |
| edit (control) | 0.99 | | | 0.99 |

- On the Xeon, store-buffer stalls on affine2p CIGAR dropped from 4.63 G to
  0.43 G cycles, about C's level. Against C: affine2p CIGAR 1.01× → 0.89×,
  banded 1.28× → 0.83×, 10 kb 0.97× → 0.87×.
- On the EPYC it is about the same as just removing the prefetches (2–6%
  faster on most affine2p cases, the rest within noise). ARM is unchanged.
