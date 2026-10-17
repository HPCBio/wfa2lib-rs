# wfa2lib-rs

**A pure Rust port of [WFA2-lib](https://github.com/Marco-Sola-Lab/WFA2-lib)**, the wavefront alignment library for pairwise sequence alignment.

This project was created with [Claude Code](https://claude.ai/claude-code).

## Overview

wfa2lib-rs is a from-scratch Rust implementation of the Wavefront Alignment (WFA) algorithm. WFA exploits the structure of sequence similarity to achieve O(ns) time complexity — proportional to the sequence length (n) and the alignment score (s) — making it significantly faster than traditional O(n²) dynamic programming approaches for similar sequences.

This implementation covers the full WFA feature set:

- **Distance metrics**: indel (LCS), edit (Levenshtein), gap-linear (Needleman-Wunsch), gap-affine (Smith-Waterman-Gotoh), and gap-affine 2-piece (concave penalties)
- **Alignment modes**: end-to-end (global), ends-free (semi-global), and extension alignment
- **BiWFA**: bidirectional WFA for O(s) memory alignment with full CIGAR output
- **Heuristics**: WF-adaptive, X-drop, Z-drop, static banding, and adaptive banding for pruning divergent sequences
- **Sequence representations**: ASCII byte sequences, packed 2-bit DNA encoding, and a lambda callback interface for custom sequence access
- **SIMD acceleration**: NEON-vectorized compute kernels on aarch64 (4 diagonals/iteration), AVX2 on x86_64 (8 diagonals/iteration), with scalar fallback
- **Score and alignment**: score-only mode for speed, or full CIGAR traceback with support for match rewards via Eizenga's formula

## Performance

Rust time / C time (below 1: Rust is faster), median over the cases in each
group, measured with `bench/sweep.py` at commit `fc51088` against C WFA2-lib
`3282f83` built without `-march=native`. Output is identical to C in all 301
cases on every machine.

| Cases | Apple M2 Pro | AMD EPYC 7713 | Intel Xeon Gold 6150 |
|---|---|---|---|
| Score only, high memory | 0.97 | 0.71 | 0.85 |
| CIGAR, high memory | 0.98 | 0.73 | 0.92 |
| CIGAR, med/low memory | 0.52 | 0.44 | 0.51 |
| CIGAR, BiWFA (ultralow memory) | 0.95 | 0.75 | 0.79 |
| Ends-free / extension, high memory | 0.84 | 0.60 | 0.71 |
| Heuristics (adaptive, X/Z-drop, banded) | 0.90 | 0.78 | 0.89 |
| 150 bp reads, med/low memory | 1.24 | 1.20 | 1.34 |
| **All 301 cases** | **0.87** | **0.69** | **0.80** |

Rust gains most with med/low memory, ends-free and extension alignments, and
gap-affine 2-piece BiWFA. Known gaps: 150 bp reads with med/low memory are
slower than C everywhere. On the Xeon, an older design (Skylake-SP, 2017),
both tools are slower than on the EPYC, but Rust more so (median 1.47× its
EPYC time, against 1.29× for C), and gap-affine 2-piece CIGAR is about at
parity with C. High-memory CIGAR on long pairs can use up to 1.9× C's peak
memory. Details, datasets, per-algorithm tables and how to reproduce them are
in [`docs/performance.md`](docs/performance.md); ongoing work is tracked in
HPCBio/wfa2lib-rs#31.

## Building

Requires Rust 1.91+ (edition 2024).

```bash
# Debug build
cargo build

# Optimized release build (recommended)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Run tests (unit, integration, and CLI tests)
cargo test
```

### Chip-specific optimizations

On Skylake-family Intel CPUs (6th–10th generation Core, Skylake-SP and
Cascade Lake Xeon), the microcode fix for the
[JCC erratum](https://www.intel.com/content/www/us/en/content-details/841076/intel-mitigations-for-jump-conditional-code-erratum.html)
makes run time depend on code alignment. On a Xeon Gold 6150, building with
LLVM's mitigation made edit, indel, gap-linear and X-drop/Z-drop alignments
7–17% faster, but gap-affine 4% slower; output is unchanged:

```bash
RUSTFLAGS="-C target-cpu=native -C llvm-args=-x86-branches-within-32B-boundaries" cargo build --release
```

Other CPUs don't need it. For example: on an AMD EPYC 7713 (Zen 3) it was
neutral to 7% slower. `RUSTFLAGS` applies to the whole build, so a crate that
depends on wfa2lib-rs sets it in its own build. To reference specific
measurements, see HPCBio/wfa2lib-rs#31.

### Tests, Parity, Performance, and CI

The integration tests (`tests/test_vs_reference.rs`) and the C-parity script
(`tests/wfa_utest.sh`) compare against the reference outputs shipped with the
C library, so they expect a [WFA2-lib](https://github.com/Marco-Sola-Lab/WFA2-lib)
checkout at `./WFA2-lib` (a symlink works; the path is gitignored):

```bash
git clone https://github.com/Marco-Sola-Lab/WFA2-lib   # or: ln -s /path/to/WFA2-lib WFA2-lib
cargo test
bash tests/wfa_utest.sh   # Rust align_benchmark vs WFA2-lib/tests/wfa.utest.check
```

Without it, `tests/test_vs_reference.rs` fails with "failed to read sequence
file". The CLI tests (`tests/cli_output.rs`, `cli` feature) need no checkout.

For broader parity checks, `tests/differential/run.sh` builds C's
`align_benchmark` at a pinned commit and compares it with Rust's on about
1300 configurations (metrics, spans, heuristics, memory modes); see
[`tests/differential/README.md`](tests/differential/README.md).

For performance, `bench/sweep.py` builds both tools (C also with
`-march=native`), generates fixed datasets and compares time, instructions,
memory and output across metrics, memory modes, spans and heuristics; it
writes `target/bench-sweep/<host>-<commit>/summary.md`. Run it on an idle
machine (`--quick` for a short check). `bench/report.py` turns sweeps
from several machines into the tables in
[`docs/performance.md`](docs/performance.md). Ongoing work is tracked in
HPCBio/wfa2lib-rs#31.

CI (`.github/workflows/ci.yml`) runs clippy, the tests (release, and debug
with overflow checks and debug assertions), a library-only build,
`wfa_utest.sh` and the differential harness on Linux (x86-64) and macOS
(arm64) for every pull request and push to `main`. Running the harness on both
checks that Rust's x86 and ARM SIMD paths give the same output. Linux also
runs the release tests and the harness with
`RUSTFLAGS="-C target-cpu=x86-64-v3"`, so the AVX2 kernels are tested both
when selected at run time and when enabled at compile time (as in
`target-cpu=native` builds).

## Usage

### As a library

```rust
use wfa2lib_rs::aligner::{AlignmentScope, WavefrontAligner};
use wfa2lib_rs::penalties::{AffinePenalties, DistanceMetric, WavefrontPenalties};

// Configure gap-affine penalties (mismatch=4, gap_open=6, gap_extend=2)
let penalties = WavefrontPenalties::new_affine(AffinePenalties {
    match_: 0,
    mismatch: 4,
    gap_opening: 6,
    gap_extension: 2,
});

let mut aligner = WavefrontAligner::new(penalties);
aligner.alignment_scope = AlignmentScope::ComputeAlignment;

let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
let text    = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

let score = aligner.align_end2end(pattern, text);
println!("Score: {}", score);
println!("CIGAR: {}", aligner.cigar);
```

### CLI benchmark tool

The included `align_benchmark` binary is compatible with the WFA2-lib reference interface:

```bash
# Gap-affine alignment
cargo run --release --bin align_benchmark -- \
    -a gap-affine-wfa \
    -g 0,4,6,2 \
    -i sequences.seq

# BiWFA (memory-efficient)
cargo run --release --bin align_benchmark -- \
    -a gap-affine-wfa \
    -g 0,4,6,2 \
    --wfa-memory ultralow \
    -i sequences.seq

# Score-only mode
cargo run --release --bin align_benchmark -- \
    -a gap-affine-wfa \
    -g 0,4,6,2 \
    --wfa-score-only \
    -i sequences.seq
```

`--wfa-memory` selects `high` (default), `med`, `low` or `ultralow`, as in C.
`med` and `low` keep only the wavefronts still needed and record the
alignment alongside them ("piggyback" backtrace); for co-optimal alignments
they can choose differently from `high`, exactly as C does. In the library:
`set_memory_mode(MemoryMode::Med)`.

Input files use the WFA2-lib format: each pair is two lines, `>pattern` then
`<text`. Output is one `score<TAB>CIGAR` line per pair (`-` for no CIGAR), written
as each pair is aligned. Like the C tool, it also supports `--wfa-span`
(`global`, `extension`, `ends-free,P0,Pf,T0,Tf`) and `--wfa-heuristic` with
`--wfa-heuristic-parameters` (`banded-static`, `banded-adaptive`,
`wfa-adaptive`, `xdrop`, `zdrop`).

By default output is identical to C, including results C gets wrong. With
`--no-c-compat` (library: `set_compat(Compat::Corrected)`), known C bugs are
fixed instead. Currently this covers ends-free alignment with a match reward
(`match < 0`), where C can return a suboptimal alignment (corrected mode
returns the optimal one), X-drop with `steps_between_cutoffs > 1`, where C
can return an invalid CIGAR, and Z-drop, where C backtraces the partial
alignment from the wrong score (empty or invalid CIGARs, and results that
differ between memory modes). With `match == 0` (the default penalties), it
also fixes the sign of the score X-drop and Z-drop use, which in C makes
X-drop abandon most divergent pairs. See `Compat` in the API docs.

## Project structure

```
src/
  lib.rs              # Public module exports
  aligner.rs          # Main WavefrontAligner struct and alignment loops
  penalties.rs        # Distance metrics and penalty configurations
  sequences.rs        # Sequence representations (ASCII, packed 2-bit, lambda)
  wavefront.rs        # Wavefront data structure
  wavefront_set.rs    # Grouped wavefront sets (M, I, D components)
  slab.rs             # Wavefront memory pool with arena allocation
  arena.rs            # Byte-oriented bump allocator for wavefront arrays
  compute/            # WFA compute kernels (edit, linear, affine, affine2p;
                      #   piggyback.rs: med/low-memory kernels and offloading)
  extend.rs           # Wavefront extend with SIMD acceleration
  backtrace.rs        # CIGAR traceback from wavefronts
  bialign.rs          # Bidirectional WFA (BiWFA)
  cigar.rs            # CIGAR operations and scoring
  heuristic.rs        # Pruning heuristics
  offset.rs           # Wavefront offset type and coordinate helpers
  pcigar.rs           # Packed CIGAR (piggyback backtrace, med/low memory)
  bt_buffer.rs        # Backtrace buffer and compaction (med/low memory)
  components.rs       # Wavefront component management
  termination.rs      # Alignment termination conditions
  bin/
    align_benchmark.rs  # CLI binary
tests/
  cli_output.rs       # CLI output tests (expected outputs from C)
  test_vs_reference.rs  # Library vs C reference outputs (needs ./WFA2-lib)
  wfa_utest.sh        # Stock C unit tests on the Rust CLI (needs ./WFA2-lib)
  data/               # Fixtures for the CLI tests
  differential/       # Rust-vs-C differential harness (pinned C commit)
bench/                # Rust-vs-C timing scripts (sweep.py, report.py)
docs/
  performance.md      # Rust-vs-C results on three machines
.github/workflows/
  ci.yml              # CI: clippy, tests, wfa_utest.sh, differential harness
```

## License

BSD 3-Clause. See [LICENSE](LICENSE).

## Citations

If you use this software, please cite the original WFA papers:

**WFA — Wavefront Alignment**:
> Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa. "Fast gap-affine pairwise alignment using the wavefront algorithm." *Bioinformatics*, Volume 37, Issue 4, February 2021, Pages 456–463. https://doi.org/10.1093/bioinformatics/btaa777

**BiWFA — Bidirectional WFA**:
> Santiago Marco-Sola, Jordan M Eizenga, Andrea Guarracino, Benedict Paten, Erik Garrison, Miquel Moreto. "Optimal gap-affine alignment in O(s) space." *Bioinformatics*, Volume 39, Issue 2, February 2023, btad074. https://doi.org/10.1093/bioinformatics/btad074
