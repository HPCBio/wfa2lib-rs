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

Benchmarked against the C reference implementation (WFA2-lib), 50,000 sequence pairs of 1000 bp at 5% divergence:

### x86_64 (Intel Broadwell, AVX2)

| Benchmark | Rust | C | Ratio |
|---|---|---|---|
| Edit distance (score) | 405 ms | 446 ms | **0.91x** |
| Gap-linear (score) | 682 ms | 1090 ms | **0.63x** |
| Gap-affine (score) | 3980 ms | 3970 ms | **1.00x** |
| Gap-affine (CIGAR) | 5440 ms | 5380 ms | 1.01x |
| Gap-affine 2-piece (score) | 16800 ms | 19710 ms | **0.85x** |
| Gap-affine 2-piece (CIGAR) | 20620 ms | 22070 ms | **0.93x** |
| BiWFA (ultralow memory) | 7133 ms | 6660 ms | 1.07x |

### aarch64 (Apple M-series, NEON)

| Benchmark | Rust vs C | Notes |
|---|---|---|
| Edit distance (score) | **0.56x** | 1.8× faster |
| Gap-linear (score) | **0.68x** | 1.5× faster |
| Gap-affine (score) | **0.85x** | 1.2× faster |
| Gap-affine (CIGAR) | **0.90x** | 1.1× faster |
| Gap-affine 2-piece (score) | **0.88x** | 1.1× faster |
| Gap-affine 2-piece (CIGAR) | **0.86x** | 1.2× faster |
| BiWFA (ultralow memory) | **0.86x** | 1.2× faster |

Ratio < 1.0 means Rust is faster. Key optimizations include hand-written AVX2/NEON SIMD compute kernels (8 or 4 diagonals/iteration) with software prefetching, pre-centered wavefront offset pointers (matching C's `offsets = offsets_mem - min_lo` pattern), separate per-component wavefront arrays for optimal cache stride, arena bump allocation with cache coloring to eliminate L1 conflict misses, conditional kernel dispatch (affine fallback when O2/E2 inputs are null), wavefront reuse between alignments, and blockwise 64-bit sequence extension.

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

Other CPUs don't need it. For example: on a AMD EPYC 7713 (Zen 3) it was
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
machine (`--quick` for a short check). Results are tracked in
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
returns the optimal one), and edit/indel ends-free with begin gaps, where C
can return an invalid CIGAR. See `Compat` in the API docs.

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
bench/                # Rust-vs-C timing scripts
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
