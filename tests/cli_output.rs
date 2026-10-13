//! Integration tests for the `align_benchmark` output format, which must match
//! the C driver's (`benchmark_print_output`) so outputs can be diffed directly.
#![cfg(feature = "cli")]

use std::fs;
use std::path::PathBuf;
use std::process::Command;

/// Run `align_benchmark` on `pairs` with `args`; return the output lines.
fn run(name: &str, pairs: &[(&str, &str)], args: &[&str]) -> Vec<String> {
    run_status(name, pairs, args).expect("align_benchmark failed")
}

/// Like `run`, but return `None` if `align_benchmark` exits unsuccessfully.
fn run_status(name: &str, pairs: &[(&str, &str)], args: &[&str]) -> Option<Vec<String>> {
    let dir = std::env::temp_dir();
    let input: PathBuf = dir.join(format!("wfa2_cli_{name}_{}.seq", std::process::id()));
    let output: PathBuf = dir.join(format!("wfa2_cli_{name}_{}.alg", std::process::id()));
    let seq: String = pairs.iter().map(|(p, t)| format!(">{p}\n<{t}\n")).collect();
    fs::write(&input, seq).unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_align_benchmark"))
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(args)
        .output()
        .expect("failed to run align_benchmark")
        .status;

    let lines = status.success().then(|| {
        fs::read_to_string(&output)
            .unwrap()
            .lines()
            .map(str::to_string)
            .collect()
    });
    let _ = fs::remove_file(&input);
    let _ = fs::remove_file(&output);
    lines
}

#[test]
fn empty_cigar_prints_minus_one_and_dash() {
    // End diagonal (k = -30) outside the [-5, 5] band: no alignment is found.
    // C prints "-1\t-" for an empty CIGAR; previously this printed "0\t".
    let pattern = "ACGTTGCAAGCT".repeat(10);
    let text = &pattern[..90];
    for algo in [
        "edit-wfa",
        "gap-linear-wfa",
        "gap-affine-wfa",
        "gap-affine2p-wfa",
    ] {
        let lines = run(
            "unreachable",
            &[(&pattern, text)],
            &[
                "-a",
                algo,
                "--wfa-heuristic=banded-static",
                "--wfa-heuristic-parameters=-5,5",
            ],
        );
        assert_eq!(lines, ["-1\t-"], "{algo}");
    }
}

#[test]
fn completed_alignment_prints_score_and_cigar() {
    let lines = run(
        "completed",
        &[("ACGTACGT", "ACGTACGT"), ("ACGTACGT", "ACGAACGT")],
        &["-a", "gap-affine-wfa"],
    );
    assert_eq!(lines, ["0\t8M", "-4\t3M1X4M"]);
}

#[test]
fn score_only_prints_dash_cigar() {
    let lines = run(
        "score_only",
        &[("ACGTACGT", "ACGAACGT")],
        &["-a", "gap-affine-wfa", "--wfa-score-only"],
    );
    assert_eq!(lines, ["-4\t-"]);
}

/// Pairs and `--wfa-span` outputs below are from C WFA2-lib's `align_benchmark`
/// (16969b5b and 3282f83 agree) with `-a gap-affine-wfa`.
const SPAN_PAIRS: [(&str, &str); 2] = [
    ("ACGTACGTACGTTTTTTTTT", "ACGTACGTACGTGGGGGGGG"),
    ("ACGTACGTAC", "ACGTACGTACGGAAA"),
];

fn span(name: &str, span: &str, score_only: bool) -> Vec<String> {
    let span_arg = format!("--wfa-span={span}");
    let mut args = vec!["-a", "gap-affine-wfa", span_arg.as_str()];
    if score_only {
        args.push("--wfa-score-only");
    }
    run(name, &SPAN_PAIRS, &args)
}

#[test]
fn wfa_span_matches_c() {
    // Previously `--wfa-span` was parsed but ignored (always global).
    assert_eq!(span("ext", "extension", false), ["0\t12M", "0\t10M"]);
    assert_eq!(span("ext_s", "extension", true), ["-22\t-", "0\t-"]);
    // Text end free: absolute count (5) and fraction of |T| (0.5).
    assert_eq!(
        span("ef_abs", "ends-free,0,0,0,5", false),
        ["-44\t12M5D3X5I", "-16\t10M5I"]
    );
    assert_eq!(
        span("ef_abs_s", "ends-free,0,0,0,5", true),
        ["-28\t-", "0\t-"]
    );
    assert_eq!(
        span("ef_frac", "ends-free,0,0,0,0.5", false),
        ["-44\t12M8D8I", "-16\t10M5I"]
    );
    assert_eq!(
        span("ef_frac_s", "ends-free,0,0,0,0.5", true),
        ["-22\t-", "0\t-"]
    );
    assert_eq!(
        span("global", "global", false),
        ["-32\t12M8X", "-16\t10M5I"]
    );
}

#[test]
fn wfa_span_rejects_what_c_rejects() {
    // BiWFA does not support ends-free/extension; free ends larger than the
    // sequences are an error. C's driver exits with status 1 for both.
    let biwfa = run_status(
        "biwfa_ext",
        &SPAN_PAIRS,
        &[
            "-a",
            "gap-affine-wfa",
            "--wfa-memory=ultralow",
            "--wfa-span=extension",
        ],
    );
    assert!(biwfa.is_none());
    let oversize = run_status(
        "oversize",
        &SPAN_PAIRS,
        &["-a", "gap-affine-wfa", "--wfa-span=ends-free,0,0,0,30"],
    );
    assert!(oversize.is_none());
}

#[test]
fn ends_free_begin_gaps_match_c() {
    // Begin gaps (P0 = 5, T0 = 3) were ignored, giving the global result.
    // Expected outputs are from C's driver (16969b5b and 3282f83 agree).
    let pairs = [
        ("TTTTTACGTACGTACGT", "ACGTACGTACGT"),
        ("ACGTACGTACGT", "GGGACGTACGTACGT"),
        ("ACGTTCGTACGA", "CCACGTACGTACGA"),
    ];
    let span = "--wfa-span=ends-free,5,0,3,0";
    // match == 0: begin-free cells seeded at score 0.
    let lines = run("bf_m0", &pairs, &["-a", "gap-affine-wfa", span]);
    assert_eq!(lines, ["-16\t5D12M", "-12\t3I12M", "-14\t2I4M1X7M"]);
    let lines = run(
        "bf_m0_s",
        &pairs,
        &["-a", "gap-affine-wfa", "--wfa-score-only", span],
    );
    assert_eq!(lines, ["0\t-", "0\t-", "-4\t-"]);
    // match < 0: begin-free cells added as the score reaches them
    // (global score-only would be 23, 27, 21).
    let m = "--affine-penalties=-3,5,3,2";
    let lines = run(
        "bf_mneg_s",
        &pairs,
        &["-a", "gap-affine-wfa", m, "--wfa-score-only", span],
    );
    assert_eq!(lines, ["36\t-", "36\t-", "28\t-"]);
}
