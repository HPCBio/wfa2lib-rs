//! Integration tests for the `align_benchmark` output format, which must match
//! the C driver's (`benchmark_print_output`) so outputs can be diffed directly.
#![cfg(feature = "cli")]

use std::fs;
use std::path::PathBuf;
use std::process::Command;

/// Run `align_benchmark` on `pairs` with `args`; return the output lines.
fn run(name: &str, pairs: &[(&str, &str)], args: &[&str]) -> Vec<String> {
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
    assert!(status.success(), "align_benchmark exited with {status}");

    let lines = fs::read_to_string(&output)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
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
