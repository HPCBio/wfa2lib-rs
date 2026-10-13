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

/// wfa.utest.seq pairs 101 and 103; expected outputs are from C's driver
/// (3282f83) with `-a gap-affine-wfa`.
const DROP_PAIRS: [(&str, &str); 2] = [
    (
        "CACGATCAGGAAGCTGCATCCTTATCTAGGGTGTGTGATCGACAGCTTGGCTGCTCGTCCCACATATGAGTTGCAAGTAACGAATACGTTTCTACGAACC",
        "CTAGAATCAGAAGCATGCCTCCTATCCAGGGTGTGTGATCACAGCTCGCTGCTCGTCCCCATAATGAGTTCACAGTACGAATAAATTTCTACAAACC",
    ),
    (
        "TCGCGTCCAACTGAGTCAAGAGATAATGTCGGTAACTCTTCATACCAGCACGTCCAGCGGTAAGCCAAGAAACTCGAGCGTGAGAGGGTCCTCCGGTTTG",
        "TCGAGTCGCAACTGAGTCAAGAGATTAAGTCGGTAAATCCCTCATACCAGCAGTCTAGCGGTAAGCCAAGAAGCTTGAGCGGGAGAGGGTCCTCCGGTTTG",
    ),
];

fn drop(name: &str, heuristic: &str, params: &str, extra: &[&str]) -> Vec<String> {
    let h = format!("--wfa-heuristic={heuristic}");
    let p = format!("--wfa-heuristic-parameters={params}");
    let mut args = vec!["-a", "gap-affine-wfa", h.as_str(), p.as_str()];
    args.extend_from_slice(extra);
    run(name, &DROP_PAIRS, &args)
}

#[test]
fn xdrop_zdrop_match_c() {
    // X-drop previously had no effect (output equalled the unbanded alignment).
    assert_eq!(drop("xd", "xdrop", "20,1", &[]), ["-1\t-", "-1\t-"]);
    // Z-drop ends at the best-scoring cell: a max-trimmed partial alignment.
    assert_eq!(
        drop("zd", "zdrop", "20,1", &[]),
        ["0\t1M", "-20\t3M1X3M1I17M1X1M1X8M"]
    );
    assert_eq!(
        drop("zd_ext", "zdrop", "20,1", &["--wfa-span=extension"]),
        ["0\t1M", "-20\t3M1X3M1I17M1X1M1X8M"]
    );
    // Score-only with a match reward: C evaluates the classic score at the
    // end position (the Z-drop cell, or the unset position when unreachable).
    let m = ["--affine-penalties=-3,5,3,2", "--wfa-score-only"];
    assert_eq!(drop("xd_s", "xdrop", "20,1", &m), ["163\t-", "213\t-"]);
    assert_eq!(drop("zd_s", "zdrop", "20,1", &m), ["15\t-", "57\t-"]);
}

#[test]
fn drops_reject_what_c_rejects() {
    // C rejects X-drop/Z-drop with edit/indel, and requires P1,P2.
    let edit = run_status(
        "xd_edit",
        &DROP_PAIRS,
        &[
            "-a",
            "edit-wfa",
            "--wfa-heuristic=xdrop",
            "--wfa-heuristic-parameters=20,1",
        ],
    );
    assert!(edit.is_none());
    let one_param = run_status(
        "xd_p1",
        &DROP_PAIRS,
        &[
            "-a",
            "gap-affine-wfa",
            "--wfa-heuristic=xdrop",
            "--wfa-heuristic-parameters=20",
        ],
    );
    assert!(one_param.is_none());
}

#[test]
fn streams_to_stdout_and_keeps_results_written_before_an_error() {
    // Without -o, results go to stdout. Output is streamed as in C's driver,
    // so results written before an error exit (here: the second pair's
    // ends-free parameters exceed its length) are kept, as with C's exit().
    let input = std::env::temp_dir().join(format!("wfa2_cli_stream_{}.seq", std::process::id()));
    fs::write(&input, ">ACGTACGTAC\n<ACGTACGTACGGAAA\n>ACGT\n<ACG\n").unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_align_benchmark"))
        .arg("-i")
        .arg(&input)
        .args(["-a", "gap-affine-wfa", "--wfa-span=ends-free,0,0,0,5"])
        .output()
        .expect("failed to run align_benchmark");
    let _ = fs::remove_file(&input);
    assert!(!out.status.success());
    assert_eq!(String::from_utf8(out.stdout).unwrap(), "-16\t10M5I\n");
}

/// 300 bp pairs (DROP_PAIRS concatenated), long enough for BiWFA to split at
/// breakpoints instead of falling back to unidirectional WFA.
fn biwfa_pairs() -> [(String, String); 2] {
    let [(p1, t1), (p2, t2)] = DROP_PAIRS;
    [
        ([p1, p2, p1].concat(), [t1, t2, t1].concat()),
        ([p2, p1, p2].concat(), [t2, t1, t2].concat()),
    ]
}

fn biwfa(name: &str, heuristic: &str, params: &str, extra: &[&str]) -> Vec<String> {
    let pairs = biwfa_pairs();
    let pairs: Vec<(&str, &str)> = pairs
        .iter()
        .map(|(p, t)| (p.as_str(), t.as_str()))
        .collect();
    let h = format!("--wfa-heuristic={heuristic}");
    let p = format!("--wfa-heuristic-parameters={params}");
    let mut args = vec!["-a", "gap-affine-wfa", "--wfa-memory=ultralow", &h, &p];
    args.extend_from_slice(extra);
    run(name, &pairs, &args)
}

#[test]
fn biwfa_heuristics_match_c() {
    // Expected outputs are from C's driver (3282f83). BiWFA previously ignored
    // heuristics (pair 1 here aligned fully, -292).
    let full = "3M1X3M1I17M1X1M1X8M1X2M1I1X11M1D3M1X16M1X2M1X5M1X20M2X2M1I4M1X1M3X3M1X4M1D3M1X13M1D6M1X1M1D11M1D4M1I6M3X4M1D6M2X7M1X7M1X3M1I17M1X1M1X8M1X2M1I1X11M1D3M1X16M1X2M1X5M1X19M";
    let banded = ["-1\t-".to_string(), format!("-236\t{full}")];
    assert_eq!(biwfa("bi_bs", "banded-static", "-2,2", &[]), banded);
    assert_eq!(biwfa("bi_ba", "banded-adaptive", "-2,2,1", &[]), banded);
    assert_eq!(biwfa("bi_xd", "xdrop", "20,1", &[]), ["-1\t-", "-1\t-"]);
    // Score-only BiWFA only searches for the breakpoint; on failure C prints
    // INT32_MIN.
    assert_eq!(
        biwfa("bi_xd_s", "xdrop", "20,1", &["--wfa-score-only"]),
        ["-2147483648\t-", "-2147483648\t-"]
    );
    let m = ["--affine-penalties=-3,5,3,2", "--wfa-score-only"];
    assert_eq!(
        biwfa("bi_xd_ms", "xdrop", "100,1", &m),
        ["539\t-", "589\t-"]
    );
}

#[test]
fn biwfa_match_reward_tie_breaks_match_c() {
    // wfa.utest.seq pair 298. Sub-problems of at most 100 bp used to fall back
    // to unidirectional WFA at every level; C checks the length only at the
    // top level, which can pick a different co-optimal CIGAR (here 3D1X
    // instead of 1X3D). Expected output from C's driver (3282f83).
    let lines = run(
        "bi_tie",
        &[(
            "AATGCCTCGGCTTCAGTAAGCAATGCCGGCCGTTTAACTGGTATGAAGAAATGATAGCTCGAGACTAAAGAGGGGGAACCTGGTTTTAAGGGCTAAACACTTGGG",
            "AAGCTTTCGTTCGGTGAAACTTGTTGCGTGGTTACGACAGAAATGACGAGATCGTCCATAGCAGGCAGGACTAGGTGACGGCGTACACCAGACGGTAGGT",
        )],
        &[
            "-a",
            "gap-affine-wfa",
            "--affine-penalties=-3,5,3,2",
            "--wfa-memory=ultralow",
        ],
    );
    assert_eq!(
        lines,
        [
            "-3\t2M1D2M1X1M1I2M2D3M1X2M2D1M1X2M2I2M2X1M2D3M6D3M1I1M1X2M1I8M5D5M1I1M1I1M2I1M1X2M1I3M1I1X3M2D2M1I3M1X3D1M1X3M1I2M1X2M3I2M2I1M1X2M1X"
        ]
    );
}

/// Run `align_benchmark` on a file in `tests/data/`; return the output.
fn run_file(name: &str, args: &[&str]) -> String {
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let out = Command::new(env!("CARGO_BIN_EXE_align_benchmark"))
        .arg("-i")
        .arg(data.join(name))
        .args(args)
        .output()
        .expect("failed to run align_benchmark");
    assert!(out.status.success());
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn memory_med_low_match_c() {
    // wfa.utest.seq pair 201 with X-drop every 10 steps. With `med`/`low`
    // (piggyback backtrace) C gives a valid alignment, different from its
    // `high` result (whose CIGAR is invalid, HPCBio/wfa2lib-rs#21). Expected
    // outputs are from C's driver (3282f83). `--wfa-memory=med|low` used to
    // behave like `high`.
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let expected = |f: &str| fs::read_to_string(data.join(f)).unwrap();
    for (mem, alg) in [
        ("med", "xdrop_med.med_low.alg"),
        ("low", "xdrop_med.med_low.alg"),
        ("high", "xdrop_med.high.alg"),
    ] {
        let mem = format!("--wfa-memory={mem}");
        let out = run_file(
            "xdrop_med.seq",
            &[
                "-a",
                "gap-affine-wfa",
                &mem,
                "--wfa-heuristic=xdrop",
                "--wfa-heuristic-parameters=50,10",
            ],
        );
        assert_eq!(out, expected(alg), "{mem}");
    }
}
