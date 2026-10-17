//! `Compat::Corrected` X-drop with `steps_between_cutoffs > 1`
//! (HPCBio/wfa2lib-rs#21): C's high-memory backtrace reads wavefronts the
//! cutoff trimmed to empty, giving invalid CIGARs; corrected mode is valid,
//! with the same score.

use wfa2lib_rs::aligner::{AlignmentScope, Compat, WavefrontAligner};
use wfa2lib_rs::cigar::Cigar;
use wfa2lib_rs::heuristic::HeuristicStrategy;
use wfa2lib_rs::penalties::{
    Affine2pPenalties, AffinePenalties, LinearPenalties, WavefrontPenalties,
};

/// `wfa.utest.seq` pair 201 (C's default penalties for each metric).
fn pair() -> (Vec<u8>, Vec<u8>) {
    let data = std::fs::read_to_string("tests/data/issue21_xdrop.seq").unwrap();
    let mut lines = data.lines();
    let pat = lines.next().unwrap().trim_start_matches('>').as_bytes().to_vec();
    let txt = lines.next().unwrap().trim_start_matches('<').as_bytes().to_vec();
    (pat, txt)
}

/// Align with X-drop in the given mode; returns the score and the CIGAR.
fn run<const N: usize>(
    penalties: WavefrontPenalties,
    xdrop: (i32, i32),
    compat: Compat,
    pat: &[u8],
    txt: &[u8],
) -> (i32, Cigar) {
    let mut a = WavefrontAligner::<N>::new(penalties);
    a.set_alignment_scope(AlignmentScope::ComputeAlignment);
    a.set_heuristic(HeuristicStrategy::XDrop {
        xdrop: xdrop.0,
        steps_between_cutoffs: xdrop.1,
    });
    a.set_compat(compat);
    let score = a.align_end2end(pat, txt);
    (score, a.cigar().clone())
}

/// C mode reproduces the invalid CIGAR (when `c_invalid`); corrected mode
/// is valid, with the same score. The CIGAR can rescore below that score
/// (never above): with pruning, the WFA path can open a gap right after
/// another, which the CIGAR merges into one (as in C's valid CIGARs).
fn check<const N: usize>(
    penalties: WavefrontPenalties,
    xdrop: (i32, i32),
    c_invalid: bool,
    rescore: impl Fn(&Cigar) -> i32,
) {
    let (pat, txt) = pair();
    let (c_score, c_cigar) = run::<N>(penalties.clone(), xdrop, Compat::C, &pat, &txt);
    assert_eq!(c_cigar.check_alignment(&pat, &txt).is_err(), c_invalid);
    let (score, cigar) = run::<N>(penalties, xdrop, Compat::Corrected, &pat, &txt);
    assert_eq!(score, c_score);
    cigar.check_alignment(&pat, &txt).unwrap();
    assert!(-rescore(&cigar) <= score);
}

#[test]
fn issue_21_affine() {
    let p = AffinePenalties {
        match_: 0,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 2,
    };
    for xdrop in [(50, 10), (100, 5)] {
        check::<3>(WavefrontPenalties::new_affine(p), xdrop, true, |c| {
            c.score_gap_affine(&p)
        });
    }
}

#[test]
fn issue_21_linear() {
    let p = LinearPenalties {
        match_: 0,
        mismatch: 4,
        indel: 2,
    };
    for xdrop in [(50, 10), (100, 5)] {
        check::<1>(WavefrontPenalties::new_linear(p), xdrop, true, |c| {
            c.score_gap_linear(&p)
        });
    }
}

#[test]
fn issue_21_affine2p() {
    let p = Affine2pPenalties {
        match_: 0,
        mismatch: 4,
        gap_opening1: 6,
        gap_extension1: 2,
        gap_opening2: 24,
        gap_extension2: 1,
    };
    // 50,10 is already valid in C for this pair.
    for (xdrop, c_invalid) in [((50, 10), false), ((100, 5), true)] {
        check::<5>(WavefrontPenalties::new_affine2p(p), xdrop, c_invalid, |c| {
            c.score_gap_affine2p(&p)
        });
    }
}
