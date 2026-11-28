//! X-drop with `steps_between_cutoffs > 1` (HPCBio/wfa2lib-rs#21): the
//! backtrace read wavefronts the cutoff trimmed to empty, giving invalid
//! CIGARs. C fixed it in 0285b52 (Marco-Sola-Lab/WFA2-lib#146); both compat
//! modes are valid and agree.

use wfa2lib_rs::aligner::{AlignmentScope, Compat, WavefrontAligner};
use wfa2lib_rs::cigar::Cigar;
use wfa2lib_rs::heuristic::HeuristicStrategy;
use wfa2lib_rs::penalties::{
    Affine2pPenalties, AffinePenalties, LinearPenalties, WavefrontPenalties,
};

/// `wfa.utest.seq` pair 305 (78 × 100 bp).
fn pair() -> (Vec<u8>, Vec<u8>) {
    let data = std::fs::read_to_string("tests/data/issue21_xdrop.seq").unwrap();
    let mut lines = data.lines();
    let pat = lines
        .next()
        .unwrap()
        .trim_start_matches('>')
        .as_bytes()
        .to_vec();
    let txt = lines
        .next()
        .unwrap()
        .trim_start_matches('<')
        .as_bytes()
        .to_vec();
    (pat, txt)
}

/// Align with X-drop in the given mode (BiWFA if `biwfa`); returns the
/// score and the CIGAR.
fn run<const N: usize>(
    penalties: WavefrontPenalties,
    xdrop: (i32, i32),
    compat: Compat,
    biwfa: bool,
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
    let score = if biwfa {
        a.align_biwfa(pat, txt)
    } else {
        a.align_end2end(pat, txt)
    };
    (score, a.cigar().clone())
}

/// Both modes give the same valid CIGAR and score, unidirectional and BiWFA.
/// The CIGAR scores at least as well as that score says (it can do better):
/// with pruning, the WFA path can open a gap right after another, which the
/// CIGAR merges into one.
fn check<const N: usize>(penalties: WavefrontPenalties, rescore: impl Fn(&Cigar) -> i32) {
    let xdrop = (20, 2);
    let (pat, txt) = pair();
    for biwfa in [false, true] {
        let (c_score, c_cigar) = run::<N>(penalties.clone(), xdrop, Compat::C, biwfa, &pat, &txt);
        c_cigar.check_alignment(&pat, &txt).unwrap();
        let (score, cigar) = run::<N>(
            penalties.clone(),
            xdrop,
            Compat::Corrected,
            biwfa,
            &pat,
            &txt,
        );
        assert_eq!(score, c_score, "biwfa={biwfa}");
        assert_eq!(
            cigar.to_string_rle(true),
            c_cigar.to_string_rle(true),
            "biwfa={biwfa}"
        );
        let sw = penalties.wf_score_to_sw_score(pat.len() as i32, txt.len() as i32, score);
        assert!(
            rescore(&cigar) >= sw,
            "biwfa={biwfa}: {} < {sw}",
            rescore(&cigar)
        );
    }
}

#[test]
fn issue_21_affine() {
    let p = AffinePenalties {
        match_: -3,
        mismatch: 5,
        gap_opening: 3,
        gap_extension: 2,
    };
    check::<3>(WavefrontPenalties::new_affine(p), |c| {
        c.score_gap_affine(&p)
    });
}

#[test]
fn issue_21_linear() {
    let p = LinearPenalties {
        match_: -1,
        mismatch: 3,
        indel: 2,
    };
    check::<1>(WavefrontPenalties::new_linear(p), |c| {
        c.score_gap_linear(&p)
    });
}

#[test]
fn issue_21_affine2p() {
    let p = Affine2pPenalties {
        match_: -1,
        mismatch: 4,
        gap_opening1: 6,
        gap_extension1: 2,
        gap_opening2: 24,
        gap_extension2: 1,
    };
    check::<5>(WavefrontPenalties::new_affine2p(p), |c| {
        c.score_gap_affine2p(&p)
    });
}
