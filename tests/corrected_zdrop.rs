//! `Compat::Corrected` Z-drop ends (HPCBio/wfa2lib-rs#22): C backtraces the
//! best cell at the current score, where it needn't exist (high mode: empty
//! or unchecked CIGARs; med/low: a dead cell). Corrected mode backtraces it
//! at the score it was found at, so every memory mode gives a valid partial
//! alignment with the same score.

use wfa2lib_rs::aligner::{AlignmentScope, Compat, MemoryMode, WavefrontAligner};
use wfa2lib_rs::cigar::Cigar;
use wfa2lib_rs::heuristic::HeuristicStrategy;
use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

/// `long.seq` pair 1 (tests/differential/data), 2 kb, 8% error.
fn pair() -> (Vec<u8>, Vec<u8>) {
    let data = std::fs::read_to_string("tests/data/issue22_zdrop.seq").unwrap();
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

const PEN: AffinePenalties = AffinePenalties {
    match_: -3,
    mismatch: 5,
    gap_opening: 3,
    gap_extension: 2,
};

/// Align with Z-drop `100,5`; returns the CIGAR.
fn run(compat: Compat, memory: MemoryMode) -> Cigar {
    let (pat, txt) = pair();
    let mut a = WavefrontAligner::<3>::new(WavefrontPenalties::new_affine(PEN));
    a.set_alignment_scope(AlignmentScope::ComputeAlignment);
    a.set_heuristic(HeuristicStrategy::ZDrop {
        zdrop: 100,
        steps_between_cutoffs: 5,
    });
    a.set_memory_mode(memory);
    a.set_compat(compat);
    a.align_end2end(&pat, &txt);
    a.cigar().clone()
}

/// Whether the CIGAR is a valid partial alignment: every `M` a match and
/// every `X` a mismatch, within both sequences (trailing bases may be left).
fn valid_prefix(cigar: &Cigar, pat: &[u8], txt: &[u8]) -> bool {
    let (mut i, mut j) = (0, 0);
    for &op in cigar.operations_slice() {
        match op {
            b'M' | b'X' => {
                if i >= pat.len() || j >= txt.len() || (pat[i] == txt[j]) != (op == b'M') {
                    return false;
                }
                i += 1;
                j += 1;
            }
            b'D' => i += 1,
            b'I' => j += 1,
            _ => return false,
        }
    }
    i <= pat.len() && j <= txt.len()
}

#[test]
fn issue_22_c_mode_high_is_invalid() {
    let (pat, txt) = pair();
    let cigar = run(Compat::C, MemoryMode::High);
    assert_eq!(cigar.to_string_rle(true), "2I58M");
    assert!(!valid_prefix(&cigar, &pat, &txt));
}

#[test]
fn issue_22_corrected_modes_agree() {
    let (pat, txt) = pair();
    let high = run(Compat::Corrected, MemoryMode::High);
    assert!(!high.operations_slice().is_empty());
    let score = high.score_gap_affine(&PEN);
    for memory in [MemoryMode::High, MemoryMode::Med, MemoryMode::Low] {
        let cigar = run(Compat::Corrected, memory);
        assert!(valid_prefix(&cigar, &pat, &txt), "{memory:?}");
        assert_eq!(cigar.score_gap_affine(&PEN), score, "{memory:?}");
    }
}

/// The best cell's backtrace is kept through backtrace-buffer compaction:
/// a 3 kb prefix with 5% mismatches, then 20 kb of unrelated sequence, so
/// the best cell is found early and the Z-drop comes after score 12000,
/// past two compaction probes (every 6000 scores). Without the pin its
/// block index goes stale.
#[test]
fn issue_22_corrected_survives_compaction() {
    let mut x: u64 = 7;
    let mut base = || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        b"ACGT"[(x >> 33) as usize % 4]
    };
    let prefix: Vec<u8> = (0..3000).map(|_| base()).collect();
    let mut pat = prefix.clone();
    let mut txt: Vec<u8> = prefix
        .iter()
        .enumerate()
        .map(|(i, &c)| match (i % 20, c) {
            (7, b'A') => b'C',
            (7, _) => b'A',
            _ => c,
        })
        .collect();
    pat.extend((0..20000).map(|_| base()));
    txt.extend((0..20000).map(|_| base()));
    let pen = AffinePenalties {
        match_: -1,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 2,
    };
    let align = |memory: MemoryMode, compact: bool| {
        let mut a = WavefrontAligner::<3>::new(WavefrontPenalties::new_affine(pen));
        a.set_alignment_scope(AlignmentScope::ComputeAlignment);
        a.set_heuristic(HeuristicStrategy::ZDrop {
            zdrop: 4000,
            steps_between_cutoffs: 1,
        });
        a.set_memory_mode(memory);
        if compact {
            a.set_max_memory_compact(0);
        }
        a.set_compat(Compat::Corrected);
        let score = a.align_end2end(&pat, &txt);
        assert!(score > 12000);
        a.cigar().clone()
    };
    let high = align(MemoryMode::High, false);
    assert!(valid_prefix(&high, &pat, &txt));
    for memory in [MemoryMode::Med, MemoryMode::Low] {
        let plain = align(memory, false);
        let compacted = align(memory, true);
        assert_eq!(
            compacted.to_string_rle(true),
            plain.to_string_rle(true),
            "{memory:?}"
        );
        assert_eq!(
            compacted.score_gap_affine(&pen),
            high.score_gap_affine(&pen),
            "{memory:?}"
        );
    }
}
