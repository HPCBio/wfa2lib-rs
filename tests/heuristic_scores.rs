//! With heuristics, pruning can leave a path that closes a gap and opens one
//! of the same type, which the CIGAR merges, so the CIGAR scores better than
//! the WFA score. Since C f4e0d0d/5713de8, a heuristic gap-affine/affine2p
//! alignment (and a heuristic BiWFA alignment of any metric) reports its
//! CIGAR's score. In Rust: the returned score and `cigar().score`, in WFA
//! units. A Z-dropped alignment's max-trimmed CIGAR scores its kept part.

use wfa2lib_rs::aligner::{AlignStatus, AlignmentScope, MemoryMode, WavefrontAligner};
use wfa2lib_rs::cigar::Cigar;
use wfa2lib_rs::heuristic::HeuristicStrategy;
use wfa2lib_rs::penalties::{Affine2pPenalties, AffinePenalties, WavefrontPenalties};

/// Small deterministic generator (no extra crates).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> usize {
        (self.next() % n) as usize
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4)]
    }
    /// A random sequence and a copy with ~15% mismatches and indels.
    fn pair(&mut self, max_len: u64) -> (Vec<u8>, Vec<u8>) {
        let len = 20 + self.below(max_len);
        let pattern: Vec<u8> = (0..len).map(|_| self.base()).collect();
        let mut text = Vec::with_capacity(len + 16);
        for &b in &pattern {
            match self.below(20) {
                0 => text.push(self.base()),
                1 => text.extend([self.base(), b]),
                2 => {}
                _ => text.push(b),
            }
        }
        (pattern, text)
    }
}

fn penalty_sets() -> Vec<WavefrontPenalties> {
    let affine = |match_, mismatch, gap_opening, gap_extension| {
        WavefrontPenalties::new_affine(AffinePenalties {
            match_,
            mismatch,
            gap_opening,
            gap_extension,
        })
    };
    vec![
        affine(0, 4, 6, 2),
        affine(-3, 5, 3, 2),
        WavefrontPenalties::new_affine2p(Affine2pPenalties {
            match_: 0,
            mismatch: 4,
            gap_opening1: 6,
            gap_extension1: 2,
            gap_opening2: 24,
            gap_extension2: 1,
        }),
    ]
}

fn heuristics() -> Vec<HeuristicStrategy> {
    vec![
        HeuristicStrategy::WfAdaptive {
            min_wavefront_length: 10,
            max_distance_threshold: 50,
            steps_between_cutoffs: 1,
        },
        HeuristicStrategy::XDrop {
            xdrop: 20,
            steps_between_cutoffs: 1,
        },
        HeuristicStrategy::ZDrop {
            zdrop: 20,
            steps_between_cutoffs: 1,
        },
        HeuristicStrategy::BandedAdaptive {
            min_k: -5,
            max_k: 5,
            steps_between_cutoffs: 1,
        },
    ]
}

/// The WFA score of `cigar`, spanning `v` pattern and `h` text bases.
fn wfa_score(pen: &WavefrontPenalties, cigar: &Cigar, v: i32, h: i32) -> i32 {
    let classic = match (&pen.affine_penalties, &pen.affine2p_penalties) {
        (Some(p), _) => cigar.score_gap_affine(p),
        (_, Some(p)) => cigar.score_gap_affine2p(p),
        _ => unreachable!(),
    };
    if pen.match_ < 0 {
        -pen.match_ * (v + h) - 2 * classic
    } else {
        -classic
    }
}

fn check<const N: usize>(pen: &WavefrontPenalties, p: &[u8], t: &[u8]) {
    let (plen, tlen) = (p.len() as i32, t.len() as i32);
    for h in heuristics() {
        let zdrop = matches!(h, HeuristicStrategy::ZDrop { .. });
        for memory in [MemoryMode::High, MemoryMode::Med] {
            let ctx = format!("{h:?} {memory:?} {:?} {:?}", p.len(), t.len());
            let mut a = WavefrontAligner::<N>::new(pen.clone());
            a.set_memory_mode(memory);
            a.set_alignment_scope(AlignmentScope::ComputeAlignment);
            a.set_heuristic(h.clone());
            let score = a.align_end2end(p, t);
            let cigar = a.cigar();
            match a.status() {
                AlignStatus::Completed => {
                    let expected = wfa_score(pen, cigar, plen, tlen);
                    assert_eq!(score, expected, "{ctx}");
                    assert_eq!(cigar.score, expected, "{ctx}");
                }
                AlignStatus::Partial if zdrop && !cigar.operations_slice().is_empty() => {
                    // Max-trimmed: the kept part's score.
                    let expected = wfa_score(pen, cigar, cigar.end_v, cigar.end_h);
                    assert_eq!(cigar.score, expected, "{ctx}");
                }
                _ => {}
            }
        }
        if zdrop || matches!(h, HeuristicStrategy::XDrop { .. }) {
            continue; // Not supported with BiWFA (C 7110ee3)
        }
        let mut a = WavefrontAligner::<N>::new(pen.clone());
        a.set_heuristic(h.clone());
        let score = a.align_biwfa(p, t);
        if a.status() == AlignStatus::Completed {
            let expected = wfa_score(pen, a.cigar(), plen, tlen);
            assert_eq!(score, expected, "BiWFA {h:?}");
            assert_eq!(a.cigar().score, expected, "BiWFA {h:?}");
        }
    }
}

#[test]
fn heuristic_scores_are_cigar_scores() {
    let mut rng = Lcg(0x6865_7572);
    for i in 0..120 {
        // Long enough pairs for BiWFA to split (above its fallback length).
        let (p, t) = rng.pair(if i % 4 == 0 { 2000 } else { 300 });
        for pen in penalty_sets() {
            if pen.affine_penalties.is_some() {
                check::<3>(&pen, &p, &t);
            } else {
                check::<5>(&pen, &p, &t);
            }
        }
    }
}

/// A pair whose high-memory X-drop CIGAR merges two deletions (`2D3M3D`):
/// the WFA score is 228, the CIGAR's 222 (found by random search). Med's
/// backtrace keeps them apart (`2D3M1D1M2D`, as C's), scoring 228.
#[test]
fn xdrop_merged_gaps_report_cigar_score() {
    let p = b"TATAGGCGATTTCGCGTTTGAGGGAAACTGACCTCATATAATACACATCAAGCCGGTGCTTGTTCTCCACATTCATCAGCCTGGCCGAACGGCTTTCAATG";
    let t = b"TATAGGCGATTTCGCGTTTGACGGAAACTCGACCTCATCATATACACATCAAGTATCGCTGTGTCCACATTCATACCTGGCCGAACGGCGTCTCAAATG";
    let pen = &penalty_sets()[1]; // -3,5,3,2
    for (memory, expected) in [(MemoryMode::High, 222), (MemoryMode::Med, 228)] {
        let mut a = WavefrontAligner::<3>::new(pen.clone());
        a.set_memory_mode(memory);
        a.set_alignment_scope(AlignmentScope::ComputeAlignment);
        a.set_heuristic(HeuristicStrategy::XDrop {
            xdrop: 10,
            steps_between_cutoffs: 1,
        });
        assert_eq!(a.align_end2end(p, t), expected, "{memory:?}");
        assert_eq!(a.cigar().score, expected, "{memory:?}");
        assert_eq!(wfa_score(pen, a.cigar(), 101, 99), expected, "{memory:?}");
    }
}
