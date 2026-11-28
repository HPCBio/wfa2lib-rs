//! Singletrack backtrace (C a823d54): high-memory gap-affine/affine2p
//! alignments trace back through the M wavefronts only. Its CIGARs can
//! differ from the classic backtrace's in tie-breaks, but must be valid
//! alignments with the same score. Med memory (classic piggyback backtrace)
//! is the reference here; the differential suite checks C parity.

use wfa2lib_rs::aligner::{AlignmentScope, MemoryMode, WavefrontAligner};
use wfa2lib_rs::cigar::Cigar;
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
    /// A random sequence and a copy with mismatches and indels (some long).
    fn pair(&mut self, max_len: u64) -> (Vec<u8>, Vec<u8>) {
        let len = 1 + self.below(max_len);
        let pattern: Vec<u8> = (0..len).map(|_| self.base()).collect();
        let mut text = Vec::with_capacity(len + 16);
        let mut i = 0;
        while i < pattern.len() {
            match self.below(40) {
                0 => text.push(self.base()), // mismatch (maybe)
                1 => {
                    let max = if self.below(4) == 0 { 12 } else { 3 };
                    let n = 1 + self.below(max);
                    text.extend((0..n).map(|_| self.base())); // insertion
                    text.push(pattern[i]);
                }
                2 => {
                    let max = if self.below(4) == 0 { 12 } else { 3 };
                    i += self.below(max); // deletion
                }
                _ => text.push(pattern[i]),
            }
            i += 1;
        }
        (pattern, text)
    }
}

fn penalty_sets() -> Vec<(&'static str, WavefrontPenalties)> {
    let affine = |match_, mismatch, gap_opening, gap_extension| {
        WavefrontPenalties::new_affine(AffinePenalties {
            match_,
            mismatch,
            gap_opening,
            gap_extension,
        })
    };
    let affine2p = |match_| {
        WavefrontPenalties::new_affine2p(Affine2pPenalties {
            match_,
            mismatch: 4,
            gap_opening1: 6,
            gap_extension1: 2,
            gap_opening2: 24,
            gap_extension2: 1,
        })
    };
    vec![
        ("affine 0,4,6,2", affine(0, 4, 6, 2)),
        ("affine 0,1,2,1", affine(0, 1, 2, 1)),
        ("affine -3,5,3,2", affine(-3, 5, 3, 2)),
        ("affine2p 0", affine2p(0)),
        ("affine2p -1", affine2p(-1)),
    ]
}

fn aligner<const N: usize>(pen: &WavefrontPenalties, mode: MemoryMode) -> WavefrontAligner<N> {
    let mut a = WavefrontAligner::<N>::new(pen.clone());
    a.set_memory_mode(mode);
    a.set_alignment_scope(AlignmentScope::ComputeAlignment);
    a
}

/// The classic score of `cigar` under `pen`'s penalties.
fn rescore(pen: &WavefrontPenalties, cigar: &Cigar) -> i32 {
    match (&pen.affine_penalties, &pen.affine2p_penalties) {
        (Some(p), _) => cigar.score_gap_affine(p),
        (_, Some(p)) => cigar.score_gap_affine2p(p),
        _ => unreachable!(),
    }
}

/// End-to-end: valid CIGAR scoring exactly the alignment's score, which
/// equals the classic backtrace's (med); lambda input gives the same CIGAR.
fn check_end2end<const N: usize>(name: &str, pen: &WavefrontPenalties, p: &[u8], t: &[u8]) {
    let ctx = || {
        format!(
            "{name} {:?} {:?}",
            String::from_utf8_lossy(p),
            String::from_utf8_lossy(t)
        )
    };
    let mut high = aligner::<N>(pen, MemoryMode::High);
    let score = high.align_end2end(p, t);
    let cigar = high.cigar().clone();
    cigar
        .check_alignment(p, t)
        .unwrap_or_else(|e| panic!("{e}: {}", ctx()));
    let mut med = aligner::<N>(pen, MemoryMode::Med);
    assert_eq!(med.align_end2end(p, t), score, "{}", ctx());
    let sw = if pen.match_ < 0 {
        pen.wf_score_to_sw_score(p.len() as i32, t.len() as i32, score)
    } else {
        -score
    };
    assert_eq!(rescore(pen, &cigar), sw, "{}", ctx());

    let (pl, tl) = (p.to_vec(), t.to_vec());
    let mut lambda = aligner::<N>(pen, MemoryMode::High);
    let lscore = lambda.align_lambda(
        Box::new(move |i, j| pl[i as usize] == tl[j as usize]),
        p.len() as i32,
        t.len() as i32,
    );
    assert_eq!(lscore, score, "lambda {}", ctx());
    assert_eq!(
        lambda.cigar().to_string_rle(true),
        cigar.to_string_rle(true),
        "lambda {}",
        ctx()
    );
}

/// Ends-free (begin and end gaps): valid CIGAR, score equal to the classic
/// backtrace's (med).
fn check_endsfree<const N: usize>(name: &str, pen: &WavefrontPenalties, p: &[u8], t: &[u8]) {
    let ctx = || {
        format!(
            "{name} ends-free {:?} {:?}",
            String::from_utf8_lossy(p),
            String::from_utf8_lossy(t)
        )
    };
    let free = |len: usize| (len / 10) as i32;
    let (pf, tf) = (free(p.len()), free(t.len()));
    let mut high = aligner::<N>(pen, MemoryMode::High);
    high.set_alignment_free_ends(pf, pf, tf, tf);
    let score = high.align_endsfree(p, t);
    high.cigar()
        .check_alignment(p, t)
        .unwrap_or_else(|e| panic!("{e}: {}", ctx()));
    let mut med = aligner::<N>(pen, MemoryMode::Med);
    med.set_alignment_free_ends(pf, pf, tf, tf);
    assert_eq!(med.align_endsfree(p, t), score, "{}", ctx());
}

#[test]
fn singletrack_random_pairs() {
    let mut rng = Lcg(0x5149_7e7a);
    for _ in 0..300 {
        let (p, t) = rng.pair(300);
        for (name, pen) in penalty_sets() {
            if pen.affine_penalties.is_some() {
                check_end2end::<3>(name, &pen, &p, &t);
                check_endsfree::<3>(name, &pen, &p, &t);
            } else {
                check_end2end::<5>(name, &pen, &p, &t);
                check_endsfree::<5>(name, &pen, &p, &t);
            }
        }
    }
}

#[test]
fn singletrack_edge_cases() {
    // Identical, one-sided gaps only, all mismatches, single bases.
    let cases: [(&[u8], &[u8]); 7] = [
        (b"ACGTACGT", b"ACGTACGT"),
        (b"ACGTACGT", b"ACGT"),
        (b"ACGT", b"ACGTACGT"),
        (b"AAAAAAAA", b"CCCCCCCC"),
        (b"A", b"C"),
        (b"A", b"A"),
        (b"ACGTTTTTTTTTTACGT", b"ACGTACGT"),
    ];
    for (p, t) in cases {
        for (name, pen) in penalty_sets() {
            if pen.affine_penalties.is_some() {
                check_end2end::<3>(name, &pen, p, t);
            } else {
                check_end2end::<5>(name, &pen, p, t);
            }
        }
    }
}
