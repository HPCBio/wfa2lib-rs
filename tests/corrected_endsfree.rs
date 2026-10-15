//! `Compat::Corrected` ends-free alignment with a match reward
//! (HPCBio/wfa2lib-rs#37, #38, smarco/WFA2-lib#102), checked against a
//! dynamic-programming oracle.

use wfa2lib_rs::aligner::{AlignmentScope, Compat, MemoryMode, WavefrontAligner};
use wfa2lib_rs::heuristic::HeuristicStrategy;
use wfa2lib_rs::penalties::{
    Affine2pPenalties, AffinePenalties, LinearPenalties, WavefrontPenalties,
};

/// Penalties in C's convention (match <= 0); a gap of length `l` costs
/// `min(o1 + l * e1, o2 + l * e2)`.
#[derive(Debug, Clone, Copy)]
struct Pen {
    match_: i32,
    mismatch: i32,
    o1: i32,
    e1: i32,
    o2: i32,
    e2: i32,
}

impl Pen {
    fn linear(match_: i32, mismatch: i32, indel: i32) -> Self {
        Pen {
            match_,
            mismatch,
            o1: 0,
            e1: indel,
            o2: 0,
            e2: indel,
        }
    }
    fn affine(match_: i32, mismatch: i32, o: i32, e: i32) -> Self {
        Pen {
            match_,
            mismatch,
            o1: o,
            e1: e,
            o2: o,
            e2: e,
        }
    }
    fn gap(&self, l: i64) -> i64 {
        (self.o1 as i64 + l * self.e1 as i64).min(self.o2 as i64 + l * self.e2 as i64)
    }
}

/// Free ends: (pattern begin, pattern end, text begin, text end).
type Ends = (i32, i32, i32, i32);

/// Optimal ends-free penalty (lower is better) under WFA's model: an
/// alignment may start at (v, 0) with v <= pattern begin-free or (0, h)
/// with h <= text begin-free, and end at (v, |T|) with |P| - v <= pattern
/// end-free or (|P|, h) with |T| - h <= text end-free.
fn oracle(p: &Pen, pattern: &[u8], text: &[u8], (pb, pe, tb, te): Ends) -> i64 {
    const INF: i64 = i64::MAX / 4;
    let (n, m) = (pattern.len(), text.len());
    let pieces = [(p.o1 as i64, p.e1 as i64), (p.o2 as i64, p.e2 as i64)];
    let mut h = vec![vec![INF; m + 1]; n + 1];
    let mut ins = vec![vec![[INF; 2]; m + 1]; n + 1]; // consumes text
    let mut del = vec![vec![[INF; 2]; m + 1]; n + 1]; // consumes pattern
    for v in 0..=n {
        for j in 0..=m {
            let mut best = INF;
            if (j == 0 && v as i32 <= pb) || (v == 0 && j as i32 <= tb) {
                best = 0;
            }
            if v > 0 && j > 0 && h[v - 1][j - 1] < INF {
                let c = if pattern[v - 1] == text[j - 1] {
                    p.match_
                } else {
                    p.mismatch
                };
                best = best.min(h[v - 1][j - 1] + c as i64);
            }
            for (i, &(o, e)) in pieces.iter().enumerate() {
                if j > 0 {
                    let x = (h[v][j - 1] + o + e).min(ins[v][j - 1][i] + e);
                    ins[v][j][i] = x;
                    best = best.min(x);
                }
                if v > 0 {
                    let x = (h[v - 1][j] + o + e).min(del[v - 1][j][i] + e);
                    del[v][j][i] = x;
                    best = best.min(x);
                }
            }
            h[v][j] = best;
        }
    }
    let mut best = INF;
    for (v, row) in h.iter().enumerate() {
        if (n - v) as i32 <= pe {
            best = best.min(row[m]);
        }
    }
    for (j, &x) in h[n].iter().enumerate() {
        if (m - j) as i32 <= te {
            best = best.min(x);
        }
    }
    best
}

/// Ends-free penalty of a CIGAR: up to the free-end limits, a leading or
/// trailing gap run is free; the rest of it is charged as a gap.
fn cigar_penalty(p: &Pen, ops: &[u8], (pb, pe, tb, te): Ends) -> i64 {
    let mut runs: Vec<(u8, i64)> = Vec::new();
    for &op in ops {
        match runs.last_mut() {
            Some((o, l)) if *o == op => *l += 1,
            _ => runs.push((op, 1)),
        }
    }
    let free = |op: u8, begin: bool| -> i64 {
        (match (op, begin) {
            (b'D', true) => pb,
            (b'I', true) => tb,
            (b'D', false) => pe,
            (b'I', false) => te,
            _ => 0,
        }) as i64
    };
    let last = runs.len().saturating_sub(1);
    let mut total = 0;
    for (i, &(op, l)) in runs.iter().enumerate() {
        let mut charged = l;
        if i == 0 {
            charged -= charged.min(free(op, true));
        }
        if i == last {
            charged -= charged.min(free(op, false));
        }
        total += match op {
            b'M' => l * p.match_ as i64,
            b'X' => l * p.mismatch as i64,
            _ if charged > 0 => p.gap(charged),
            _ => 0,
        };
    }
    total
}

/// Penalty kinds: 0 linear, 1 affine, 2 affine2p, 3 edit, 4 indel.
fn penalties(p: &Pen, kind: usize) -> WavefrontPenalties {
    match kind {
        3 => WavefrontPenalties::new_edit(),
        4 => WavefrontPenalties::new_indel(),
        0 => WavefrontPenalties::new_linear(LinearPenalties {
            match_: p.match_,
            mismatch: p.mismatch,
            indel: p.e1,
        }),
        1 => WavefrontPenalties::new_affine(AffinePenalties {
            match_: p.match_,
            mismatch: p.mismatch,
            gap_opening: p.o1,
            gap_extension: p.e1,
        }),
        _ => WavefrontPenalties::new_affine2p(Affine2pPenalties {
            match_: p.match_,
            mismatch: p.mismatch,
            gap_opening1: p.o1,
            gap_extension1: p.e1,
            gap_opening2: p.o2,
            gap_extension2: p.e2,
        }),
    }
}

/// Align with `compat`; returns the real penalty from the WFA
/// score, and the CIGAR when `cigar` is set.
#[allow(clippy::too_many_arguments)]
fn align_corrected<const N: usize>(
    wp: &WavefrontPenalties,
    pattern: &[u8],
    text: &[u8],
    ends: Ends,
    memory: MemoryMode,
    cigar: bool,
    compat: Compat,
) -> (i64, Option<Vec<u8>>) {
    let mut a = WavefrontAligner::<N>::new(wp.clone());
    a.set_compat(compat);
    a.set_memory_mode(memory);
    if cigar {
        a.set_alignment_scope(AlignmentScope::ComputeAlignment);
    }
    a.set_alignment_free_ends(ends.0, ends.1, ends.2, ends.3);
    let w = a.align_endsfree(pattern, text);
    // With a match reward, the WFA score is transformed; else it's the penalty.
    let swg = if wp.match_ < 0 {
        wp.wf_score_to_sw_score(pattern.len() as i32, text.len() as i32, w)
    } else {
        -w
    };
    if cigar {
        let c = a.cigar();
        c.check_alignment(pattern, text)
            .unwrap_or_else(|e| panic!("invalid CIGAR {}: {e}", c.to_string_rle(true)));
        (-swg as i64, Some(c.operations_slice().to_vec()))
    } else {
        (-swg as i64, None)
    }
}

#[allow(clippy::too_many_arguments)]
fn run_kind(
    kind: usize,
    wp: &WavefrontPenalties,
    pattern: &[u8],
    text: &[u8],
    ends: Ends,
    memory: MemoryMode,
    cigar: bool,
    compat: Compat,
) -> (i64, Option<Vec<u8>>) {
    match kind {
        0 | 3 | 4 => align_corrected::<1>(wp, pattern, text, ends, memory, cigar, compat),
        1 => align_corrected::<3>(wp, pattern, text, ends, memory, cigar, compat),
        _ => align_corrected::<5>(wp, pattern, text, ends, memory, cigar, compat),
    }
}

/// Check score-only and every memory mode's CIGAR against the oracle.
fn check(kind: usize, p: &Pen, pattern: &[u8], text: &[u8], ends: Ends) {
    let wp = penalties(p, kind);
    let want = oracle(p, pattern, text, ends);
    let ctx = || {
        format!(
            "kind={kind} {p:?} ends={ends:?} P={} T={}",
            String::from_utf8_lossy(pattern),
            String::from_utf8_lossy(text)
        )
    };
    let (score, _) = run_kind(
        kind,
        &wp,
        pattern,
        text,
        ends,
        MemoryMode::High,
        false,
        Compat::Corrected,
    );
    assert_eq!(score, want, "score-only: {}", ctx());
    for memory in [MemoryMode::High, MemoryMode::Med, MemoryMode::Low] {
        let (score, ops) = run_kind(
            kind,
            &wp,
            pattern,
            text,
            ends,
            memory,
            true,
            Compat::Corrected,
        );
        assert_eq!(score, want, "{memory:?} score: {}", ctx());
        let ops = ops.unwrap();
        assert_eq!(
            cigar_penalty(p, &ops, ends),
            want,
            "{memory:?} CIGAR {}: {}",
            String::from_utf8_lossy(&ops),
            ctx()
        );
    }
}

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
    fn below(&mut self, n: u64) -> i32 {
        (self.next() % n) as i32
    }
    /// A random sequence of length below `max_len`.
    fn seq(&mut self, alphabet: &[u8], max_len: u64) -> Vec<u8> {
        let len = self.below(max_len);
        (0..len)
            .map(|_| alphabet[self.below(alphabet.len() as u64) as usize])
            .collect()
    }
}

#[test]
fn issue_37_repeat_shift() {
    // C (and C-compatible mode) shift by one repeat unit: IIIMMMMMMMMDDDD, 40.
    let p = Pen::affine(-5, 4, 0, 8);
    let (pat, txt) = (b"ACGTACGTACGT", b"ACGACGTACGT");
    let ends = (12, 12, 11, 11);
    assert_eq!(oracle(&p, pat, txt, ends), -47);
    let wp = penalties(&p, 1);
    let (_, c_ops) = run_kind(1, &wp, pat, txt, ends, MemoryMode::High, true, Compat::C);
    assert_eq!(c_ops.unwrap(), b"IIIMMMMMMMMDDDD");
    let (score, ops) = run_kind(
        1,
        &wp,
        pat,
        txt,
        ends,
        MemoryMode::High,
        true,
        Compat::Corrected,
    );
    assert_eq!(score, -47);
    assert_eq!(ops.unwrap(), b"MMMDMMMMMMMM");
    check(1, &p, pat, txt, ends);
}

#[test]
fn issue_38_whole_length_free_ends() {
    // C (and C-compatible mode) return the all-gap alignment.
    let p = Pen::affine(-5, 4, 0, 8);
    let pat = b"AAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGT";
    let txt = b"AGTCGTGACTCTTTTTGCATGTACCTTATGAATAAGGACCGGCTAATTCCG";
    let (n, m) = (pat.len() as i32, txt.len() as i32);
    let ends = (n, n, m, m);
    let wp = penalties(&p, 1);
    let (_, c_ops) = run_kind(1, &wp, pat, txt, ends, MemoryMode::High, true, Compat::C);
    assert!(c_ops.unwrap().iter().all(|&op| op == b'I' || op == b'D'));
    assert_eq!(oracle(&p, pat, txt, ends), -78);
    check(1, &p, pat, txt, ends);
    check(1, &p, b"CAGA", b"GGGA", (4, 4, 4, 4));
}

#[test]
fn match_zero_identical_to_c() {
    let mut rng = Lcg(7);
    for _ in 0..300 {
        let kind = rng.below(3) as usize;
        let p = Pen {
            match_: 0,
            mismatch: 1 + rng.below(6),
            o1: rng.below(6),
            e1: 1 + rng.below(4),
            o2: rng.below(12),
            e2: 1 + rng.below(2),
        };
        let wp = penalties(&p, kind);
        let (pat, txt) = (rng.seq(b"ACGT", 16), rng.seq(b"ACGT", 16));
        let (n, m) = (pat.len() as u64 + 1, txt.len() as u64 + 1);
        let ends = (rng.below(n), rng.below(n), rng.below(m), rng.below(m));
        for memory in [MemoryMode::High, MemoryMode::Med] {
            let c = run_kind(kind, &wp, &pat, &txt, ends, memory, true, Compat::C);
            let fixed = run_kind(kind, &wp, &pat, &txt, ends, memory, true, Compat::Corrected);
            assert_eq!(c, fixed);
        }
    }
}

#[test]
fn random_against_oracle() {
    let mut rng = Lcg(42);
    for iter in 0..3000 {
        let kind = rng.below(3) as usize;
        let match_ = -1 - rng.below(6);
        let p = match kind {
            0 => Pen::linear(match_, 1 + rng.below(8), 1 + rng.below(8)),
            1 => Pen::affine(match_, 1 + rng.below(8), rng.below(8), 1 + rng.below(8)),
            _ => {
                let (o1, e1) = (rng.below(6), 2 + rng.below(6));
                Pen {
                    match_,
                    mismatch: 1 + rng.below(8),
                    o1,
                    e1,
                    o2: o1 + 2 + rng.below(10),
                    e2: 1,
                }
            }
        };
        // Small alphabets and repeats make periodic, shifted solutions likely.
        let alphabet: &[u8] = if iter % 2 == 0 { b"AC" } else { b"ACGT" };
        let base = rng.seq(alphabet, 6);
        let mut pat = rng.seq(alphabet, 20);
        let mut txt = rng.seq(alphabet, 20);
        if iter % 3 == 0 {
            pat = base.repeat(1 + rng.below(4) as usize);
            txt = base.repeat(1 + rng.below(4) as usize);
            if !txt.is_empty() {
                let i = rng.below(txt.len() as u64) as usize;
                txt.remove(i);
            }
        }
        let (n, m) = (pat.len() as i32, txt.len() as i32);
        let ends = match iter % 4 {
            0 => (n, n, m, m),
            1 => (0, n, 0, m), // extension-like
            _ => (
                rng.below(n as u64 + 1),
                rng.below(n as u64 + 1),
                rng.below(m as u64 + 1),
                rng.below(m as u64 + 1),
            ),
        };
        check(kind, &p, &pat, &txt, ends);
    }
}

/// Longer pairs (a shared core with edits, random flanks, unequal lengths),
/// where the first end found is far from optimal and the pruning of
/// diagonals that can't beat it trims most of each wavefront.
#[test]
fn random_long_against_oracle() {
    let mut rng = Lcg(7);
    for iter in 0..300 {
        let kind = rng.below(3) as usize;
        let match_ = -1 - rng.below(4);
        let p = match kind {
            0 => Pen::linear(match_, 2 + rng.below(6), 1 + rng.below(6)),
            1 => Pen::affine(match_, 2 + rng.below(6), rng.below(8), 1 + rng.below(4)),
            _ => Pen {
                match_,
                mismatch: 2 + rng.below(6),
                o1: 2,
                e1: 4,
                o2: 12,
                e2: 1,
            },
        };
        let core = rng.seq(b"ACGT", 120);
        let mut pat = Vec::new();
        for &c in &core {
            match rng.below(20) {
                0 => {}
                1 => pat.push(b"ACGT"[rng.below(4) as usize]),
                2 => pat.extend([c, b"ACGT"[rng.below(4) as usize]]),
                _ => pat.push(c),
            }
        }
        let mut txt = core;
        // Random flanks, each present half the time
        for s in [&mut pat, &mut txt] {
            if rng.below(2) == 0 {
                let f = rng.seq(b"ACGT", 40);
                s.splice(0..0, f);
            }
            if rng.below(2) == 0 {
                let f = rng.seq(b"ACGT", 40);
                s.extend(f);
            }
        }
        let (n, m) = (pat.len() as i32, txt.len() as i32);
        let ends = match iter % 3 {
            0 => (n, n, m, m),
            1 => (0, n, 0, m), // extension-like
            _ => (
                rng.below(n as u64 + 1),
                rng.below(n as u64 + 1),
                rng.below(m as u64 + 1),
                rng.below(m as u64 + 1),
            ),
        };
        check(kind, &p, &pat, &txt, ends);
    }
}

/// Corrected mode with `BandedStatic` where the end diagonal lies outside
/// the band: the band and the pruning empty every wavefront long before the
/// optimum, a begin-free seed outside the band at score 754, is reached.
/// Null steps must not end the search while a candidate is pending (they
/// did, returning a worse end).
#[test]
fn banded_static_pending_seeds() {
    let data = std::fs::read_to_string("tests/data/corrected_banded_seeds.seq").unwrap();
    let mut lines = data.lines();
    let pat = lines.next().unwrap().trim_start_matches('>').as_bytes();
    let txt = lines.next().unwrap().trim_start_matches('<').as_bytes();
    let p = Pen::linear(-2, 4, 3);
    let (n, m) = (pat.len() as i32, txt.len() as i32);
    let ends = (n, n, m, m);
    let want = oracle(&p, pat, txt, ends);
    let wp = penalties(&p, 0);
    for memory in [MemoryMode::High, MemoryMode::Med, MemoryMode::Low] {
        let mut a = WavefrontAligner::<1>::new(wp.clone());
        a.set_compat(Compat::Corrected);
        a.set_memory_mode(memory);
        a.set_heuristic(HeuristicStrategy::BandedStatic {
            min_k: -10,
            max_k: 200,
        });
        a.set_alignment_scope(AlignmentScope::ComputeAlignment);
        a.set_alignment_free_ends(ends.0, ends.1, ends.2, ends.3);
        a.align_endsfree(pat, txt);
        let c = a.cigar();
        c.check_alignment(pat, txt).unwrap();
        assert_eq!(
            cigar_penalty(&p, c.operations_slice(), ends),
            want,
            "{memory:?}"
        );
    }
}

/// Med/low with the backtrace buffer compacted at every probe (every 6000
/// scores), while an end candidate found before the compaction is pending:
/// its chain must survive (it is pinned) and give the same alignment as
/// high memory.
#[test]
fn compaction_keeps_pending_candidate() {
    let mut rng = Lcg(3);
    let p = Pen::affine(-1, 6, 2, 2);
    let wp = penalties(&p, 1);
    // Pattern ~ text prefix with ~20% mismatches; a random text tail is left
    // to the free end, so the candidate (pattern consumed) is found well
    // before its total is reached.
    let mut pat = Vec::new();
    let mut txt = Vec::new();
    for _ in 0..3000 {
        let c = b"ACGT"[rng.below(4) as usize];
        pat.push(c);
        txt.push(if rng.below(5) == 0 {
            b"ACGT"[rng.below(4) as usize]
        } else {
            c
        });
    }
    txt.extend((0..1500).map(|_| b"ACGT"[rng.below(4) as usize]));
    let ends = (0, 0, 0, 1500);
    let (want, ops) = run_kind(
        1,
        &wp,
        &pat,
        &txt,
        ends,
        MemoryMode::High,
        true,
        Compat::Corrected,
    );
    assert_eq!(cigar_penalty(&p, &ops.unwrap(), ends), want);
    for memory in [MemoryMode::Med, MemoryMode::Low] {
        let mut a = WavefrontAligner::<3>::new(wp.clone());
        a.set_compat(Compat::Corrected);
        a.set_memory_mode(memory);
        a.set_max_memory_compact(0);
        a.set_alignment_scope(AlignmentScope::ComputeAlignment);
        a.set_alignment_free_ends(ends.0, ends.1, ends.2, ends.3);
        a.align_endsfree(&pat, &txt);
        let c = a.cigar();
        c.check_alignment(&pat, &txt).unwrap();
        assert_eq!(
            cigar_penalty(&p, c.operations_slice(), ends),
            want,
            "{memory:?}"
        );
    }
}

/// HPCBio/wfa2lib-rs#38 pair B: two full-length 16S HiFi reads with
/// whole-length free ends. C returns the all-gap alignment.
#[test]
fn issue_38_pair_b_16s() {
    let data = std::fs::read_to_string("tests/data/issue38_pair_b.seq").unwrap();
    let mut lines = data.lines();
    let pat = lines.next().unwrap().trim_start_matches('>').as_bytes();
    let txt = lines.next().unwrap().trim_start_matches('<').as_bytes();
    let p = Pen::affine(-5, 4, 0, 8);
    let (n, m) = (pat.len() as i32, txt.len() as i32);
    let ends = (n, n, m, m);
    let wp = penalties(&p, 1);
    let (_, c_ops) = run_kind(1, &wp, pat, txt, ends, MemoryMode::High, true, Compat::C);
    assert!(c_ops.unwrap().iter().all(|&op| op == b'I' || op == b'D'));
    check(1, &p, pat, txt, ends);
}

/// Edit/indel oracle penalties: unit costs; indel forbids mismatches.
fn edit_pen(kind: usize) -> Pen {
    Pen::linear(0, if kind == 3 { 1 } else { 1000 }, 1)
}

/// HPCBio/wfa2lib-rs#18 (smarco/WFA2-lib#142): edit/indel ends-free with
/// begin gaps returns an invalid CIGAR in C; corrected mode is valid.
#[test]
fn issue_18_edit_begin_gaps() {
    let (pat, txt) = (b"ACGTTCGTACGA".as_slice(), b"CCACGTACGTACGA".as_slice());
    let ends = (0, 0, 3, 0);
    let mut a = WavefrontAligner::<1>::new(WavefrontPenalties::new_edit());
    a.set_alignment_scope(AlignmentScope::ComputeAlignment);
    a.set_alignment_free_ends(ends.0, ends.1, ends.2, ends.3);
    assert_eq!(a.align_endsfree(pat, txt), 1);
    assert_eq!(a.cigar().to_string_rle(true), "2I12M");
    assert!(a.cigar().check_alignment(pat, txt).is_err());
    check(3, &edit_pen(3), pat, txt, ends);
}

#[test]
fn edit_indel_random_against_oracle() {
    let mut rng = Lcg(11);
    for iter in 0..3000 {
        let kind = 3 + (iter % 2);
        let alphabet: &[u8] = if iter % 3 == 0 { b"AC" } else { b"ACGT" };
        let pat = rng.seq(alphabet, 20);
        let txt = rng.seq(alphabet, 20);
        let (n, m) = (pat.len() as u64 + 1, txt.len() as u64 + 1);
        let ends = (rng.below(n), rng.below(n), rng.below(m), rng.below(m));
        check(kind, &edit_pen(kind), &pat, &txt, ends);
    }
}
