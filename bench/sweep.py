#!/usr/bin/env python3
"""Rust-vs-C benchmark sweep (HPCBio/wfa2lib-rs#31).

Builds Rust's align_benchmark and C WFA2-lib's at the pinned commit
(generic, as the differential harness builds it, and with -march=native),
generates deterministic datasets, and runs every configuration with each
tool. For each run it records the best-of-N alignment time (the tools'
`Time.Alignment`), max RSS, instructions (Linux `perf stat`, macOS
`time -l`; whole process, so input parsing and output are included), and
whether the output is identical to generic C's.

Usage (from the repository root):
    bench/sweep.py                     # build, generate data, run everything
    bench/sweep.py --quick             # 1/10 of the pairs, 1 run each
    bench/sweep.py --suites core span  # core, span, heur
    bench/sweep.py --pin 3             # Linux: run on CPU 3 (taskset)
    bench/sweep.py --offline           # node without network (see below)

Offline (e.g. a compute node): with network, clone this repository, clone
https://github.com/Marco-Sola-Lab/WFA2-lib into ./WFA2-lib, and run `cargo fetch`.
On the node, run with --offline. Builds happen on the node: Rust and the
native C build target the CPU they're built on, so they live in per-host
directories and are never reused from another machine (a login node).

Results go to target/bench-sweep/<host>/: results.tsv (one row per tool
and case) and summary.md (ratios vs C, machine info). Single-threaded:
run it on an otherwise idle machine (or an exclusive node) for stable
times. Needs Python 3.6+, cargo, make and a C compiler.
"""

import argparse
import hashlib
import os
import platform
import random
import re
import shutil
import socket
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ALG = {
    "edit": ["-a", "edit-wfa"],
    "indel": ["-a", "indel-wfa"],
    "linear": ["-a", "gap-linear-wfa"],
    "affine": ["-a", "gap-affine-wfa"],
    "affine2p": ["-a", "gap-affine2p-wfa"],
}
MODES = {
    "cigar": ["--wfa-memory", "high"],
    "score": ["--wfa-memory", "high", "--wfa-score-only"],
    "med": ["--wfa-memory", "med"],
    "low": ["--wfa-memory", "low"],
    "ultralow": ["--wfa-memory", "ultralow"],
}
# name: (pairs, length, error rate, extra indels (count, max length), length
# of random flanks added to both ends of both sequences)
DATASETS = {
    "s150_2": (50000, 150, 0.02, (0, 0), 0),
    "m1k_5": (20000, 1000, 0.05, (0, 0), 0),
    "m1k_20": (5000, 1000, 0.20, (0, 0), 0),
    "m1k_indel": (2000, 1000, 0.05, (5, 50), 0),
    "l10k_5": (300, 10000, 0.05, (0, 0), 0),
    "m2k_flank": (300, 2000, 0.05, (0, 0), 200),
}


def cases(suites):
    """(suite, name, args, dataset) for every configuration."""
    out = []
    if "core" in suites:
        for data in ["s150_2", "m1k_5", "m1k_20", "m1k_indel", "l10k_5"]:
            for alg in ALG:
                for mode in MODES:
                    out.append(("core", "%s %s %s" % (data, alg, mode), ALG[alg] + MODES[mode], data))
    if "span" in suites:
        pens = [
            ("edit", []),
            ("linear", []),
            ("affine", []),
            ("affine m-1", ["-g", "-1,4,6,2"]),
            ("affine2p", []),
            ("affine2p m-1", ["--affine2p-penalties", "-1,4,6,2,24,1"]),
        ]
        for data in ["m1k_5", "m2k_flank"]:
            for name, pen in pens:
                alg = ALG[name.split()[0]]
                for span in ["ends-free,0.1,0.1,0.1,0.1", "extension"]:
                    for mode in ["cigar", "score", "med"]:
                        out.append(("span", "%s %s %s %s" % (data, name, span.split(",")[0], mode),
                                    alg + pen + ["--wfa-span", span] + MODES[mode], data))
    if "heur" in suites:
        heur = [("wfa-adaptive", "10,50,1"), ("xdrop", "100,1"), ("zdrop", "100,1"),
                ("banded-static", "-50,50"), ("banded-adaptive", "-50,50,1")]
        for data in ["m1k_5", "m1k_20"]:
            for alg in ["edit", "affine", "affine2p"]:
                for h, par in heur:
                    if alg == "edit" and h in ("xdrop", "zdrop"):
                        continue  # gap-affine only, rejected by both tools
                    for mode in ["cigar", "score", "med", "ultralow"]:
                        out.append(("heur", "%s %s %s %s" % (data, alg, h, mode),
                                    ALG[alg] + ["--wfa-heuristic", h, "--wfa-heuristic-parameters", par] + MODES[mode],
                                    data))
    return out


def rust_args(args):
    """Rust's CLI (clap) reads a value starting with '-' as a flag, so pass
    `--opt=value` there; C's getopt would read the '=' into a short option's
    value, so C gets them separately."""
    out, i = [], 0
    while i < len(args):
        if args[i].startswith("-") and i + 1 < len(args) and args[i + 1].startswith("-") \
                and re.match(r"^-\d", args[i + 1]):
            out.append(args[i] + "=" + args[i + 1])
            i += 2
        else:
            out.append(args[i])
            i += 1
    return out


def generate(path, pairs, length, err, indels, flank, seed):
    """Pairs in the WFA2-lib format: a random pattern, and a text with
    `err * length` random mismatches/insertions/deletions, plus optional
    extra indels and random flanks. Deterministic for a given seed."""
    rng = random.Random(seed)
    acgt = "ACGT"
    with open(path, "w") as f:
        for _ in range(pairs):
            p = [rng.choice(acgt) for _ in range(length)]
            t = list(p)
            for _ in range(int(err * length)):
                pos = rng.randrange(len(t) + 1)
                kind = rng.randrange(3)
                if kind == 0 and pos < len(t):
                    t[pos] = rng.choice([c for c in acgt if c != t[pos]])
                elif kind == 1:
                    t.insert(pos, rng.choice(acgt))
                elif pos < len(t):
                    del t[pos]
            for _ in range(indels[0]):
                n = rng.randint(1, indels[1])
                pos = rng.randrange(len(t) + 1)
                if rng.randrange(2):
                    t[pos:pos] = [rng.choice(acgt) for _ in range(n)]
                else:
                    del t[pos:pos + n]
            if flank:
                fl = lambda: "".join(rng.choice(acgt) for _ in range(rng.randint(0, flank)))
                p = list(fl()) + p + list(fl())
                t = list(fl()) + t + list(fl())
            f.write(">" + "".join(p) + "\n<" + "".join(t) + "\n")


def machine_info():
    info = ["host: %s" % socket.gethostname(), "os: %s %s" % (platform.system(), platform.release())]
    cpu = platform.processor() or platform.machine()
    try:
        if sys.platform == "darwin":
            cpu = subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string"]).decode().strip()
        else:
            for line in open("/proc/cpuinfo"):
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
            flags = open("/proc/cpuinfo").read()
            info.append("simd: " + " ".join(f for f in ["avx2", "avx512f", "avx512bw"] if " %s" % f in flags))
    except Exception:
        pass
    info.insert(1, "cpu: %s" % cpu)
    for cmd in (["rustc", "--version"], ["cc", "--version"]):
        try:
            info.append("%s: %s" % (cmd[0], subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode().splitlines()[0]))
        except Exception:
            pass
    try:
        info.append("rust commit: %s" % subprocess.check_output(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"]).decode().strip())
    except Exception:
        pass
    info.append("load at start: %s" % " ".join("%.2f" % x for x in os.getloadavg()))
    return info


def perf_works():
    if not sys.platform.startswith("linux") or not shutil.which("perf"):
        return False
    try:
        r = subprocess.run(["perf", "stat", "-x,", "-e", "instructions:u", "true"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=30)
        return r.returncode == 0 and b"not supported" not in r.stderr and b"<not" not in r.stderr
    except Exception:
        return False


def run_once(binary, args, inp, out, pin, use_perf):
    """(alignment ms, instructions, max RSS MB, return code)."""
    cmd = [binary, "-i", inp, "-o", out] + args
    perf_out = out + ".perf"
    if use_perf:
        cmd = ["perf", "stat", "-x,", "-e", "instructions:u", "-o", perf_out, "--"] + cmd
    elif sys.platform == "darwin":
        cmd = ["/usr/bin/time", "-l"] + cmd
    if pin is not None:
        cmd = ["taskset", "-c", str(pin)] + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    text = p.stdout.read().decode(errors="replace")
    _, status, ru = os.wait4(p.pid, 0)
    p.returncode = os.waitstatus_to_exitcode(status) if hasattr(os, "waitstatus_to_exitcode") else (status >> 8)
    rss = ru.ru_maxrss / (1048576.0 if sys.platform == "darwin" else 1024.0)
    ms = None
    m = re.search(r"Time\.Alignment\s+([\d.]+)\s*(us|ms|s|m)\b", text)
    if m:
        ms = float(m.group(1)) * {"us": 0.001, "ms": 1.0, "s": 1000.0, "m": 60000.0}[m.group(2)]
    ins = None
    if use_perf and os.path.exists(perf_out):
        for line in open(perf_out):
            parts = line.split(",")
            if len(parts) > 2 and "instructions" in parts[2] and parts[0].strip().isdigit():
                ins = int(parts[0])
    else:
        m = re.search(r"(\d+)\s+instructions retired", text)
        if m:
            ins = int(m.group(1))
    return ms, ins, rss, p.returncode


def md5(path):
    if not os.path.exists(path):
        return "-"
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def host_name():
    return socket.gethostname().split(".")[0]


def rust_target_dir():
    """Per-host: Rust is built with target-cpu=native, which cargo doesn't
    track, so a binary built on another machine (e.g. a login node) would be
    reused here."""
    return os.path.join(ROOT, "target", "bench-sweep", "build-" + host_name())


def c_dirs():
    commit = open(os.path.join(ROOT, "tests", "differential", "C_COMMIT")).read().strip()
    base = os.path.join(ROOT, "target", "differential", "c-" + commit)
    # The native build is per host too (see build_c.sh)
    return commit, base, base + "-native-" + host_name()


def build(no_c_native, offline):
    def sh(cmd, env=None):
        print("+ " + " ".join(cmd), file=sys.stderr)
        lines = subprocess.check_output(cmd, cwd=ROOT, env=env).decode().strip().splitlines()
        return lines[-1] if lines else ""
    commit, _, _ = c_dirs()
    if offline and not os.environ.get("WFA2_LIB_REPO"):
        repo = os.path.join(ROOT, "WFA2-lib")
        ok = os.path.isdir(os.path.join(repo, ".git")) and subprocess.call(
            ["git", "-C", repo, "cat-file", "-e", commit + "^{commit}"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
        if not ok:
            sys.exit("--offline: needs a WFA2-lib git checkout containing %s at ./WFA2-lib (or "
                     "WFA2_LIB_REPO=<path>). With network: git clone https://github.com/Marco-Sola-Lab/WFA2-lib"
                     % commit)
    env = dict(os.environ, CARGO_TARGET_DIR=rust_target_dir())
    sh(["cargo", "build", "--release", "--bin", "align_benchmark"] + (["--offline"] if offline else []), env=env)
    tools = {"Rust": os.path.join(rust_target_dir(), "release", "align_benchmark")}
    tools["C"] = sh(["bash", "tests/differential/build_c.sh"])
    if not no_c_native:
        try:
            tools["C-native"] = sh(["bash", "tests/differential/build_c.sh"], env=dict(os.environ, C_NATIVE="1"))
        except subprocess.CalledProcessError:
            print("C -march=native build failed; skipping C-native", file=sys.stderr)
    return tools


def ratio(a, b):
    return "%.2f" % (a / b) if a and b else "—"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--suites", nargs="+", default=["core", "span", "heur"], choices=["core", "span", "heur"])
    ap.add_argument("--reps", type=int, default=None, help="runs per case and tool (default 3, --quick 1)")
    ap.add_argument("--quick", action="store_true", help="1/10 of the pairs and 1 run each")
    ap.add_argument("--filter", default=None, help="only cases whose name matches this regex")
    ap.add_argument("--pin", type=int, default=None, help="Linux: pin runs to this CPU (taskset)")
    ap.add_argument("--no-c-native", action="store_true", help="skip the -march=native C build")
    ap.add_argument("--no-build", action="store_true", help="use existing builds (from this host)")
    ap.add_argument("--offline", action="store_true",
                    help="no network: cargo --offline (run `cargo fetch` beforehand) and C from ./WFA2-lib")
    ap.add_argument("--out", default=None, help="output directory")
    a = ap.parse_args()
    reps = a.reps or (1 if a.quick else 3)
    host = host_name()
    outdir = a.out or os.path.join(ROOT, "target", "bench-sweep", host + ("-quick" if a.quick else ""))
    datadir = os.path.join(ROOT, "target", "bench-sweep", "data" + ("-quick" if a.quick else ""))
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(datadir, exist_ok=True)

    if a.no_build:
        _, c_generic, c_native = c_dirs()
        tools = {"Rust": os.path.join(rust_target_dir(), "release", "align_benchmark")}
        for name, d in (("C", c_generic), ("C-native", c_native)):
            path = os.path.join(d, "bin", "align_benchmark")
            if os.path.exists(path) and not (name == "C-native" and a.no_c_native):
                tools[name] = path
        missing = [t for t in ("Rust", "C") if not os.path.exists(tools.get(t, ""))]
        if missing:
            sys.exit("--no-build: no %s binary built on this host; run without --no-build" % " or ".join(missing))
    else:
        tools = build(a.no_c_native, a.offline)
    use_perf = perf_works()

    todo = cases(a.suites)
    if a.filter:
        todo = [c for c in todo if re.search(a.filter, c[1])]
    for name in sorted({c[3] for c in todo}):
        pairs, length, err, indels, flank = DATASETS[name]
        if a.quick:
            pairs = max(10, pairs // 10)
        path = os.path.join(datadir, name + ".seq")
        if not os.path.exists(path):
            print("generating %s (%d pairs)" % (name, pairs), file=sys.stderr)
            generate(path + ".tmp", pairs, length, err, indels, flank, seed=sum(map(ord, name)))
            os.rename(path + ".tmp", path)

    info = machine_info() + ["instructions: %s" % ("perf stat" if use_perf else ("time -l" if sys.platform == "darwin" else "unavailable")),
                             "tools: %s" % ", ".join(sorted(tools)), "runs per case: %d" % reps]
    print("\n".join(info), file=sys.stderr)
    tsv = open(os.path.join(outdir, "results.tsv"), "w")
    tsv.write("suite\tcase\ttool\tms\tinstructions\tmax_rss_mb\toutput_md5\trc\n")
    results = {}
    t0 = time.time()
    for i, (suite, name, args, data) in enumerate(todo):
        inp = os.path.join(datadir, data + ".seq")
        for tool, binary in sorted(tools.items()):
            targs = rust_args(args) if tool == "Rust" else args
            out = os.path.join(outdir, "out_%s.txt" % tool)
            runs = [run_once(binary, targs, inp, out, a.pin, use_perf) for _ in range(reps)]
            ms = min((r[0] for r in runs if r[0] is not None), default=None)
            ins = min((r[1] for r in runs if r[1] is not None), default=None)
            rss = max(r[2] for r in runs)
            rc = runs[0][3]
            h = md5(out)
            results[(name, tool)] = (ms, ins, rss, h, rc)
            tsv.write("%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%d\n" % (suite, name, tool, "" if ms is None else "%.1f" % ms,
                                                              "" if ins is None else ins, rss, h, rc))
            tsv.flush()
        r, c = results[(name, "Rust")], results.get((name, "C"))
        print("[%d/%d %.0fs] %s: Rust/C time %s" % (i + 1, len(todo), time.time() - t0, name,
                                                  ratio(r[0], c[0] if c else None)), file=sys.stderr)
    tsv.close()

    # Summary: ratios vs generic C (and C-native), identical output vs generic C
    md = ["# Rust vs C sweep: %s" % host, ""] + ["- " + x for x in info] + [""]
    md += ["Ratios are Rust / C (below 1: Rust better). \"Output\" compares Rust with generic C; "
           "C-native may legitimately differ on ends-free (HPCBio/wfa2lib-rs#34).", ""]
    native = "C-native" in tools
    for suite in a.suites:
        rows = [c for c in todo if c[0] == suite]
        if not rows:
            continue
        md += ["## %s" % suite, ""]
        hdr = "| Case | time vs C | instr. vs C |" + (" time vs C-native |" if native else "") + \
              " C ms | Rust ms | C MB | Rust MB | Output |"
        md += [hdr, "|" + "---|" * (hdr.count("|") - 1)]
        for _, name, _, _ in rows:
            r, c = results[(name, "Rust")], results[(name, "C")]
            if r[4] != 0 and c[4] != 0:
                same = "both rejected" if r[4] == c[4] else "rc %d/%d" % (c[4], r[4])
            else:
                same = "same" if r[3] == c[3] else "**DIFF**"
            line = "| %s | %s | %s |" % (name, ratio(r[0], c[0]), ratio(r[1], c[1]))
            if native:
                line += " %s |" % ratio(r[0], results[(name, "C-native")][0])
            line += " %s | %s | %.0f | %.0f | %s |" % ("%.0f" % c[0] if c[0] else "—", "%.0f" % r[0] if r[0] else "—",
                                                       c[2], r[2], same)
            md.append(line)
        md.append("")
    with open(os.path.join(outdir, "summary.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    print("wrote %s/results.tsv and summary.md (%.0f s)" % (outdir, time.time() - t0), file=sys.stderr)


if __name__ == "__main__":
    main()
