#!/usr/bin/env python3
"""Summarize bench/sweep.py runs from several machines as Markdown tables.

    python3 bench/report.py RUN_DIR [RUN_DIR ...] > tables.md

Each RUN_DIR is a sweep output directory (results.tsv and summary.md). Prints
the machines, Rust/C time by group of cases and per algorithm, peak memory
(Linux runs only: macOS reports RSS differently) and an overall summary, one
column per run. docs/performance.md is built from this output.
"""
import collections
import csv
import os
import re
import statistics as st
import sys


def meta(d):
    m = {}
    for line in open(os.path.join(d, "summary.md")):
        mm = re.match(r"- (\w[\w ]*): (.*)", line)
        if mm:
            m[mm.group(1)] = mm.group(2).strip()
    cpu = re.sub(r"\(R\)|\(TM\)| CPU| @ [\d.]+GHz| \d+-Core Processor", "", m.get("cpu", "?"))
    m["short"] = cpu.strip()
    m["linux"] = m.get("os", "").startswith("Linux")
    return m


def cc_name(cc):
    g = re.search(r"gcc ([\d.]+)", cc)
    if g:
        return "GCC " + g.group(1)
    a = re.search(r"(Apple )?clang version ([\d.]+)", cc)
    return ("%sclang %s" % (a.group(1) or "", a.group(2))) if a else cc


def load(d):
    rows = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(d, "results.tsv")), delimiter="\t"):
        rows[(r["suite"], r["case"])][r["tool"]] = r
    return rows


def big(w):
    return w[0] != "s150_2"


GROUPS = [  # (label, predicate on suite and case words)
    ("Score only, high memory", lambda s, w: s == "core" and big(w) and w[-1] == "score"),
    ("CIGAR, high memory", lambda s, w: s == "core" and big(w) and w[-1] == "cigar"),
    ("CIGAR, med/low memory", lambda s, w: s == "core" and big(w) and w[-1] in ("med", "low")),
    ("CIGAR, BiWFA (ultralow memory)", lambda s, w: s == "core" and big(w) and w[-1] == "ultralow"),
    ("Ends-free / extension, high memory", lambda s, w: s == "span" and w[-1] in ("score", "cigar")),
    ("Ends-free / extension, med memory", lambda s, w: s == "span" and w[-1] == "med"),
    ("Heuristics (adaptive, X/Z-drop, banded)", lambda s, w: s == "heur"),
    ("150 bp: high memory and BiWFA", lambda s, w: s == "core" and not big(w) and w[-1] in ("score", "cigar", "ultralow")),
    ("150 bp: med/low memory", lambda s, w: s == "core" and not big(w) and w[-1] in ("med", "low")),
]

MEM_GROUPS = [
    ("CIGAR, high memory", lambda s, w: s == "core" and big(w) and w[-1] == "cigar"),
    ("CIGAR, med/low memory", lambda s, w: s == "core" and big(w) and w[-1] in ("med", "low")),
    ("CIGAR, BiWFA (ultralow memory)", lambda s, w: s == "core" and big(w) and w[-1] == "ultralow"),
    ("Ends-free / extension, CIGAR", lambda s, w: s == "span" and w[-1] == "cigar"),
]

ALGS = [("edit", "Edit"), ("indel", "Indel"), ("linear", "Gap-linear"),
        ("affine", "Gap-affine"), ("affine2p", "Gap-affine 2-piece")]
MODES = [("score", "score only"), ("cigar", "CIGAR"), ("med", "CIGAR, med memory"), ("ultralow", "CIGAR, BiWFA")]


def ratio(t, tool, key):
    return float(t["Rust"][key]) / float(t[tool][key])


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out + [""]


def main(dirs):
    runs = [(meta(d), load(d)) for d in dirs]
    names = [m["short"] for m, _ in runs]
    out = ["### Machines", ""]
    out += table(["Machine", "OS", "Rust commit", "rustc", "C compiler"],
                 [[m["short"], m.get("os", "?").split(" ")[0], "`%s`" % m.get("rust commit", "?"),
                   m.get("rustc", "?").split(" (")[0].replace("rustc ", ""), cc_name(m.get("cc", "?"))]
                  for m, _ in runs])

    out += ["### Time by group", "",
            "Rust time / C time, median (range) over the cases in each group; below 1 means Rust is faster.", ""]
    rows = []
    for label, pred in GROUPS:
        cells, n = [], 0
        for _, rr in runs:
            v = [ratio(t, "C", "ms") for (s, c), t in rr.items() if pred(s, c.split())]
            n = max(n, len(v))
            cells.append("%.2f (%.2f–%.2f)" % (st.median(v), min(v), max(v)) if v else "–")
        rows.append([label, str(n)] + cells)
    out += table(["Cases", "n"] + names, rows)

    out += ["### Time per algorithm", "",
            "Rust time / C time on `m1k_5` (20,000 pairs of 1 kb, 5% error).", ""]
    rows = []
    for alg, aname in ALGS:
        for mode, mlabel in MODES:
            cells = []
            for _, rr in runs:
                t = rr.get(("core", "m1k_5 %s %s" % (alg, mode)))
                cells.append("%.2f" % ratio(t, "C", "ms") if t else "–")
            rows.append([aname, mlabel] + cells)
    out += table(["Algorithm", "Mode"] + names, rows)

    linux = [(m, rr) for m, rr in runs if m["linux"]]
    if linux:
        out += ["### Peak memory", "",
                "Rust peak RSS / C peak RSS, median (max) over the cases; Linux runs only "
                "(macOS reports RSS differently).", ""]
        rows = []
        for label, pred in MEM_GROUPS:
            cells = []
            for _, rr in linux:
                v = [ratio(t, "C", "max_rss_mb") for (s, c), t in rr.items() if pred(s, c.split())]
                cells.append("%.2f (%.2f)" % (st.median(v), max(v)) if v else "–")
            rows.append([label] + cells)
        out += table(["Cases"] + [m["short"] for m, _ in linux], rows)

    out += ["### Overall", ""]
    lines = collections.OrderedDict()
    for _, rr in runs:
        ts = [ratio(t, "C", "ms") for t in rr.values()]
        tn = [ratio(t, "C-native", "ms") for t in rr.values() if "C-native" in t]
        same = sum(t["Rust"]["output_md5"] == t["C"]["output_md5"] for t in rr.values())
        lines.setdefault("Cases faster than C", []).append("%d / %d" % (sum(x < 1 for x in ts), len(ts)))
        lines.setdefault("Median Rust / C", []).append("%.2f" % st.median(ts))
        lines.setdefault("Median Rust / C built with `-march=native`", []).append(
            "%.2f" % st.median(tn) if tn else "–")
        lines.setdefault("Output identical to C", []).append("%d / %d" % (same, len(rr)))
    out += table([""] + names, [[k] + v for k, v in lines.items()])
    print("\n".join(out).rstrip())


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
