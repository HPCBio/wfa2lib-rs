#!/bin/bash
# Differential test: run every case from cases.sh on Rust's and C's
# align_benchmark and compare the outputs byte for byte.
#
#   tests/differential/run.sh [SUITE...]     (default: all suites)
#
# Environment:
#   JOBS     parallel cases (default: half the CPUs; high-memory runs on the
#            10 kb utest pairs need ~2 GB each)
#   TIMEOUT  seconds per run (default: 300)
#   C_BIN    C align_benchmark (default: built by build_c.sh at C_COMMIT)
#   RUST_BIN Rust align_benchmark (default: cargo build --release)
#
# Writes target/differential/results.tsv. Exits non-zero if a case differs
# and isn't listed in known_diffs.tsv.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
WORK="$ROOT/target/differential"
ncpu=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
JOBS="${JOBS:-$(( ncpu > 1 ? ncpu / 2 : 1 ))}"
export TIMEOUT="${TIMEOUT:-300}"

if [ -z "${RUST_BIN:-}" ]; then
  (cd "$ROOT" && cargo build --release --features cli --quiet) || exit 1
  RUST_BIN="$ROOT/target/release/align_benchmark"
fi
if [ -z "${C_BIN:-}" ]; then
  C_BIN="$("$HERE/build_c.sh")" || exit 1
fi
export RUST_BIN C_BIN

# Datasets: the stock C utest set and 300 x 2 kb pairs at 8% error with
# long indels (C generate_dataset -n 300 -l 2000 -e 0.08 --indels 3,40;
# generate_dataset seeds from the clock, so the file is committed).
export DATA="$WORK/data"
mkdir -p "$DATA"
C_TREE="$(cd "$(dirname "$C_BIN")/.." && pwd)"
cp "$C_TREE/tests/wfa.utest.seq" "$DATA/utest.seq"
gunzip -c "$HERE/data/long.seq.gz" > "$DATA/long.seq"

CASES="$WORK/cases.tsv"
if [ $# -gt 0 ]; then
  pattern="^($(IFS='|'; echo "$*"))	"
  bash "$HERE/cases.sh" | grep -E "$pattern" > "$CASES"
else
  bash "$HERE/cases.sh" > "$CASES"
fi
total=$(wc -l < "$CASES" | tr -d ' ')
echo "Running $total cases ($JOBS jobs); C: $C_BIN" >&2

OUT="$WORK/out"
rm -rf "$OUT"
mkdir -p "$OUT"
export OUT

# Run one case (line number $1 of $CASES): print the result row.
run_case() {
  local n=$1 line suite ds args rr rc
  line=$(sed -n "${n}p" "$CASES")
  suite=$(printf '%s' "$line" | cut -f1)
  ds=$(printf '%s' "$line" | cut -f2)
  args=$(printf '%s' "$line" | cut -f3)
  local r="$OUT/$n.rust" c="$OUT/$n.c"
  # shellcheck disable=SC2086
  perl -e 'alarm shift; exec @ARGV' "$TIMEOUT" "$RUST_BIN" -i "$DATA/$ds.seq" -o "$r" $args >/dev/null 2>&1; rr=$?
  # shellcheck disable=SC2086
  perl -e 'alarm shift; exec @ARGV' "$TIMEOUT" "$C_BIN" -i "$DATA/$ds.seq" -o "$c" $args >/dev/null 2>&1; rc=$?
  local status
  # perl's alarm kills with SIGALRM: exit 142.
  if [ $rr -eq 142 ] || [ $rc -eq 142 ]; then
    status="timeout (rust=$rr c=$rc)"
  elif [ $rr -ne 0 ] || [ $rc -ne 0 ]; then
    if [ $rr -eq $rc ]; then status="same-error"; else status="differ (exit rust=$rr c=$rc)"; fi
  elif cmp -s "$r" "$c"; then
    status="identical"
  else
    status="differ ($(paste "$r" "$c" | awk -F'\t' '$1!=$3{s++} $1==$3&&$2!=$4{g++}
      END{printf "%d scores, %d cigars of %d", s, g, NR}'))"
  fi
  rm -f "$r" "$c"
  printf '%s\t%s\t%s\t%s\n' "$suite" "$ds" "$args" "$status"
}
export -f run_case
export CASES

RESULTS="$WORK/results.tsv"
seq 1 "$total" | xargs -P "$JOBS" -I{} bash -c 'run_case {}' > "$RESULTS.tmp"
sort "$RESULTS.tmp" > "$RESULTS"
rm -f "$RESULTS.tmp"

# Summary per suite.
awk -F'\t' '{n[$1]++; s=$4; sub(/ .*/, "", s); c[$1" "s]++}
  END{for (k in n) printf "%-9s %4d cases: %4d identical, %3d same-error, %3d differ, %d timeout\n",
      k, n[k], c[k" identical"], c[k" same-error"], c[k" differ"], c[k" timeout"]}' "$RESULTS" | sort

# Compare differences with the known list (suite, dataset, args, reason).
KNOWN="$HERE/known_diffs.tsv"
grep -v '^#' "$KNOWN" | cut -f1-3 | sort > "$WORK/known.keys"
awk -F'\t' '$4 !~ /^(identical|same-error)/ {print $1"\t"$2"\t"$3}' "$RESULTS" | sort > "$WORK/differ.keys"
cut -f1-3 "$CASES" | sort > "$WORK/run.keys"
unexpected=$(comm -23 "$WORK/differ.keys" "$WORK/known.keys")
# Known differences that were run and now match.
fixed=$(comm -13 "$WORK/differ.keys" "$WORK/known.keys" | comm -12 - "$WORK/run.keys")

if [ -n "$fixed" ]; then
  echo
  echo "Known differences that now match C (remove from known_diffs.tsv):"
  echo "$fixed"
fi
if [ -n "$unexpected" ]; then
  echo
  echo "UNEXPECTED differences from C:"
  printf '%s\n' "$unexpected" | awk -F'\t' 'NR==FNR {u[$1"\t"$2"\t"$3]; next}
    ($1"\t"$2"\t"$3) in u' - "$RESULTS"
  exit 1
fi
echo
echo "OK: no unexpected differences ($(wc -l < "$WORK/differ.keys" | tr -d ' ') known). Results: $RESULTS"
