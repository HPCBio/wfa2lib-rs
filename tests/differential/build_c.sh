#!/bin/bash
# Build C WFA2-lib's align_benchmark at the pinned commit (C_COMMIT), or at
# the commit given as $1, into target/differential/c-<commit>/.
# With C_NATIVE=1, keep C's -march=native (its AVX2/AVX-512 kernels) and
# build into c-<commit>-native-<host>/ instead (per host: it targets the
# build machine's CPU); for benchmarks (bench/sweep.py), not the harness.
#
# Source: $WFA2_LIB_REPO (default: ./WFA2-lib if it is a git checkout,
# else https://github.com/Marco-Sola-Lab/WFA2-lib).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
COMMIT="${1:-$(cat "$HERE/C_COMMIT")}"
if [ -z "${WFA2_LIB_REPO:-}" ]; then
  if [ -d "$ROOT/WFA2-lib/.git" ]; then
    WFA2_LIB_REPO="$ROOT/WFA2-lib"
  else
    WFA2_LIB_REPO="https://github.com/Marco-Sola-Lab/WFA2-lib"
  fi
fi
DEST="$ROOT/target/differential/c-$COMMIT${C_NATIVE:+-native-$(hostname -s)}"
BIN="$DEST/bin/align_benchmark"

# Bump when the build recipe below changes, so cached builds are redone.
RECIPE=2
STAMP="$DEST/.build_c-recipe"

if [ -x "$BIN" ] && [ "$(cat "$STAMP" 2>/dev/null)" = "$RECIPE" ]; then
  echo "$BIN"
  exit 0
fi

# Build in a private directory and move it into place when complete, so
# concurrent runs (e.g. on several hosts sharing a checkout) don't clobber
# each other: DEST only ever holds a finished build. A DEST without the
# binary is left over from an interrupted older run.
# A DEST from an older recipe is rebuilt too.
[ -d "$DEST" ] && { [ ! -x "$BIN" ] || [ "$(cat "$STAMP" 2>/dev/null)" != "$RECIPE" ]; } && rm -rf "$DEST"
FINAL="$DEST"
DEST="$FINAL.tmp.$(hostname -s).$$"
trap 'rm -rf "$DEST"' EXIT
mkdir -p "$(dirname "$DEST")"
git clone --quiet "$WFA2_LIB_REPO" "$DEST" >&2
git -C "$DEST" checkout --quiet "$COMMIT"

# Drop -march=native so C uses its scalar kernels on every machine. With
# AVX2, C's ends-free extend checks termination in a different diagonal
# order and can return a different, equally optimal alignment; the Rust
# port follows the scalar order. See HPCBio/wfa2lib-rs#34.
# Apple clang has no OpenMP: also drop -fopenmp (unused here).
drop_flag() { # pattern
  grep -rl --include=Makefile -e "$1" "$DEST" | while read -r mk; do
    sed "s/ *$1//g" "$mk" > "$mk.tmp" && mv "$mk.tmp" "$mk"
  done
}
if [ -z "${C_NATIVE:-}" ]; then
  drop_flag '-march=native'
fi
if [ "$(uname)" = Darwin ]; then
  drop_flag '-fopenmp'
  # C's profiler reads the clock on macOS with three Mach calls
  # (host_get_clock_service, clock_get_time, mach_port_deallocate), about
  # 1.3 us per alignment inside the timed interval: most of the measured
  # time for short pairs. Use clock_gettime, as on Linux (macOS 10.12+).
  # Benchmark timing only; alignments are unchanged.
  timer="$DEST/system/profiler_timer.c"
  grep -q '^#ifdef __MACH__ // OS X does not have clock_gettime' "$timer" || {
    echo "build_c.sh: macOS timer patch does not apply to $timer" >&2; exit 1; }
  sed 's|^#ifdef __MACH__ // OS X does not have clock_gettime.*|#if 0 // build_c.sh: use clock_gettime on macOS too|' \
    "$timer" > "$timer.tmp" && mv "$timer.tmp" "$timer"
fi

make -C "$DEST" >"$DEST/build.log" 2>&1 || {
  cp "$DEST/build.log" "$FINAL.build-failed.log" 2>/dev/null
  echo "C build failed; see $FINAL.build-failed.log" >&2
  exit 1
}
[ -x "$DEST/bin/align_benchmark" ] || { echo "no binary after the C build" >&2; exit 1; }
echo "$RECIPE" > "$DEST/.build_c-recipe"
# Move into place unless another run finished first (then keep theirs).
# If that race is lost between the check and the move, mv either fails
# (BSD: directory not empty) or nests our copy inside theirs (GNU), which is
# removed below; either way their finished build is used.
if [ ! -e "$FINAL" ]; then
  mv "$DEST" "$FINAL" 2>/dev/null || true
fi
rm -rf "$FINAL/$(basename "$DEST")"
[ -x "$BIN" ] || { echo "no $BIN after the C build" >&2; exit 1; }
echo "$BIN"
