#!/bin/bash
# Build C WFA2-lib's align_benchmark at the pinned commit (C_COMMIT), or at
# the commit given as $1, into target/differential/c-<commit>/.
# With C_NATIVE=1, keep C's -march=native (its AVX2/AVX-512 kernels) and
# build into c-<commit>-native-<host>/ instead (per host: it targets the
# build machine's CPU); for benchmarks (bench/sweep.py), not the harness.
#
# Source: $WFA2_LIB_REPO (default: ./WFA2-lib if it is a git checkout,
# else https://github.com/smarco/WFA2-lib).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
COMMIT="${1:-$(cat "$HERE/C_COMMIT")}"
if [ -z "${WFA2_LIB_REPO:-}" ]; then
  if [ -d "$ROOT/WFA2-lib/.git" ]; then
    WFA2_LIB_REPO="$ROOT/WFA2-lib"
  else
    WFA2_LIB_REPO="https://github.com/smarco/WFA2-lib"
  fi
fi
DEST="$ROOT/target/differential/c-$COMMIT${C_NATIVE:+-native-$(hostname -s)}"
BIN="$DEST/bin/align_benchmark"

if [ -x "$BIN" ]; then
  echo "$BIN"
  exit 0
fi

rm -rf "$DEST"
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
fi

make -C "$DEST" >"$DEST/build.log" 2>&1 || {
  echo "C build failed; see $DEST/build.log" >&2
  exit 1
}
[ -x "$BIN" ] || { echo "no $BIN after build" >&2; exit 1; }
echo "$BIN"
