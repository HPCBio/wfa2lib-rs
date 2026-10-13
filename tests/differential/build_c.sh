#!/bin/bash
# Build C WFA2-lib's align_benchmark at the pinned commit (C_COMMIT), or at
# the commit given as $1, into target/differential/c-<commit>/.
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
DEST="$ROOT/target/differential/c-$COMMIT"
BIN="$DEST/bin/align_benchmark"

if [ -x "$BIN" ]; then
  echo "$BIN"
  exit 0
fi

rm -rf "$DEST"
mkdir -p "$(dirname "$DEST")"
git clone --quiet "$WFA2_LIB_REPO" "$DEST" >&2
git -C "$DEST" checkout --quiet "$COMMIT"

# Apple clang has no OpenMP: drop -fopenmp (the parallel code is unused here).
if [ "$(uname)" = Darwin ]; then
  grep -rl --include=Makefile -e '-fopenmp' "$DEST" | while read -r mk; do
    sed -i '' 's/ *-fopenmp//g' "$mk"
  done
fi

make -C "$DEST" >"$DEST/build.log" 2>&1 || {
  echo "C build failed; see $DEST/build.log" >&2
  exit 1
}
[ -x "$BIN" ] || { echo "no $BIN after build" >&2; exit 1; }
echo "$BIN"
