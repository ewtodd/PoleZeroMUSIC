#!/usr/bin/env bash
# The reaction events MUSIC's 37Cl analysis tagged and its pure-beam sample,
# with and without the pole-zero correction, for the paper's reaction figure:
# output/reaction_events/reaction_events.root (see scripts/reaction_events.cpp).
#
# Reads MUSIC's own outputs only: the production strip-sum-scatter cache (with
# run numbers, cache v46) and the no-correction events of
# scripts/nocorrection.sh, so run that first. Builds libmusic.a from the local
# ~/MUSIC tooling with MUSIC's own Makefile into a temporary directory and
# links the program against it.
# Run from the project root:  bash scripts/reaction_events.sh
set -euo pipefail
cd "$(dirname "$0")/.."
MUSIC=${MUSIC:-$HOME/MUSIC}
CACHE=${CACHE:-$MUSIC/analysis/37Cl/root_files/StripSumScatter_cache.root}
NOCORRECTION=${NOCORRECTION:-$PWD/output/nocorrection}
OUT=$PWD/output/reaction_events
JOBS=${JOBS:-8}
mkdir -p "$OUT"
WORK=$(mktemp -d "$OUT/work.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

PROJECT=$PWD
# MUSIC's shell takes its dataset paths from where it starts, so start it there.
export MUSIC CACHE NOCORRECTION OUT WORK JOBS PROJECT
cd "$MUSIC"
nix develop -c bash -c '
  set -euo pipefail
  make -s -C "$MUSIC/tooling" MUSIC_DATASET=37Cl WITH_GPU=0 \
    BUILD_DIR="$WORK/build" "$WORK/build/libmusic.a" -j"$JOBS"
  g++ -O2 -std=c++17 $(root-config --cflags) -I"$MUSIC/tooling/include" \
    "-DR__ADD_INCLUDE_PATH(...)=" "$PROJECT/scripts/reaction_events.cpp" \
    "$WORK/build/libmusic.a" -Wl,--gc-sections $(root-config --glibs) \
    -lSpectrum -lMinuit -l:libanalysis-utils.so -ldl -lpthread \
    -o "$WORK/reaction_events"
  "$WORK/reaction_events" "$CACHE" "$NOCORRECTION" \
    "$OUT/reaction_events.root" 2>&1 | grep -v "^Info in" |
    tee "$OUT/reaction_events.log"'
