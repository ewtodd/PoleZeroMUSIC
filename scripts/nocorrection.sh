#!/usr/bin/env bash
# MUSIC's 37Cl pipeline and reaction search with the pole-zero correction off,
# for the late epoch (runs 97-137), into output/nocorrection: built events,
# calibration, summaries, plots, and the strip-sum-scatter cache with MUSIC's
# own tags on uncorrected data (see scripts/nocorrection.cpp). The production
# results in ~/MUSIC are not touched: MUSIC_RESULTS_DIR points here instead,
# and the one file the pipeline writes into the dataset directory, its log,
# is moved here afterwards and the production log put back.
#
# About as long as a production pipeline run and about 25 GB of disk. Then run
# scripts/reaction_events.sh for the reaction figures.
# Run from the project root:  bash scripts/nocorrection.sh
set -euo pipefail
cd "$(dirname "$0")/.."
MUSIC=${MUSIC:-$HOME/MUSIC}
OUT=$PWD/output/nocorrection
JOBS=${JOBS:-8}
mkdir -p "$OUT"
WORK=$(mktemp -d "$OUT/work.XXXXXX")
LOG=$MUSIC/analysis/37Cl/pipeline_fused.log
if [ -f "$LOG" ]; then cp -p "$LOG" "$WORK/production_pipeline_fused.log"; fi
restore() {
  if [ -f "$LOG" ]; then mv "$LOG" "$OUT/pipeline_fused.log"; fi
  if [ -f "$WORK/production_pipeline_fused.log" ]; then
    mv "$WORK/production_pipeline_fused.log" "$LOG"
  fi
  rm -rf "$WORK"
}
trap restore EXIT

PROJECT=$PWD
# MUSIC's shell takes its dataset paths from where it starts, so start it there.
export MUSIC OUT WORK JOBS PROJECT
cd "$MUSIC"
nix develop -c bash -c '
  set -euo pipefail
  make -s -C "$MUSIC/tooling" MUSIC_DATASET=37Cl WITH_GPU=0 \
    BUILD_DIR="$WORK/build" "$WORK/build/libmusic.a" -j"$JOBS"
  g++ -O2 -std=c++17 $(root-config --cflags) -I"$MUSIC/tooling/include" \
    "-DR__ADD_INCLUDE_PATH(...)=" "$PROJECT/scripts/nocorrection.cpp" \
    "$WORK/build/libmusic.a" -Wl,--gc-sections $(root-config --glibs) \
    -lSpectrum -lMinuit -l:libanalysis-utils.so -ldl -lpthread \
    -o "$WORK/nocorrection"
  MUSIC_RESULTS_DIR="$OUT" "$WORK/nocorrection" 2>&1 | grep -v "^Info in" |
    tee "$OUT/nocorrection.log"'
