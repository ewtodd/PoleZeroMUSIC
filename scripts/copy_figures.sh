#!/usr/bin/env bash
# Copy the PDFs and tables written by figures into the paper's figures directory.
#   bash scripts/copy_figures.sh [source-dir]
# The source defaults to doc/figs; set DEST to override the paper directory.
set -euo pipefail

cd "$(dirname "$0")/.."

SOURCE=${1:-doc/figs}
DEST=${DEST:-/home/e-work/Papers/Writing/WIP-PZListMode-NIMA/figures}

if ! compgen -G "$SOURCE/*.pdf" > /dev/null; then
  echo "no PDFs in $SOURCE" >&2
  exit 1
fi

mkdir -p "$DEST"
shopt -s nullglob
for file in "$SOURCE"/*.pdf "$SOURCE"/*.tex; do
  cp -v "$file" "$DEST/"
done
