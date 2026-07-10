#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GENOME="${CHOPOFF_HUMAN_GENOME:-/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"
GUIDES="${CHOPOFF_HUMAN_GUIDES:-$ROOT_DIR/test/local_human/data/guides_for_tests.txt}"
DISTANCE="${CHOPOFF_HUMAN_DISTANCE:-3}"
MODE="${CHOPOFF_RUST_SASSY_MODE:-range}"
THREADS="${CHOPOFF_RUST_SASSY_THREADS:-8}"
OUT_PARENT="${CHOPOFF_RUST_SASSY_OUT_PARENT:-$ROOT_DIR/test/local_human/outputs}"
STAMP="$(date -u +%Y%m%d_%H%M%S)"
OUT_DIR="$OUT_PARENT/rust_v2_${MODE}_$STAMP"

mkdir -p "$OUT_DIR"

(
  cd "$ROOT_DIR/sassy"
  EXTRA_ARGS=()
  if [[ "$MODE" == "range" ]]; then
    EXTRA_ARGS+=(--range-only)
  fi
  cargo run --release --bin chopoff_batch_crispr -- \
    --guides "$GUIDES" \
    --distance "$DISTANCE" \
    --threads "$THREADS" \
    --output "$OUT_DIR/rust_sassy_v2.csv" \
    "${EXTRA_ARGS[@]}" \
    "$GENOME" 2> "$OUT_DIR/summary.csv"
)

printf 'output_dir,%s\n' "$OUT_DIR" | tee -a "$OUT_DIR/summary.csv"
