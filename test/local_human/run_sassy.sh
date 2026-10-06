#!/usr/bin/env bash
# Stock Sassy baseline on GRCh38, no CHOPOFF-specific code.
# CHOPOFF_SASSY_MODE=search: PAMless `sassy search` on the guides.
# CHOPOFF_SASSY_MODE=crispr: `sassy crispr` on guide+PAM; k excludes the PAM.
# Writes timings.csv (config,distance,run,wall_s,rows) and one TSV per config/distance.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SASSY="${CHOPOFF_SASSY_BIN:-$ROOT_DIR/sassy/target/release/sassy}"
MODE="${CHOPOFF_SASSY_MODE:-search}"
GENOME="${CHOPOFF_SASSY_GENOME:-/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"
GUIDES="${CHOPOFF_SASSY_GUIDES:-$ROOT_DIR/test/local_human/data/guides_for_tests.txt}"
PAM="${CHOPOFF_SASSY_PAM:-NGG}"
OUT="${CHOPOFF_SASSY_OUT:-$ROOT_DIR/test/local_human/outputs/sassy_${MODE}_$(date -u +%Y%m%d)}"
THREADS="${CHOPOFF_SASSY_THREADS:-24}"
DISTANCES="${CHOPOFF_SASSY_DISTANCES:-0 1 2 3}"
RUNS="${CHOPOFF_SASSY_RUNS:-1}"
# name:extra-args, "_" stands for a space; max-n-frac 0 matches CHOPOFF ambig_max=0
case "$MODE" in
  search) DEFAULT_CONFIGS="v1_iupac: v1_dna:-a_dna v2_iupac:--v2" ;;
  crispr) DEFAULT_CONFIGS="crispr:" ;;
  *) echo "CHOPOFF_SASSY_MODE must be search or crispr" >&2; exit 1 ;;
esac
CONFIGS="${CHOPOFF_SASSY_CONFIGS:-$DEFAULT_CONFIGS}"

mkdir -p "$OUT"
[[ -f "$OUT/timings.csv" ]] || echo "config,distance,run,wall_s,rows" > "$OUT/timings.csv"
(cd "$ROOT_DIR/sassy" && git log -1 --format='sassy %h %cd %s') > "$OUT/version.txt"
if [[ "$MODE" == crispr ]]; then
  # sassy crispr expects the PAM at the 3' end of each guide
  sed "/^$/d; s/\$/$PAM/" "$GUIDES" > "$OUT/guides_with_pam.txt"
fi

for d in $DISTANCES; do
  for config in $CONFIGS; do
    name="${config%%:*}"
    extra="${config#*:}"
    extra="${extra//_/ }"
    for run in $(seq 1 "$RUNS"); do
      tsv="$OUT/${name}_d${d}.tsv"
      log="$OUT/${name}_d${d}.log"
      start=$EPOCHREALTIME
      # shellcheck disable=SC2086
      if [[ "$MODE" == search ]]; then
        "$SASSY" search -l "$GUIDES" -k "$d" -j "$THREADS" --max-n-frac 0 $extra \
          "$GENOME" > "$tsv" 2> "$log"
      else
        # crispr prints its log to stdout, so the TSV goes through -o
        "$SASSY" crispr -g "$OUT/guides_with_pam.txt" -k "$d" -j "$THREADS" \
          --max-n-frac 0 -o "$tsv" $extra "$GENOME" > "$log" 2>&1
      fi
      wall=$(awk -v a="$start" -v b="$EPOCHREALTIME" 'BEGIN{printf "%.3f", b-a}')
      rows=$(($(wc -l < "$tsv") - 1))
      echo "$name,$d,$run,$wall,$rows" | tee -a "$OUT/timings.csv"
    done
  done
done
