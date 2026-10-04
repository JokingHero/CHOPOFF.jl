#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT_DIR/scripts/dev_env.sh"

printf "Using Julia: %s\n" "$JULIA_BIN"
printf "Using depot: %s\n" "$JULIA_DEPOT_PATH"

cd "$ROOT_DIR/test/verification/julia_test"
"$JULIA_BIN" --project=../../.. check_encoding.jl
"$JULIA_BIN" --project=../../.. run_check.jl
"$JULIA_BIN" --project=../../.. run_search_check.jl
"$JULIA_BIN" --project=../../.. run_search_check_avx512.jl

cd "$ROOT_DIR"
echo "Sassy verification init completed."
