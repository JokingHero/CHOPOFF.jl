#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JULIA_BIN="${JULIA_BIN:-/home/ai/.julia/juliaup/julia-1.10.10+0.x64.linux.gnu/bin/julia}"

# Prefer writable depot first to avoid lockfile issues in restricted environments.
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/tmp/julia-depot:${HOME}/.julia}"

printf "Using Julia: %s\n" "$JULIA_BIN"
printf "Using depot: %s\n" "$JULIA_DEPOT_PATH"

cd "$ROOT_DIR/test/verification/julia_test"
"$JULIA_BIN" --project=../../.. check_encoding.jl
"$JULIA_BIN" --project=../../.. run_check.jl
"$JULIA_BIN" --project=../../.. run_search_check.jl
"$JULIA_BIN" --project=../../.. run_search_check_avx512.jl

cd "$ROOT_DIR"
echo "Sassy verification init completed."
