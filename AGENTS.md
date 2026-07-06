# AGENTS.md

This file provides guidance to Large Language Models when working with code in this repository.

## Project Overview

CHOPOFF.jl is a Julia package for **sensitive and fast CRISPR off-target detection**. It provides multiple alignment algorithms optimized for CRISPR guide RNA searches, supporting mismatches, bulges, ambiguous bases, and VCF files for personalized searches.

**Minimum Julia version:** 1.10.10

## Local Environment

This workspace uses a full Julia install outside the repo:

```bash
export JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot:
JULIA=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia
```

Do not use `build/bin/julia` for development. It belongs to the standalone app output and is missing files needed by `Pkg`.

If Git reports dubious ownership, use a process-local override or configure this repo as safe:

```bash
git config --global --add safe.directory /home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl
```

## Common Commands

### Development

```bash
# Run all tests (test runner assumes cwd=test)
(cd test && JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
  /home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=.. runtests.jl)

# Run specific test file
JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
  /home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=. -e 'include("test/src/utils.jl")'

# Build standalone application
./build_standalone.sh
# Skip precompilation for faster build
./build_standalone.sh --noprecompile

# Build documentation
./compute_coverage.sh
```

### Speed Benchmark (Sassy vs prefixHashDB)

```bash
# Runs both thread modes:
#   1) single-thread (JULIA_NUM_THREADS=1)
#   2) current environment thread count
julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl

# Optional controls
CHOPOFF_BENCH_RUNS=11 JULIA_NUM_THREADS=8 julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl
CHOPOFF_BENCH_MODE=single julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl
CHOPOFF_BENCH_OUT=/tmp/chopoff_speed_summary.csv julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl
```

Environment variables used by the benchmark:
- `CHOPOFF_BENCH_RUNS` (default `7`) - timed repetitions per algorithm
- `CHOPOFF_BENCH_DISTANCE` (default `3`) - search distance
- `CHOPOFF_BENCH_MODE` (`both`, `single`, `env`; default `both`)
- `CHOPOFF_BENCH_OUT` - path to write combined summary CSV
- `CHOPOFF_BENCH_KEEP_TMP` (`0`/`1`) - keep temp benchmark artifacts

### Standalone Binary Usage (after build)

```bash
# Set thread count (important for performance)
export JULIA_NUM_THREADS=6  # use 1 for build, more for search

# Build database
CHOPOFF build --name Cas9_hg38 --genome genome.fa -o output_dir --distance 3 --motif Cas9 prefixHashDB

# Search off-targets
CHOPOFF search --database output_dir --guides guides.txt --output results.csv --distance 2 prefixHashDB

# List available motifs
CHOPOFF list

# Filter overlapping off-targets
CHOPOFF filter --detail_file results.csv --output filtered.csv --distance 3
```

## Architecture

### Module Structure (`src/CHOPOFF.jl`)

Core modules are loaded in order:
- `persistence.jl` - save/load functionality
- `distance_metrics.jl` - Hamming, Levenshtein distances
- `motif.jl` - Motif definition and handling
- `db_info.jl` - database metadata (`GenomeInfo`, `DBInfo`, `Offtarget`, `Loc`)
- `motif_path_templates.jl` - predefined motif templates
- `find_offtargets.jl` - core off-target finding logic
- `db_helpers.jl` - shared database helper functions
- Database implementations (`db_*.jl` files)
- `sassy.jl` - SIMD-Myers transposed algorithm

### Database Types

Each database has `build_*DB` and `search_*DB` functions:

| Database | Build Function | Search Function | Characteristics |
|----------|---------------|-----------------|-----------------|
| linearDB | `build_linearDB` | `search_linearDB` | Simple prefix-based, early stopping supported |
| treeDB | `build_treeDB` | `search_treeDB` | Vantage Point tree with triangle inequality |
| prefixHashDB | `build_prefixHashDB` | `search_prefixHashDB` | **Default choice**, fast with hash filtering |
| motifDB | `build_motifDB` | `search_motifDB` | Q-gram filtering |
| hashDB | `build_hashDB` | `search_hashDB` | Distance=1 only, estimation |
| dictDB | `build_dictDB` | `search_dictDB` | Simple dictionary of unique guides |
| vcfDB | `build_vcfDB` | `search_vcfDB` | VCF/personalized search |
| sassy | N/A (direct search) | `search_sassy` | SIMD-optimized |

### Key Data Structures

**`Motif`** (`src/motif.jl`): Defines search pattern
- `fwd`, `rve`: LongDNA{4} sequences for forward/reverse
- `pam_loci_fwd`, `pam_loci_rve`: PAM positions
- `distance`: max edit distance for search
- `extends5`: alignment direction (true for Cas9)
- `ambig_max`: max ambiguous bases allowed

**`DBInfo`** (`src/db_info.jl`): Database metadata
- Contains `GenomeInfo` with chromosome list, position type
- Stores motif, VCF path, creation date

**`Offtarget`** (`src/db_info.jl`): Search result
- `loc`: `Loc{T,K}` with chromosome index, position, strand
- `dist`: edit distance
- `aln_guide`, `aln_ref`: aligned sequences

**`Loc{T,K}`**: Genome location with type parameters for chromosome (T) and position (K) - chosen based on genome size to minimize memory.

### Distance Metrics

- `hamming(seq1, seq2)` - mismatch count only
- `levenshtein(seq1, seq2)` - full edit distance with SIMD optimization
- `align(seq1, seq2)` - returns `Aln` struct with distance and alignments
- `isinclusive(motif)` - checks if distance metric is inclusive of gaps

### Threading

Set `JULIA_NUM_THREADS` environment variable:
- Build: 1 thread recommended (to avoid file handle limits)
- Search: Use multiple threads (6+ for large genomes)

### File Limits

Some algorithms create many files (e.g., prefix 7 = 4^7 = 16384 files). Check `ulimit -n` and increase if needed.

## Testing

Tests are auto-included from `test/src/` by `test/runtests.jl`. Sample data in `test/sample_data/`.

**Test categories:**
- `db.jl` - database algorithm tests
- `distance_metrics.jl` - alignment verification
- `utils.jl` - utility functions
- `argparse.jl` - CLI argument parsing
- `db_extends5_false.jl` - non-Cas9 motifs
- `motif_path_templates.jl` - motif templates

## CLI Commands

Standalone binary supports: `build`, `search`, `estimate`, `filter`, `summarize`, `list`

Each database type has subcommands with specific options (see `--help` for each).

## Genome Formats

- **FASTA** (`.fa`, `.fasta`, `.fna`): requires `.fai` index
- **2bit** (`.2bit`): more compact, no index needed

## Working on SASSY alghoritm

To replicate the specific Horizontal Bit-Parallelism (Transposed Myers) + SIMD Tiling architecture used in Sassy for CRISPR search, you can inspect these files in "sassy" folder with the original sassy implementation in rust.
These files contain the exact mathematical logic for the "Transposed" step, the memory layout for the text chunks, and the specific scoring logic used to extract edit distances from the bit-vectors.

1. src/bitpacking.rs
Why: This contains the mathematical core.
Function: compute_block_simd
Significance: This implements the "Transposed Myers" step. Unlike standard Myers, it updates vertical deltas (v) and horizontal deltas (h) for a SIMD vector of text blocks. This is the equation the Julia code must replicate exactly.

2. src/search.rs
Why: This contains the search engine and tiling logic.
Structs: Searcher, LaneState
Functions: search_internal, find_minima_with_overhang
Significance:
search_internal: Shows the main loop structure. It iterates through the text in chunks (SIMD Lanes), and then iterates through the Pattern characters (Horizontal parallelism).
LaneState: Shows how text is buffered into 64-byte blocks.
find_minima_with_overhang: Shows the complex logic required to convert the final bit-vectors (Vp, Vm) back into an actual integer edit distance (score) for every position in the text.

3. src/profiles/iupac.rs
Why: This contains the SIMD-optimized text profiling.
Function: encode_ref
Significance: Because Sassy iterates over the Pattern, it must pre-process the Text chunks. This file shows how a 64-byte text block is converted into bitmasks (e.g., "which bits in this u64 correspond to 'A'?"). It handles the IUPAC logic (N, R, Y, etc.) efficiently.

4. bin/crispr.rs
Why: This contains the CRISPR-specific high-level logic.
Significance: It shows how to invoke the searcher with a filter_fn. This is crucial for CRISPR because you need to check the PAM sequence (e.g., NGG) strictly after finding a candidate match. It shows how to combine the approximate search (edit distance) with the strict suffix check.
