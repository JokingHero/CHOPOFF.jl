# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CHOPOFF.jl is a Julia package for **sensitive and fast CRISPR off-target detection**. It provides multiple alignment algorithms optimized for CRISPR guide RNA searches, supporting mismatches, bulges, ambiguous bases, and VCF files for personalized searches.

**Minimum Julia version:** 1.10.10

## Common Commands

### Development

```bash
# Run all tests
julia --project=. test/runtests.jl

# Run specific test file
julia --project=. -e 'include("test/src/utils.jl")'

# Build standalone application
./build_standalone.sh
# Skip precompilation for faster build
./build_standalone.sh --noprecompile

# Build documentation
./compute_coverage.sh
```

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
- `FMidx/FMindexes.jl` - FM-index implementation
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
| fmiDB | `build_fmiDB` | `search_fmiDB` | FM-index based |
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
