```@meta
CollapsedDocStrings = true
```

# Sassy search

## What this mode is

`search_sassy` runs off-target search directly on a genome sequence (FASTA/2bit)
using the SASSY reimplementation in CHOPOFF.

Compared to index-based methods (`prefixHashDB`, `linearDB`, etc.), this mode:

- does not require building a CHOPOFF search database first,
- uses SIMD bit-parallel candidate generation,
- performs full DP traceback for accepted candidates.

## High-level pipeline

For each guide and strand:

1. Encode guide pattern and run `search_sassy_impl` to collect candidate end positions.
2. Expand each candidate in a `-k:k` shift window.
3. Apply strict PAM checks (default behavior).
4. Build local alignment windows.
5. Run traceback backend (`:custom` by default, `:align` optional).
6. Emit deduplicated results in CHOPOFF output format.

Primary implementation:

- `src/sassy/core.jl`
- `src/sassy/interface.jl`

## Traceback backends

`search_sassy` supports:

- `traceback_backend = :custom` (default):
  full DP traceback implementation in `src/sassy/interface.jl`.
- `traceback_backend = :align`:
  delegates traceback to legacy CHOPOFF `align`.

Notes:

- Both backends return full alignment strings.
- Distances and core hit identity are expected to match.
- Exact gap placement can differ when multiple optimal paths exist.

## Strict PAM behavior

Default is `strict_pam = true`.

In strict mode:

- SASSY candidate generation is guide-centric.
- PAM is validated at candidate coordinates before expensive traceback.
- Reference ambiguity is limited by `motif.ambig_max` within the candidate
  motif window (guide + PAM), before traceback.
- Ambiguous guide bases are still supported through IUPAC pattern matching and
  are not counted against `motif.ambig_max`.

This is the default parity mode used for CHOPOFF SASSY verification.

## API

```@docs
search_sassy
```

## Command-line usage

Search with SASSY directly:

```bash
julia --threads 8 --project=. src/CHOPOFF.jl search sassy \
  --guides test/sample_data/guides.txt \
  --genome test/sample_data/genome/semirandom.fa \
  --output /tmp/sassy_hits.csv \
  --distance 3
```

Notes:

- `--force_safe_minima` disables BMI2/PEXT minima and forces the safe fallback.

## Current parity contract

Current acceptance parity with `linearDB` is enforced on:

- `guide`
- `distance`
- `chromosome`
- `start`
- `strand`

Alignment path-string identity is diagnostic-only and not a hard failure when
multiple optimal alignments exist.

## Verification and benchmarks

Useful checks:

```bash
julia --project=. -e 'include("test/src/test_sassy_correctness.jl")'
julia --project=. -e 'include("test/src/test_sassy_traceback_parity.jl")'
julia --project=. -e 'ENV["CHOPOFF_VERIFY_TRACEBACK_BACKENDS"]="1"; include("test/src/verify_sassy_core.jl")'
```

Useful benchmarks:

```bash
julia --project=. scripts/benchmark_sassy_minima_backend.jl
julia --project=. scripts/benchmark_sassy_traceback.jl
julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl
```


## Human Profiling

Human-scale SASSY profiling is kept local under `test/local_human/`. Profiling
packages are installed outside this project in:

```bash
/home/rstudio/livemount/kornel_dev/temp_upload/profiletools
```

Installed tools: `PProf`, `StatProfilerHTML`, `BenchmarkTools`, `TimerOutputs`.
They are loaded by `test/local_human/profile_human_sassy.jl`; CHOPOFF's
`Project.toml` and `Manifest.toml` are not modified.

Typical CPU profile:

```bash
CHOPOFF_PROFILE_MODE=cpu \
CHOPOFF_PROFILE_USE_AVX512=1 \
JULIA_NUM_THREADS=8 \
JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=. \
  test/local_human/profile_human_sassy.jl
```

Useful modes:

- `baseline`: timed SASSY run only.
- `cpu`: stdlib `Profile` plus `PProf` and `StatProfilerHTML` output.
- `allocs`: allocation profile via `Profile.Allocs` and `PProf.Allocs`.
- `scaling`: launches baseline runs with multiple `JULIA_NUM_THREADS` values.

Key outputs:

- `profiles/cpu.pb.gz`
- `profiles/cpu_flat.txt`
- `profiles/cpu_top.txt`
- `profiles/statprof/index.html`
- `profiles/allocs.pb.gz`

Current finding from the human profile: active CPU time is dominated by
`search_sassy_impl`, especially text block encoding (`encode_block_avx2!`) and
the Myers/minima loop. Traceback is not the main bottleneck after filtering
ambiguous reference windows. The current implementation still scans the genome
per guide/strand, so batching/pattern-tiling should be evaluated before small
traceback optimizations.
