# Local Human Genome SASSY Benchmark

Manual human-genome benchmark for SASSY. It is not included by `test/runtests.jl`.

Default genome:

```text
/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

Default guides:

```text
test/local_human/data/guides_for_tests.txt
```

Run:

```bash
JULIA_NUM_THREADS=8 \
JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=. \
test/local_human/run_human_sassy.jl
```

The runner builds/reuses prefixHashDB by default. To run SASSY only:

```bash
CHOPOFF_HUMAN_COMPARE_PREFIX=0 JULIA_NUM_THREADS=8 \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=. \
test/local_human/run_human_sassy.jl
```

Generated `data/`, `indexes/`, and `outputs/` are ignored.

## Rust SASSY v2 Batch Baseline

This local-only baseline calls Rust SASSY v2 `encode_patterns` +
`search_encoded_patterns` directly. It batches all guide+`NGG` patterns, then
applies strict PAM and reference ambiguity filtering. It is for speed ceiling
measurement, not exact CHOPOFF output parity: it does not run CHOPOFF traceback
or emit CHOPOFF's guide/PAM-side normalized coordinates.

Run after Cargo/Rust is available on `PATH`:

```bash
test/local_human/run_rust_sassy_v2.sh
```

Outputs are written to `test/local_human/outputs/rust_v2_<timestamp>/`:

- `summary.csv`: records, guide count, raw/accepted match counts, elapsed seconds.
- `rust_sassy_v2.csv`: accepted candidate rows from the Rust v2 batch search.

Compare `elapsed_s` with the Julia human runner and profiler outputs to decide
whether porting v2 pattern tiling is worth doing before micro-optimizing the
current per-guide Julia scan.

