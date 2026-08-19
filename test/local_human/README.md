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

Select a SASSY backend with `CHOPOFF_HUMAN_BACKEND`; accepted values are
`auto`, `avx512`, `avx2_pext`, and `avx2_safe`. The summary records both the
requested and resolved backend. Profiling uses the corresponding
`CHOPOFF_PROFILE_BACKEND` variable.

For deferred Zen 1/2 qualification, run this synthetic backend benchmark on
the real server first:

```bash
SASSY_BENCH_OUT=/tmp/sassy_zen_backends.csv \
JULIA_NUM_THREADS=8 \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=. \
  scripts/benchmark_sassy_minima_backend.jl
```

It checks output parity for every supported backend. On Zen 1/2, `auto` must
resolve to `avx2_safe`. Then run the human fixture with
`CHOPOFF_HUMAN_BACKEND=auto`; compare against an explicit `avx2_safe` run if
needed. Cross-compilation with `julia -C znver2` validates generated code, not
real-hardware performance.

Generated `data/`, `indexes/`, and `outputs/` are ignored.

## Generic prefixHashDB parity matrix

`benchmark_human_generic_parity.jl` qualifies the generic prefixHashScan
engine against prefixHashDB. Qualification covers full-GRCh38 `Cas9_NGA`,
`CasX`, and custom `25N_NGG` searches with 65 distributed guides, distances
0 through 4, and ambiguity limits 0 through 3. A bounded chromosome-21 tier
covers internal/PAMless motifs, the 16-base lower guide boundary, strand
subsets, FASTA, 2bit, IUPAC symbols, and small chunk boundaries. Longer custom
guides are exercised by the full-GRCh38 `25N_NGG` case; a separate 28-base
prefixHashDB case is intentionally excluded because distance-4 symbolic-path
construction is outside the practical qualification scope.

One distance-4, `ambig_max=3` prefixHashDB is reused as the comparison superset
for every lower distance and ambiguity limit. Unambiguous rows require exact
detail parity. Ambiguous differences pass only when the reference-backed
classifier assigns a documented semantic reason. Guide lengths above 28 are
outside this matrix because prefixHashDB cannot pack a distance-4 candidate
longer than 32 bases.

Run the fast synthetic smoke matrix:

```bash
CHOPOFF_GENERIC_PARITY_MODE=smoke \
CHOPOFF_GENERIC_PARITY_OUT=/tmp/chopoff_generic_parity_smoke \
JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=. \
  test/local_human/benchmark_human_generic_parity.jl
```

Launch the full qualification as a background job:

```bash
GENERIC_OUT=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl/test/local_human/outputs/generic_parity_grch38
GENERIC_LOG=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl/test/local_human/outputs/generic_parity_grch38.log

nohup setsid env \
  JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
  CHOPOFF_GENERIC_PARITY_MODE=qualification \
  CHOPOFF_GENERIC_PARITY_OUT="$GENERIC_OUT" \
  CHOPOFF_GENERIC_PARITY_THREADS=24 \
  /home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia \
  --project=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl \
  /home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl/test/local_human/benchmark_human_generic_parity.jl \
  </dev/null >"$GENERIC_LOG" 2>&1 &

PID=$!
tail --pid="$PID" -f "$GENERIC_LOG"
```

Stages `prepare`, `build`, `search`, and `compare` can be resumed independently
with `CHOPOFF_GENERIC_PARITY_STAGE`. Set
`CHOPOFF_GENERIC_PARITY_REBUILD=1` to rebuild existing indexes. The principal
outputs are `manifest.csv`, `builds.csv`, `timings.csv`,
`parity_summary.csv`, `parity_differences.csv`, and `count_parity.csv`.

## prefixHashScan ambiguity benchmark

The tuning runner measures query construction, prepared scanning, and the
public end-to-end search independently. This command sweeps Cas9 distance 3
and `ambig_max=0:3` on the 61 standard guides with eight scan threads:

```bash
CHOPOFF_TUNING_STAGE=final \
CHOPOFF_TUNING_AMBIG_MAXES=0,1,2,3 \
CHOPOFF_TUNING_THREADS=8 \
CHOPOFF_TUNING_WARMUPS=2 \
CHOPOFF_TUNING_RUNS=11 \
CHOPOFF_TUNING_ALLOCATION_RUNS=3 \
CHOPOFF_TUNING_OUT=test/local_human/outputs/prefix_hash_scan_ambiguity \
JULIA_NUM_THREADS=8 \
JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia --project=. \
  scripts/benchmark_prefix_hash_scan_tuning.jl
```

The summary reports median, mean, standard deviation, minimum, maximum,
end-to-end coefficient of variation, and ratios/deltas against `ambig_max=0`.
The ambiguity parity CSV checks that every lower-level result multiset is
contained in the next level. Runs, allocations, diagnostics, parity, and final
result CSVs are written under `CHOPOFF_TUNING_OUT`.

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

## Overnight prefixHashDB/prefixHashScan Distance Sweep

`benchmark_human_prefix_sweep.jl` compares `prefixHashDB` and
`prefixHashScan` on GRCh38 for Cas9 and Cas12a at distances 0 through 4.
It builds one distance-4 prefixHashDB per motif with one Julia thread, then
reuses those databases for every requested search distance. Existing d4/p16
symbolic path assets are reused; the runner does not generate paths.

The default search phase uses 24 threads, 61 guides per motif, one warmup, and
five timed repetitions per algorithm and distance. Timed runs alternate
algorithm order. Raw prefixHashDB detail output is the gold standard; exact
detail rows and duplicate multiplicity are compared with prefixHashScan.

Launch the default sweep:

```bash
SWEEP_OUT=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl/test/local_human/outputs/prefix_sweep_20260722
SWEEP_LOG=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl/test/local_human/outputs/prefix_sweep_20260722/overnight.log
mkdir -p "$SWEEP_OUT"

nohup setsid env \
  JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
  CHOPOFF_HUMAN_SWEEP_OUT="$SWEEP_OUT" \
  CHOPOFF_HUMAN_SWEEP_RUNS=5 \
  CHOPOFF_HUMAN_SWEEP_THREADS=24 \
  CHOPOFF_HUMAN_SWEEP_MOTIFS=Cas9,Cas12a \
  CHOPOFF_HUMAN_SWEEP_DISTANCES=0,1,2,3,4 \
  CHOPOFF_HUMAN_SWEEP_REBUILD=0 \
  /home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia \
  --project=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl \
  /home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl/test/local_human/benchmark_human_prefix_sweep.jl \
  </dev/null >"$SWEEP_LOG" 2>&1 &

PID=$!
tail --pid="$PID" -f "$SWEEP_LOG"
```

The parent process launches the database-build stage with one Julia thread and
the search stage with `CHOPOFF_HUMAN_SWEEP_THREADS`.

Outputs under `CHOPOFF_HUMAN_SWEEP_OUT`:

- `builds.csv`: build/reuse status, build time, and index size;
- `timings.csv`: every timed repetition;
- `summary.csv`: median, mean, minimum, and maximum search time;
- `parity.csv`: exact raw prefixHashDB/prefixHashScan result differences;
- `prefixhashscan_stats.csv`: untimed internal scan diagnostics;
- `cases/<motif>_d<distance>/`: final result and parity CSVs;
- `overnight.log`: progress and errors.

Useful controls:

- `CHOPOFF_HUMAN_GENOME`: FASTA path; defaults to the local GRCh38 reference;
- `CHOPOFF_HUMAN_SWEEP_INDEX_PARENT`: database parent directory;
- `CHOPOFF_HUMAN_SWEEP_CAS9_GUIDES`: Cas9 guide file;
- `CHOPOFF_HUMAN_SWEEP_CAS12A_GUIDES`: Cas12a guide file;
- `CHOPOFF_HUMAN_SWEEP_GUIDE_LIMIT`: optional smoke-test guide count;
- `CHOPOFF_HUMAN_SWEEP_RUNS`: timed repetitions per case;
- `CHOPOFF_HUMAN_SWEEP_THREADS`: search threads;
- `CHOPOFF_HUMAN_SWEEP_MOTIFS`: comma-separated motif list;
- `CHOPOFF_HUMAN_SWEEP_DISTANCES`: comma-separated distances within 0-4;
- `CHOPOFF_HUMAN_SWEEP_REBUILD=1`: rebuild existing d4 indexes.

Recompute parity from existing raw outputs without rerunning searches:

```bash
CHOPOFF_HUMAN_SWEEP_STAGE=parity \
CHOPOFF_HUMAN_SWEEP_OUT=/absolute/path/to/prefix_sweep \
CHOPOFF_HUMAN_SWEEP_MOTIFS=Cas9,Cas12a \
CHOPOFF_HUMAN_SWEEP_DISTANCES=0,1,2,3,4 \
JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia \
  --project=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl \
  test/local_human/benchmark_human_prefix_sweep.jl
```

The parity stage updates only `parity.csv`, `summary.csv`, and per-case parity
differences. Existing timing, build, and diagnostic CSVs are not rewritten.
Legacy `prefixhash_validated.csv` and `prefixhash_rejected.csv` files from older
runs are ignored.

For a short sample-genome smoke run:

```bash
SMOKE_ROOT=$(mktemp -d)
CHOPOFF_HUMAN_GENOME=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl/test/sample_data/genome/semirandom.fa \
CHOPOFF_HUMAN_SWEEP_OUT="$SMOKE_ROOT/output" \
CHOPOFF_HUMAN_SWEEP_INDEX_PARENT="$SMOKE_ROOT/indexes" \
CHOPOFF_HUMAN_SWEEP_DISTANCES=0,1 \
CHOPOFF_HUMAN_SWEEP_RUNS=1 \
CHOPOFF_HUMAN_SWEEP_THREADS=2 \
CHOPOFF_HUMAN_SWEEP_GUIDE_LIMIT=1 \
JULIA_DEPOT_PATH=/home/rstudio/livemount/kornel_dev/temp_upload/Soft/julia_depot: \
/home/rstudio/livemount/kornel_dev/temp_upload/Soft/bin/julia \
  --project=/home/rstudio/livemount/kornel_dev/temp_upload/CHOPOFF.jl \
  test/local_human/benchmark_human_prefix_sweep.jl
```
