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
