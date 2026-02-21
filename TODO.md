# Sassy Reimplementation TODO (CRISPR-focused)

## Goal

Make `CHOPOFF.Sassy` both:

1. Rust-parity correct at the low-level kernel/search stage.
2. Biologically/format correct at the CHOPOFF integration stage (`search_sassy` output parity with `linearDB` where intended).

## Verified Current Status

### Done

- `compute_block` parity fixed against Rust vectors.
  - Command: `test/verification/julia_test/run_check.jl`
  - Result: 40,000 / 40,000 passed.
- `search_sassy_impl` parity fixed against Rust-generated search vectors.
  - Command: `test/verification/julia_test/run_search_check.jl`
  - Result: passed.
- AVX-512 lane shape path validated.
  - Command: `test/verification/julia_test/run_search_check_avx512.jl`
  - Result: passed.
- Candidate explosion / apparent "runs forever" issue removed.
  - `test/src/verify_sassy_core.jl` now completes.
- Integration parity against linearDB restored for Cas9 and Cas12a for off-target identity.
  - Command: `test/src/verify_sassy_core.jl`
  - Result: full parity on `guide, distance, chromosome, start, strand` across `d=1:3`.
- Full package test suite now passes end-to-end.
  - Command: `julia --project=. -e 'cd("test"); include("runtests.jl")'`
  - Result: completed with exit code 0.
- Legacy IUPAC test API compatibility restored.
  - Added `TextBlockProfile`, `encode_text_block`, and exported `get_iupac_mask` via `CHOPOFF.Sassy`.
- `verify_sassy_core` test working directory handling made robust.
  - Uses `cd(dirname(dirname(@__FILE__)))` instead of hardcoded relative `cd("test")`.
- Alignment-path differences are now treated as diagnostics (non-blocking), not correctness failures.
- Added first-mismatch isolation helper:
  - `scripts/debug_sassy_mismatch.jl --motif Cas9|Cas12a --distance 1|2|3`
- BMI2/PEXT minima fast path reactivated with runtime feature gating.
  - Auto mode: PEXT on BMI2 x86 hosts, safe minima fallback otherwise.
  - User override retained via `force_safe_minima` (API) and `--force_safe_minima` (CLI).
- PEXT-vs-safe identity checks strengthened.
  - `test/src/test_sassy_correctness.jl` now asserts full tuple identity on
    `guide, alignment_guide, alignment_reference, distance, chromosome, start, strand`
    across mixed-distance, boundary, and no-hit fixtures.
- Minima backend benchmark/report script added.
  - `scripts/benchmark_sassy_minima_backend.jl`
  - Current host report: no regression (`auto` median <= `safe` median).
- Memory Allocation in Hot Loop (OOM Risk) resolved via zero-copy `codeunits()`.
- Redundant Thread Lock causing thread starvation removed (safe parallel appends confirmed).
- Type Instability & Dynamic Dispatch in `search_sassy_guide` resolved via parameterized `impl_func::F`.
- String Allocation on I/O mitigated via direct `print()` chunking instead of concatenation.
- Redundant SIMD Work for Sequence overlap removed by enforcing `max(1, blocks_per_chunk)`.
- `PACKED_TABLE` and `NIBBLE_TABLE` optimizations incorporated into the SIMD fallback scan path (`scan_block_minima`), yielding an 8.5x overall `search_sassy` speedup vs `prefixHashDB`.

### Still Failing

### 1. Algorithmic Logic Flaws

#### A. "Early Stopping" Isn't Actually Stopping Early

In `search_sassy`, you check `is_es` like this:

```julia
if is_es; return; end
# ...
results_fwd = search_sassy_guide(...)
process_results!(results_fwd, true)
```

The problem is `search_sassy_guide` (and underneath, `search_sassy_impl`) has no awareness of your early stopping limits. It processes the **entire chromosome** and returns all matches before `process_results!` can tally them and flip `is_es = true`.
If your limit is 10, but a guide matches 50,000 times on chromosome 1, the algorithm will compute, align, and allocate all 50,000 matches before realizing it should have stopped at 10.

**The Fix (FUTURE TODO):**
For true early stopping, you must pass the limits down into `search_sassy_impl` so the DP matrix loop can `break` early. For the initial reimplementation, we are leaving this as "Early stopping per-chromosome-strand," not a strict absolute cut-off, to avoid adding complex branching to the tight SIMD loop.

## Execution Strategy (Decision-Complete)

### M1. Kernel/Scanner Parity Lock (completed)

- Keep `compute_block` semantics aligned with Rust `bitpacking.rs::compute_block_simd`.
- Keep minima extraction state-machine (local minima with `decreasing` lane state) aligned with Rust `search.rs::find_minima_with_overhang` (non-overhang branch).
- Keep lane overlap pruning behavior aligned with Rust `search.rs::prune_lane_overlaps`.

### M2. Strict PAM Policy and Guide-Only Search (completed default behavior)

- Default: `strict_pam=true` in `search_sassy`.
- Search pattern under strict mode is guide-only (strand-normalized), not guide+PAM.
- PAM is filtered with IUPAC mask matching at candidate coordinates before alignment.
- Keep opt-out: `strict_pam=false` for fuzzy legacy behavior.

### M3. Integration Parity with CHOPOFF Output (completed)

Primary target files:

- `src/sassy/interface.jl`
- `src/db_linear.jl` (reference semantics)

Actions:

1. Coordinate convention audit:
   - Ensure `start` coordinate aligns with `linearDB` conventions for both strands.
   - Validate off-by-one and guide-span logic against known synthetic fixtures.
2. Alignment orientation audit:
   - Keep biologically correct strand/orientation behavior for `alignment_guide` and `alignment_reference`.
   - Do not require exact path-string identity when multiple optimal alignments exist.
   - Ensure Cas9/Cas12a strand handling is deterministic and documented.
3. Distance/alignment tie-break policy:
   - In shift window search (`-k:k`), define deterministic tie-break (distance first, then preferred positional rule).

Acceptance gate:

- `test/src/verify_sassy_core.jl` passes for both Cas9 and Cas12a with parity on:
  - `guide, distance, chromosome, start, strand`
- Alignment-path deltas are printed as diagnostics and are not test failures.

### M4. Expand Regression Tests for Integration Semantics (completed baseline)

Target file:

- `test/src/test_sassy_correctness.jl`

Add/keep explicit tests for:

- strict PAM acceptance/rejection (including IUPAC PAM chars).
- coordinate correctness on known placements.
- strand-specific alignment orientation for Cas9 and Cas12a.
- local-minima uniqueness (no duplicate nearby hits for one true site).

Acceptance gate:

- `test/src/test_sassy_correctness.jl` passes with strict PAM default and CHOPOFF coordinate conventions.

### M5. PEXT Fast Path Re-activation (completed)

Target file:

- `src/sassy/minima.jl`

Actions:

- Kept nibble path as correctness baseline.
- Reintroduced BMI2/PEXT minima scan path.
- Added runtime CPU gating (`bmi2`) with safe fallback and preserved `force_safe_minima`.
- Extended equivalence tests to enforce tuple identity between backends.

Acceptance gate:

- `PEXT` and nibble outputs are identical on regression corpus (pass).

### M6. Performance Hardening (final)

Actions:

- Profile `search_sassy` and identify remaining hotspots (alignment calls, profile encoding, allocation pressure).
- Optional: pre-encode reference blocks if profiler confirms meaningful wins without semantic risk.
- Ensure no debug prints remain in hot loops.

Acceptance gate:

- Runtime stable and scalable.
- Multithreaded & single threaded performance verified via `benchmark_sassy_vs_prefixhash.jl`
- Verified overall ~8.5x speedup compared to `prefixHashDB`.

## Reproducible Verification Workflow

Use:

- `scripts/init_sassy.sh`

It runs:

1. encoding parity
2. block kernel parity
3. search parity
4. AVX-512 path parity

## Speed Benchmark Workflow

Use:

- `julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl`

What it runs by default:

1. `sassy` vs `prefixHashDB` search-only timings.
2. Two thread modes: single-thread (`JULIA_NUM_THREADS=1`) and current env-thread mode.
3. Cas9 + Cas12a fixtures at distance 3.
4. Core tuple parity checks on `guide, distance, chromosome, start, strand`.

Useful options:

- `CHOPOFF_BENCH_RUNS=11` (more timing repetitions)
- `CHOPOFF_BENCH_MODE=single|env|both` (default `both`)
- `CHOPOFF_BENCH_OUT=/tmp/chopoff_speed_summary.csv` (write combined summary)
- `CHOPOFF_BENCH_KEEP_TMP=1` (keep temporary artifacts for debugging)
