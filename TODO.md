# Sassy Reimplementation Status (CRISPR-focused)

## Summary

Status as of 2026-02-23: `CHOPOFF.Sassy` is production-ready for the currently defined scope.

The migration objective was:

1. Rust-parity correctness at low-level kernel/search stages.
2. Biologically and format-correct integration in CHOPOFF (`search_sassy`) with parity against `linearDB` on core hit identity.

Both goals are currently met, with explicit non-blocking items listed below.

## Completed and verified

### Kernel and scanner parity

- `compute_block` parity against Rust vectors is fixed and verified.
- `search_sassy_impl` parity against Rust-generated search vectors is fixed and verified.
- AVX-512 lane shape path is validated.

Verification commands:

- `test/verification/julia_test/run_check.jl`
- `test/verification/julia_test/run_search_check.jl`
- `test/verification/julia_test/run_search_check_avx512.jl`

### CHOPOFF integration parity

- Integration parity against `linearDB` is restored for Cas9 and Cas12a on:
  `guide, distance, chromosome, start, strand`.
- Alignment-path differences are treated as diagnostics, not failures, when multiple optimal alignments exist.
- Core regression checks pass in current setup.

Verification commands:

- `julia --project=. -e 'include("test/src/test_sassy_correctness.jl")'`
- `julia --project=. -e 'include("test/src/test_sassy_traceback_parity.jl")'`
- `julia --project=. -e 'ENV["CHOPOFF_VERIFY_TRACEBACK_BACKENDS"]="1"; include("test/src/verify_sassy_core.jl")'`

### Traceback modernization completed

- Added selectable traceback backend in SASSY integration:
  - `:align` (legacy CHOPOFF align)
  - `:custom` (new full DP traceback)
- `:custom` performs full DP-based alignment reconstruction (not distance-only rendering).
- Deterministic traceback tie-break policy is implemented.
- Fast reject precheck is used before full traceback.
- API validates `traceback_backend` values and errors on invalid symbols.

Primary file:

- `src/sassy/interface.jl`

### Minima backend and performance hardening completed

- BMI2/PEXT minima fast path reactivated with runtime feature gating.
- Safe fallback path retained and test coverage strengthened for backend equivalence.
- Hot-loop memory and dispatch issues addressed in the SASSY integration path.

Primary files/scripts:

- `src/sassy/minima.jl`
- `scripts/benchmark_sassy_minima_backend.jl`

## Current guarantees

- Default strict PAM behavior is enabled (`strict_pam=true`).
- Core hit identity parity with `linearDB` is expected on:
  `guide, distance, chromosome, start, strand`.
- Full alignment strings are produced.
- Exact alignment path strings are not guaranteed to match `align` when multiple optimal paths exist.

## Known non-blocking limitations

1. Early stopping is currently strand/chromosome scoped in practice, not strict global cut-off.
2. End-to-end speedup from custom traceback is modest because traceback is no longer the main bottleneck.
3. Further runtime gains now primarily depend on candidate generation and allocation pressure, not DP traceback alone.

## Active optimization backlog

1. True early stopping in `search_sassy_impl`
- Pass stopping criteria into the kernel loop so scanning can terminate early.

2. Reduce traceback call overhead
- Reuse per-thread DP workspaces where safe.
- Reduce duplicated traceback calls inside shift-window candidate expansion.

3. Additional profiling pass
- Re-profile `search_sassy` end-to-end and prioritize hotspots after current traceback changes.

4. Optional future architecture experiment
- Evaluate query-bank/pattern-tiling-like pre-encoding for CHOPOFF-specific workloads if justified by benchmarks.

## Upstream Rust reference policy

Current policy: keep upstream `sassy/` temporarily as a pinned parity/debug reference during stabilization.

Planned direction: remove vendored upstream reference after one stabilization cycle if all gates below hold.

Removal gates:

1. Core parity tests pass without requiring upstream source inspection.
2. Traceback parity and correctness tests remain green.
3. Benchmark gates show no regression on maintained fixtures.
4. Internal CHOPOFF docs are sufficient for maintenance/onboarding.

## Reproducible verification workflow

Primary workflow:

- `scripts/init_sassy.sh`

What it covers:

1. Encoding parity
2. Block-kernel parity
3. Search parity
4. AVX-512 parity

Secondary benchmark workflows:

- `julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl`
- `julia --project=. scripts/benchmark_sassy_minima_backend.jl`
- `julia --project=. scripts/benchmark_sassy_traceback.jl`
