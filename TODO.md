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
  - `test/src/verify_sassy_core.jl` now completes (fails on result mismatch, not runtime lockup).

### Still Failing
- `verify_sassy_core.jl` (`linearDB vs sassy`) reports many output mismatches (coordinates/orientation/alignment semantics), across d=1..3.

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
2. Alignment string orientation audit:
   - Match `linearDB` output orientation conventions for `alignment_guide` and `alignment_reference`.
   - Ensure Cas9/Cas12a strand handling is deterministic and documented.
3. Distance/alignment tie-break policy:
   - In shift window search (`-k:k`), define deterministic tie-break (distance first, then preferred positional rule).

Acceptance gate:
- `test/src/verify_sassy_core.jl` passes (`linearDB vs sassy` parity restored for Cas9 fixture).

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

### M5. PEXT Fast Path Re-activation (after integration parity)
Target file:
- `src/sassy/minima.jl`

Actions:
- Keep nibble path as correctness baseline.
- Reintroduce true BMI2/PEXT acceleration only after proving equivalence to nibble path.
- Feature-gate cleanly and preserve `force_safe_minima` behavior.

Acceptance gate:
- `PEXT` and nibble outputs identical on regression corpus.

### M6. Performance Hardening (final)
Actions:
- Profile `search_sassy` and identify remaining hotspots (alignment calls, profile encoding, allocation pressure).
- Optional: pre-encode reference blocks if profiler confirms meaningful wins without semantic risk.
- Ensure no debug prints remain in hot loops.

Acceptance gate:
- Runtime stable and scalable on semirandom fixture plus at least one larger genome slice.

## Reproducible Verification Workflow
Use:
- `scripts/init_sassy.sh`

It runs:
1. encoding parity
2. block kernel parity
3. search parity
4. AVX-512 path parity

## Immediate Next Fixes (ordered)
1. Re-activate true BMI2/PEXT path in `src/sassy/minima.jl` with nibble-equivalence gates.
2. Profile `search_sassy` and reduce allocation pressure in alignment-heavy sections.
3. Add explicit Cas12a integration parity checks against `linearDB`.
