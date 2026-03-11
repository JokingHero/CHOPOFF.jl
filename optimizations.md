# Sassy Performance Analysis: Julia vs Rust

## Executive Summary

Comparing Julia `src/sassy/` against the reference Rust `sassy/src/`, five architectural gaps explain the observed performance difference. Two dominant issues — early termination and text encoding — have been fixed. The remaining gaps are smaller and compound multiplicatively.

| # | Bottleneck | Est. Impact | Status |
|---|-----------|-------------|--------|
| 1 | Missing early termination | 2–4x | **Fixed** |
| 2 | Scalar text encoding (no SIMD) | 2–3x | **Fixed** (3-tier: AVX2 / scalar-fast / generic) |
| 3 | Per-call heap allocations | 1.3–1.5x | Open |
| 4 | Genome string copies | 1.1–1.2x | Open |
| 5 | Dedup/sort/shift overhead | 1.05–1.1x | Open |

Remaining theoretical gap after fixes 1+2: **~1.5–2x** (bottlenecks 3-5 combined).

---

## Bottleneck 2: Text Encoding — FIXED

### What was done

Three-tier encoding replaces the original scalar loop (`core.jl:92-133`):

```
if use_avx2 && limit == BLOCK_SIZE  →  Tier 1: AVX2 vpshufb
elseif n_bases == 4                 →  Tier 2: Scalar fast path (BASE_MATCH LUT)
else                                →  Tier 3: Generic scalar (any n_bases)
```

### Tier 1: AVX2 encoding (`simd_encoding.jl`)

Direct translation of Rust `encode_ref` (profiles/iupac.rs:70-130):

1. `vload` 64 bytes as 2×32-byte `Vec{32, UInt8}` chunks
2. `vpshufb` — parallel IUPAC table lookup via nibble-indexed shuffle (`llvm.x86.avx2.pshuf.b`)
3. Nibble select — blend high/low nibbles based on character index (maps 32-entry IUPAC table through 16-entry pshufb)
4. Per-base: AND with base mask → compare nonzero → `vpmovmskb` → combine halves to UInt64
5. ~40 SIMD ops total for 4 bases × 64 positions

Fires for all full 64-byte blocks on x86-64 with AVX2 (vast majority of genome blocks). Runtime CPU detection via `can_use_avx2()` with cached result.

Key intrinsics via `Base.llvmcall`:
- `vpshufb`: `@llvm.x86.avx2.pshuf.b` — 32-byte parallel table lookup
- `vpmovmskb`: `@llvm.x86.avx2.pmovmskb` — extract MSBs to 32-bit mask

### Tier 2: Scalar fast path (`core.jl`, `constants.jl`)

For `n_bases == 4` on non-AVX2 CPUs or partial blocks (limit < 64):

- `BASE_MATCH[256, 4]` LUT: pre-expanded per-byte → per-base match (0/1)
- 4 register-local accumulators (`m1`–`m4`), no array indexing in inner loop
- Fully branchless: `UInt64(BASE_MATCH[byte+1, j]) << bit`
- 1 text load + 4 LUT reads + 4 shifts + 4 ORs per position

### Tier 3: Generic scalar (`core.jl`)

For patterns with IUPAC ambiguity (n_bases > 4), or non-AVX2 partial blocks:

- Transposed loop: outer=base, inner=bit (better for LLVM pipelining)
- Fully branchless: `UInt64((mask & bm) != 0) << bit`
- Direct assignment (`= mask`) eliminates separate zeroing loop

### Files changed

| File | Change |
|------|--------|
| `src/sassy/simd_encoding.jl` | **New.** AVX2 intrinsics, CPU detection, PACKED_NIBBLES table, `encode_block_avx2!()` |
| `src/sassy/constants.jl` | Added `BASE_MATCH` 256×4 LUT |
| `src/sassy/core.jl` | 3-tier dispatch; `using SIMD` moved to simd_encoding.jl |
| `src/sassy/sassy.jl` | Added `include("simd_encoding.jl")` |

### Impact on other bottlenecks

- **Bottleneck 3 (allocations)**: No change. The encoding uses pre-existing `current_lane_profiles` matrix — zero new allocations.
- **Bottleneck 4 (genome copies)**: No change. AVX2 path uses `vload` from the same `text::AbstractVector{UInt8}`. Eliminating copies would still help by reducing memory pressure / cache pollution.
- **ntuple closure in eq_store**: Still present (`core.jl:136`). Now that encoding is fast, this closure allocation may become a larger relative cost. Prioritize fixing.

---

## Bottleneck 3: Per-Call Heap Allocations (1.3–1.5x)

### What Rust does (`search.rs:221-247`)

`Searcher` struct owns all buffers and reuses them across calls:

```rust
pub struct Searcher<P: Profile> {
    hp: Vec<S>,           // allocated once, reused
    hm: Vec<S>,           // allocated once, reused
    lanes: [LaneState<P>; LANES],  // stack-allocated array
    matches: Vec<Match>,  // cleared + reused
}
```

### What Julia does (`core.jl:62-74`)

Fresh allocations every call to `search_sassy_impl`:

```julia
hp_store = fill(VecL(1), m)                          # heap alloc
hm_store = fill(VecL(0), m)                          # heap alloc
lane_states_decreasing = fill(true, LANES)            # heap alloc
lane_matches = [Tuple{Int,Int}[] for _ in 1:LANES]    # heap alloc × LANES
lane_end = zeros(Int, LANES)                          # heap alloc
current_lane_profiles = zeros(UInt64, LANES, n_bases)  # heap alloc
eq_store = Vector{VecL}(undef, m)                     # heap alloc
```

Called per guide × strand × chromosome. For 100 guides, 25 chroms, 2 strands = **5,000 calls → ~35,000 allocations** + GC pressure.

### Fix

Create a reusable workspace struct:

```julia
mutable struct SassyWorkspace{VecL}
    hp_store::Vector{VecL}
    hm_store::Vector{VecL}
    eq_store::Vector{VecL}
    lane_matches::Vector{Vector{Tuple{Int,Int}}}
    current_lane_profiles::Matrix{UInt64}
    lane_states_decreasing::Vector{Bool}
    lane_end::Vector{Int}
    dist_to_end::Vector{Int}
end

function reset!(ws::SassyWorkspace{VecL}, m::Int) where VecL
    resize!(ws.hp_store, m); fill!(ws.hp_store, VecL(1))
    resize!(ws.hm_store, m); fill!(ws.hm_store, VecL(0))
    resize!(ws.eq_store, m)
    for v in ws.lane_matches; empty!(v); end
    fill!(ws.lane_states_decreasing, true)
    fill!(ws.lane_end, 0)
    fill!(ws.current_lane_profiles, UInt64(0))
end
```

Allocate once in `search_sassy`, pass to all calls. Thread-local copies for threaded search.

### Priority after encoding fix: **HIGH**

Now that encoding is fast, allocation/GC overhead is a larger fraction of total time. This is the next bottleneck to address.

---

## Bottleneck 4: Genome String Copies (1.1–1.2x)

### What happens (`interface.jl:535,356-357`)

```julia
seq_str = String(seq)                    # copy 1: LongDNA → String (~250MB for chr1)
# Then inside search_sassy_guide:
genome_bytes = codeunits(genome_str)     # view (ok)
genome_seq = LongDNA{4}(genome_str)      # copy 2: String → LongDNA (~250MB)
```

Two full copies of each chromosome. For hg38 (3.1GB), that's ~6.2GB of extra allocation.

### Fix

```julia
# In search_sassy: pass seq directly, avoid String conversion
seq_bytes = Vector{UInt8}(String(seq))  # single copy, reuse everywhere

# In search_sassy_guide: accept bytes + original seq
function search_sassy_guide(guide_seq, genome_bytes::Vector{UInt8},
                            genome_seq::LongDNA{4}, ...)
```

Or better — investigate if BioSequences provides direct byte access via `unsafe_wrap` or similar.

### Priority: MEDIUM

Memory bandwidth matters more now that CPU is less bottlenecked by encoding. Reducing copies improves cache utilization for the AVX2 encoding path.

---

## Bottleneck 5: Deduplication & Post-Processing (1.05–1.1x)

### Current (`core.jl:163-177`)

```julia
matches_all = Tuple{Int, Int}[]
for lane in 1:LANES
    append!(matches_all, lane_matches[lane])
end
unique_matches = Tuple{Int, Int}[]
seen_pos = Set{Int}()
sort!(matches_all, by = x -> (x[1], x[2]))
for mt in matches_all
    if !(mt[1] in seen_pos)
        push!(unique_matches, mt)
        push!(seen_pos, mt[1])
    end
end
```

### Rust approach (`search.rs:1086-1123`)

In-place `retain()` on each lane's matches — no sort, no Set, no concatenation:

```rust
self.lanes[0].matches.retain(|&(end_pos, _)| end_pos < cur_lane_end);
for lane in 1..LANES {
    self.lanes[lane].matches.retain(|&(end_pos, _)| {
        end_pos >= prev_lane_end && (lane == LANES-1 || end_pos < cur_lane_end)
    });
}
```

### Fix

Replace with in-place filtering + direct iteration:

```julia
# Filter each lane in-place (already done in lines 149-161, but then re-collected)
# Just iterate lane_matches directly instead of creating matches_all + unique_matches
```

### Priority: LOW

Small absolute impact, but becomes easier once Bottleneck 3 (workspace struct) is done.

---

## Additional Observations

### `ntuple` closure in eq_store construction

```julia
eq_store[j] = VecL(ntuple(l -> current_lane_profiles[l, pat_idx], Val(LANES)))
```

This creates a closure capturing `current_lane_profiles` and `pat_idx` on every iteration. With m=20 and ~46M blocks = **~920M closure allocations**. Should be replaced with direct construction or a helper that avoids the closure.

**Priority after encoding fix: HIGH** — this is now one of the top remaining per-block costs.

### Shift loop in `search_sassy_guide`

```julia
for shift in -k:k  # 2k+1 iterations per match
    trial_end = global_reported_end + shift
    # ... PAM check + traceback for each shift ...
end
```

Rust integrates PAM filtering into the search via `filter_fn`, avoiding this post-hoc window scan. This is a correctness-sensitive area — the shift compensates for approximate position reporting from sassy. Consider whether better position tracking inside sassy could eliminate this loop.

### SIMD.jl Vec verification

Need to verify with `@code_llvm` that `Vec{4, UInt64}` operations in `compute_block` actually emit AVX2 instructions (`vpaddq`, `vpand`, `vpor`, etc.) rather than scalar fallbacks. If LLVM isn't vectorizing, the inner DP loop itself could be 4x slower than expected.

```julia
# Quick check:
using SIMD
f(a, b) = a + b
@code_llvm f(Vec{4, UInt64}((1,2,3,4)), Vec{4, UInt64}((5,6,7,8)))
# Should show: <4 x i64> add instruction
```

---

## Recommended Implementation Order (Updated)

### ~~Phase 1: Early Termination~~ — DONE
- **Expected gain**: 2–4x

### ~~Phase 2: Text Encoding~~ — DONE
- **Expected gain**: 2–3x (AVX2 path matches Rust encode_ref)
- **Files**: `simd_encoding.jl` (new), `constants.jl`, `core.jl`, `sassy.jl`

### Phase 3: Allocation Reuse ← NEXT
- **Files**: `src/sassy/core.jl`, `src/sassy/interface.jl`
- **Expected gain**: 1.3–1.5x
- **Risk**: Low — structural refactor, no algorithm change
- **Includes**: Fix ntuple closure in eq_store (same refactor)

### Phase 4: Genome Copy Elimination
- **Files**: `src/sassy/interface.jl`
- **Expected gain**: 1.1–1.2x
- **Risk**: Low

### Phase 5: Post-processing Cleanup
- **Files**: `src/sassy/core.jl`, `src/sassy/interface.jl`
- **Expected gain**: 1.05–1.1x
- **Risk**: Low

---

## Verification Plan

After each phase:
1. `julia --project=. -e 'include("test/src/test_sassy_correctness.jl")'` — sassy unit tests (15 groups)
2. `julia --project=. -e 'include("test/src/verify_sassy_core.jl")'` — linearDB parity
3. `CHOPOFF_BENCH_RUNS=7 julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl` — benchmark

Target: bring Julia sassy within **1.5–2x** of Rust reference performance (from current ~10x gap, with phases 1+2 already closing ~4–8x).
