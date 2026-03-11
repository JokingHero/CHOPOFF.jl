# Sassy Performance Analysis: Julia vs Rust

## Executive Summary

Comparing Julia `src/sassy/` against the reference Rust `sassy/src/`, five architectural gaps explain the observed performance difference. Three dominant issues — early termination, text encoding, and per-call allocations — have been fixed. The remaining gaps are small.

| # | Bottleneck | Est. Impact | Status |
|---|-----------|-------------|--------|
| 1 | Missing early termination | 2–4x | **Fixed** |
| 2 | Scalar text encoding (no SIMD) | 2–3x | **Fixed** (3-tier: AVX2 / scalar-fast / generic) |
| 3 | Per-call heap allocations | 1.3–1.5x | **Fixed** (SassyWorkspace reuse) |
| 4 | Genome string copies | 1.1–1.2x | Open |
| 5 | Dedup/sort/shift overhead | 1.05–1.1x | Open |

Remaining theoretical gap after fixes 1–3: **~1.15–1.3x** (bottlenecks 4-5 combined).

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
- **ntuple closure in eq_store**: Investigated — confirmed **no heap allocation**. Julia unrolls `ntuple` with `Val(LANES)` at compile time.

---

## Bottleneck 3: Per-Call Heap Allocations — FIXED

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

### What was done

Introduced `SassyWorkspace{L}` — an immutable struct holding pre-allocated mutable buffers. All ~9 per-call heap allocations (`hp_store`, `hm_store`, `eq_store`, `current_lane_profiles`, `lane_matches`, `lane_end`, `dist_to_end`, plus result-collection vectors `matches_all`, `unique_matches`, `seen_pos`) are now allocated once and reused via `reset!`.

One workspace per guide, created before the chromosome loop in `search_sassy`, reused for all chromosomes + both strands. Thread-safe: each guide gets its own workspace indexed by `guide_idx` in the `ThreadsX.foreach` closure.

```julia
struct SassyWorkspace{L}
    hp_store::Vector{Vec{L, UInt64}}
    hm_store::Vector{Vec{L, UInt64}}
    eq_store::Vector{Vec{L, UInt64}}
    current_lane_profiles::Matrix{UInt64}   # (L, max_bases=16)
    lane_states_decreasing::Vector{Bool}
    lane_matches::Vector{Vector{Tuple{Int,Int}}}
    lane_end::Vector{Int}
    dist_to_end::Vector{Int}
    matches_all::Vector{Tuple{Int,Int}}
    unique_matches::Vector{Tuple{Int,Int}}
    seen_pos::Set{Int}
end
```

`search_sassy_impl` accepts an optional `workspace` kwarg. When `nothing` (default), allocates fresh buffers for backward compatibility with tests. When provided, calls `reset!` and reuses buffers — `resize!` is effectively a no-op after the first call since all guides share the same `m`.

### Files changed

| File | Change |
|------|--------|
| `src/sassy/core.jl` | Added `SassyWorkspace{L}`, constructor, `reset!`; `search_sassy_impl` accepts `workspace` kwarg |
| `src/sassy/interface.jl` | `_search_impl` closures forward `workspace`; pre-allocates `workspaces` array; `search_sassy_guide` accepts and passes `workspace` |
| `src/sassy/sassy.jl` | Added `SassyWorkspace` export |

### Design notes

- Return value safety: `search_sassy_guide` iterates `ws.unique_matches` immediately in a `for` loop, extracting value-type `(Int, Int)` tuples. No reference escapes before next `reset!`.
- `ntuple` closure in eq_store construction (`core.jl:189`) was investigated and confirmed to **NOT allocate** — Julia unrolls `ntuple` with `Val(LANES)` at compile time. Left as-is.

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

### ~~`ntuple` closure in eq_store construction~~ — NOT AN ISSUE

```julia
eq_store[j] = VecL(ntuple(l -> current_lane_profiles[l, pat_idx], Val(LANES)))
```

Investigated: Julia unrolls `ntuple` with `Val(LANES)` at compile time. The closure is never heap-allocated — it's inlined by the compiler. No action needed.

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

### ~~Phase 3: Allocation Reuse~~ — DONE
- **Expected gain**: 1.3–1.5x
- **Files**: `src/sassy/core.jl`, `src/sassy/interface.jl`, `src/sassy/sassy.jl`
- `ntuple` closure in eq_store confirmed to NOT allocate — no fix needed

### Phase 4: Genome Copy Elimination ← NEXT
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

Target: bring Julia sassy within **1.5–2x** of Rust reference performance (from current ~10x gap, with phases 1–3 already closing the majority of the gap).
