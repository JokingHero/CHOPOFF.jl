# Sassy Performance Analysis: Julia vs Rust (~10x gap)

## Executive Summary

Comparing Julia `src/sassy/` against the reference Rust `sassy/src/`, five architectural gaps explain the observed ~10x performance difference. Two dominant issues — missing early termination and scalar text encoding — account for the bulk. The remaining gaps compound multiplicatively.

| # | Bottleneck | Est. Impact | Complexity |
|---|-----------|-------------|------------|
| 1 | Missing early termination | 2–4x | Medium |
| 2 | Scalar text encoding (no SIMD) | 2–3x | Medium-High |
| 3 | Per-call heap allocations | 1.3–1.5x | Low |
| 4 | Genome string copies | 1.1–1.2x | Low |
| 5 | Dedup/sort/shift overhead | 1.05–1.1x | Low |

Combined theoretical: **3.4–19.8x** → realistic target: **4–8x improvement**.

---

## Bottleneck 1: Missing Early Termination (2–4x)

### What Rust does (`search.rs:996-1046`)

After each pattern row `j`, Rust checks whether ANY SIMD lane still has `dist_to_end ≤ k`. If not, and `j > prev_end_last_below`, it calls `check_lanes()` which computes `min_in_lane()` (prefix minimum of vp/vm bit-vectors). If no lane is promising:

```rust
self.reset_rows(j + 1, prev_max_j);   // reset hp/hm for skipped rows
prev_end_last_below = cur_end_last_below.max(CHECK_AT_LEAST_ROWS);
continue 'text_chunk;                   // skip to next 64-byte block
```

For CRISPR search with k=3 and pattern length m=20, the vast majority of 64-byte genome blocks contain no approximate match. Rust processes ~8–10 rows before bailing; Julia processes all 20.

### What Julia does (`core.jl:114-124`)

```julia
for j in 1:m
    # ... compute_block ...
    # NO early termination check — always runs all m iterations
end
```

### Why it matters

Human genome ≈ 3 billion bases → ~46 million 64-byte blocks. With m=20:
- **Julia**: 46M × 20 = **920M** `compute_block` calls
- **Rust**: 46M × ~10 (avg with early exit) = **~460M** calls, plus many blocks fully skipped

### Fix

Add `dist_to_end` tracking, `prefix_min_val()` helper, and `reset_rows!()` to `core.jl`:

```julia
# After compute_block in inner loop:
dist_to_end += hp_out
dist_to_end -= hm_out

# Periodically check if any lane is still promising
if any_lane_below_threshold(dist_to_end, k)
    cur_end_last_below = j
end

if j > prev_end_last_below
    promising = false
    for lane in 1:LANES
        if dist_to_start[lane] + prefix_min_val(vp[lane], vm[lane]) <= k
            promising = true; prev_end_last_below = j; break
        end
    end
    if !promising
        reset_rows!(hp_store, hm_store, j+1, prev_max_j, VecL)
        prev_end_last_below = max(cur_end_last_below, CHECK_AT_LEAST_ROWS)
        prev_max_j = j
        @goto next_block
    end
end
```

New helper using existing `NIBBLE_TABLE`:

```julia
@inline function prefix_min_val(vp::UInt64, vm::UInt64)::Int
    min_val = Int(0); cur = Int(0)
    for i in 0:15
        byte = UInt8(((vp >> (i*4)) & 0xF) | (((vm >> (i*4)) & 0xF) << 4))
        @inbounds tmin, tend = NIBBLE_TABLE[byte + 1]
        min_val = min(min_val, cur + Int(tmin))
        cur += Int(tend)
    end
    return min_val
end
```

---

## Bottleneck 2: Scalar Text Encoding (2–3x)

### What Rust does (`profiles/iupac.rs:70-130`)

SIMD-vectorized encoding of 64-byte text blocks into per-base bitmasks:

```rust
// Load 32 bytes into AVX2 register
let chunk0 = u8x32::from(&b[0..32]);
let chunk1 = u8x32::from(&b[32..64]);

// IUPAC lookup via AVX2 pshufb (parallel table lookup for 32 bytes)
let shuffled0 = half_shuffle(tbl256, low4_0);  // __mm256_shuffle_epi8
// ... blend high/low nibbles, compare, extract bitmask
let low = match0.to_bitmask() as u32 as u64;   // movmskb → u64
let high = match1.to_bitmask() as u32 as u64;
*out.get_unchecked_mut(i) = (high << 32) | low;
```

Total: ~40 SIMD operations for 4 bases × 64 positions.

### What Julia does (`core.jl:82-105`)

Scalar loop with per-bit, per-base, per-lane iteration:

```julia
for lane in 1:LANES           # 4
    for bit in 0:(limit-1)     # up to 64
        c_mask = get_iupac_mask(text[pos])
        for j in 1:n_bases     # 4+
            if (bases_masks[j] & c_mask) != 0
                current_lane_profiles[lane, j] |= bit_mask
            end
        end
    end
end
```

Total: **4 × 64 × 4 = 1,024 scalar iterations** with conditional branches.

### Why it matters

Per block: ~1,024 scalar ops (Julia) vs ~40 SIMD ops (Rust) = **~25x per-block overhead**.

This runs for every block in the genome (~46M blocks), making it the dominant per-block cost.

### Fix (incremental)

**Step 1** — Unroll base loop for common ACGT case, eliminate inner branch:

```julia
for bit in 0:(limit-1)
    @inbounds c_mask = IUPAC_TABLE[text[start_pos + bit] + 1]
    c_mask == 0 && continue
    bit_mask = UInt64(1) << bit
    # Branchless for 4 standard bases
    current_lane_profiles[lane, 1] |= ifelse((c_mask & bases_masks[1]) != 0, bit_mask, zero(UInt64))
    current_lane_profiles[lane, 2] |= ifelse((c_mask & bases_masks[2]) != 0, bit_mask, zero(UInt64))
    current_lane_profiles[lane, 3] |= ifelse((c_mask & bases_masks[3]) != 0, bit_mask, zero(UInt64))
    current_lane_profiles[lane, 4] |= ifelse((c_mask & bases_masks[4]) != 0, bit_mask, zero(UInt64))
end
```

**Step 2** — Process 8 positions at a time with precomputed lookup:

```julia
# Precompute 256-entry table: byte → tuple of 4 base-match flags
const ENCODING_LUT = let
    t = Vector{NTuple{4,Bool}}(undef, 256)
    A, C, T, G = UInt8(1), UInt8(2), UInt8(4), UInt8(8)
    for i in 0:255
        m = IUPAC_TABLE[i+1]
        t[i+1] = ((m & A) != 0, (m & C) != 0, (m & T) != 0, (m & G) != 0)
    end
    t
end
```

**Step 3** (future) — Julia SIMD.jl or llvmcall for AVX2 pshufb:

```julia
# Requires careful LLVM IR generation for _mm256_shuffle_epi8
# and _mm256_movemask_epi8 equivalents
```

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
    // ...
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

---

## Additional Observations

### `ntuple` closure in eq_store construction

```julia
eq_store[j] = VecL(ntuple(l -> current_lane_profiles[l, pat_idx], Val(LANES)))
```

This creates a closure capturing `current_lane_profiles` and `pat_idx` on every iteration. With m=20 and ~46M blocks = **~920M closure allocations**. Should be replaced with direct construction or a helper that avoids the closure.

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

## Recommended Implementation Order

### Phase 1: Early Termination
- **Files**: `src/sassy/core.jl`, `src/sassy/minima.jl`
- **Expected gain**: 2–4x
- **Risk**: Medium — must preserve correctness at chunk boundaries
- **Test**: `test/src/test_sassy_correctness.jl`, `test/src/verify_sassy_core.jl`

### Phase 2: Allocation Reuse
- **Files**: `src/sassy/core.jl`, `src/sassy/interface.jl`
- **Expected gain**: 1.3–1.5x
- **Risk**: Low — structural refactor, no algorithm change
- **Test**: Same test suite

### Phase 3: Text Encoding Optimization
- **Files**: `src/sassy/core.jl`
- **Expected gain**: 1.5–2x (unrolled scalar), 2–3x (SIMD)
- **Risk**: Low for scalar unroll, high for SIMD intrinsics
- **Test**: `test/src/test_sassy_correctness.jl` encoding tests

### Phase 4: Genome Copy Elimination
- **Files**: `src/sassy/interface.jl`
- **Expected gain**: 1.1–1.2x
- **Risk**: Low
- **Test**: Full test suite

### Phase 5: Post-processing Cleanup
- **Files**: `src/sassy/core.jl`, `src/sassy/interface.jl`
- **Expected gain**: 1.05–1.1x
- **Risk**: Low
- **Test**: Full test suite

---

## Verification Plan

After each phase:
1. `julia --project=. test/runtests.jl` — full test suite
2. `julia --project=. -e 'include("test/src/test_sassy_correctness.jl")'` — sassy unit tests
3. `julia --project=. -e 'include("test/src/verify_sassy_core.jl")'` — linearDB parity
4. `CHOPOFF_BENCH_RUNS=7 julia --project=. scripts/benchmark_sassy_vs_prefixhash.jl` — benchmark

Target: bring Julia sassy within 2–3x of Rust reference performance (from current ~10x gap).
