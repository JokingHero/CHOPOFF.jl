using SIMD
include("minima.jl")

"""
    compute_block(hp_in, hm_in, vp, vm, eq_mask)

Pure function representing one step of the Myers bit-vector calculation.
Vectorized to process L lanes simultaneously.
"""
@inline function compute_block(
    hp_in,
    hm_in,
    vp,
    vm,
    eq_mask
)
    # Mirrors Rust `compute_block_simd` exactly (sassy/src/bitpacking.rs).
    vx = eq_mask | vm
    eq = eq_mask | hm_in
    hx = (((eq & vp) + vp) ⊻ vp) | eq
    hp = vm | ~(hx | vp)
    hm = vp & hx

    hp_out = hp >> 63
    hm_out = hm >> 63

    hp_shifted = (hp << 1) | hp_in
    hm_shifted = (hm << 1) | hm_in

    vp_new = hm_shifted | ~(vx | hp_shifted)
    vm_new = hp_shifted & vx

    return (hp_out, hm_out, vp_new, vm_new, hp_shifted, hm_shifted)
end

"""
    search_sassy_impl(pattern_indices, text, k, bases, ::Val{LANES}, ::Val{USE_PEXT})

Executes the Transposed Myers algorithm using SIMD with `LANES` and PEXT settings.
Generic implementation: supports Val{4}/Val{8} Lanes, and Val{true}/Val{false} PEXT.
"""
function search_sassy_impl(
    pattern_indices::Vector{Int},
    text::AbstractVector{UInt8},
    k::Int,
    bases::Vector{UInt8},
    ::Val{LANES} = Val(4),
    ::Val{USE_PEXT} = Val(true)
) where {LANES, USE_PEXT}
    m = Base.length(pattern_indices)
    text_len = Base.length(text)
    if m == 0 || text_len == 0
        return Tuple{Int, Int}[]
    end

    # Match Rust chunking for single-text/single-pattern mode.
    blocks_in_text = cld(text_len, BLOCK_SIZE)
    max_overlap_blocks = cld(m + k, BLOCK_SIZE)
    blocks_per_chunk = max(1, cld(max(blocks_in_text - max_overlap_blocks, 0), LANES))
    lane_chunk_offsets = [(l - 1) * blocks_per_chunk * BLOCK_SIZE for l in 1:LANES]

    VecL = Vec{LANES, UInt64}
    hp_store = fill(VecL(1), m)
    hm_store = fill(VecL(0), m)
    lane_states = [LaneMinimaState(true) for _ in 1:LANES]
    lane_matches = [Tuple{Int, Int}[] for _ in 1:LANES]
    lane_end = zeros(Int, LANES)

    # Outer loop over 64bp blocks for each lane/chunk.
    for i in 0:(blocks_per_chunk + max_overlap_blocks - 1)
        vp = VecL(0)
        vm = VecL(0)
        dist_to_start = VecL(0)

        # Prepare lane profiles for this block.
        current_lane_profiles = Vector{Vector{UInt64}}(undef, LANES)
        for lane in 1:LANES
            lane_end[lane] = lane_chunk_offsets[lane] + (i + 1) * BLOCK_SIZE
            start_pos = lane_chunk_offsets[lane] + i * BLOCK_SIZE + 1
            padded = fill(UInt8('X'), BLOCK_SIZE)
            if start_pos <= text_len
                end_pos = min(start_pos + BLOCK_SIZE - 1, text_len)
                span = end_pos - start_pos + 1
                @inbounds padded[1:span] .= text[start_pos:end_pos]
            end
            current_lane_profiles[lane] = encode_text_profile(padded, bases)
        end

        # Pattern loop (rows of DP matrix).
        for j in 1:m
            hp_in = hp_store[j]
            hm_in = hm_store[j]
            dist_to_start += hp_in
            dist_to_start -= hm_in

            pat_idx = pattern_indices[j]
            eq = VecL(ntuple(l -> current_lane_profiles[l][pat_idx], Val(LANES)))
            (hp_out, hm_out, vp, vm, _, _) = compute_block(hp_in, hm_in, vp, vm, eq)
            hp_store[j] = hp_out
            hm_store[j] = hm_out
        end

        dist_arr = Tuple(dist_to_start)
        vp_arr = Tuple(vp)
        vm_arr = Tuple(vm)
        for lane in 1:LANES
            block_abs_start = lane_chunk_offsets[lane] + i * BLOCK_SIZE
            if block_abs_start >= text_len
                continue
            end
            scan_block_minima(
                vp_arr[lane],
                vm_arr[lane],
                Int(dist_arr[lane]),
                k,
                block_abs_start,
                text_len,
                lane_matches[lane],
                lane_states[lane],
                Val(USE_PEXT)
            )
        end
    end

    # Prune overlapped chunk outputs to mirror Rust `prune_lane_overlaps`.
    if !isempty(lane_matches)
        lane1_end = lane_end[1]
        filter!(x -> x[1] < lane1_end, lane_matches[1])
        for lane in 2:LANES
            prev_lane_end = lane_end[lane - 1]
            cur_lane_end = lane_end[lane]
            if lane == LANES
                filter!(x -> x[1] >= prev_lane_end, lane_matches[lane])
            else
                filter!(x -> (x[1] >= prev_lane_end) && (x[1] < cur_lane_end), lane_matches[lane])
            end
        end
    end

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
    return unique_matches
end
