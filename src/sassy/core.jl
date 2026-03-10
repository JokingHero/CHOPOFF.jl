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

    blocks_in_text = cld(text_len, BLOCK_SIZE)
    max_overlap_blocks = cld(m + k, BLOCK_SIZE)
    blocks_per_chunk = max(1, cld(max(blocks_in_text - max_overlap_blocks, 0), LANES))
    lane_chunk_offsets = [(l - 1) * blocks_per_chunk * BLOCK_SIZE for l in 1:LANES]

    VecL = Vec{LANES, UInt64}
    hp_store = fill(VecL(1), m)
    hm_store = fill(VecL(0), m)
    
    # Primitives and Arrays preallocated
    lane_states_decreasing = fill(true, LANES)
    lane_matches = [Tuple{Int, Int}[] for _ in 1:LANES]
    lane_end = zeros(Int, LANES)

    n_bases = length(bases)
    bases_masks = [get_iupac_mask(b) for b in bases]
    
    current_lane_profiles = zeros(UInt64, LANES, n_bases)
    eq_store = Vector{VecL}(undef, m)

    for i in 0:(blocks_per_chunk + max_overlap_blocks - 1)
        vp = VecL(0)
        vm = VecL(0)
        dist_to_start = VecL(0)

        # 1. Direct mask generation mapping (zero allocs)
        for lane in 1:LANES
            lane_end[lane] = lane_chunk_offsets[lane] + (i + 1) * BLOCK_SIZE
            start_pos = lane_chunk_offsets[lane] + i * BLOCK_SIZE + 1
            
            for j in 1:n_bases
                current_lane_profiles[lane, j] = 0
            end

            limit = min(64, text_len - start_pos + 1)
            if limit > 0
                for bit in 0:(limit - 1)
                    pos = start_pos + bit
                    @inbounds c_mask = get_iupac_mask(text[pos])
                    if c_mask != 0
                        bit_mask = UInt64(1) << bit
                        for j in 1:n_bases
                            if (@inbounds bases_masks[j] & c_mask) != 0
                                @inbounds current_lane_profiles[lane, j] |= bit_mask
                            end
                        end
                    end
                end
            end
        end

        # 2. Precompute SIMD vectors for pattern sequence beforehand
        for j in 1:m
            @inbounds pat_idx = pattern_indices[j]
            @inbounds eq_store[j] = VecL(ntuple(l -> current_lane_profiles[l, pat_idx], Val(LANES)))
        end

        # 3. Inner DP rows processing
        for j in 1:m
            @inbounds hp_in = hp_store[j]
            @inbounds hm_in = hm_store[j]
            dist_to_start += hp_in
            dist_to_start -= hm_in

            @inbounds eq = eq_store[j]
            (hp_out, hm_out, vp, vm, _, _) = compute_block(hp_in, hm_in, vp, vm, eq)
            @inbounds hp_store[j] = hp_out
            @inbounds hm_store[j] = hm_out
        end

        # 4. Extract candidates using fast direct index SIMD.Vec
        for lane in 1:LANES
            block_abs_start = lane_chunk_offsets[lane] + i * BLOCK_SIZE
            if block_abs_start >= text_len
                continue
            end
            
            @inbounds state_dec = lane_states_decreasing[lane]
            new_dec = scan_block_minima(
                vp[lane],
                vm[lane],
                Int(dist_to_start[lane]),
                k,
                block_abs_start,
                text_len,
                lane_matches[lane],
                state_dec,
                Val(USE_PEXT)
            )
            @inbounds lane_states_decreasing[lane] = new_dec
        end
    end

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
