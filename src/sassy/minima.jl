# Optimization for Match Extraction
# Corresponds to `sassy/src/minima.rs`

# --- BMI2 PEXT Support (Fast Path for Intel/Modern AMD) ---

# From Julia's generated `features_h.jl`: JL_X86_bmi2 = 32 * 2 + 8.
const JL_X86_BMI2_FEATURE_ID = UInt32(72)
const BMI2_CACHE = Ref{Union{Nothing, Bool}}(nothing)

@inline function can_use_bmi2_pext()
    cached = BMI2_CACHE[]
    if cached !== nothing
        return cached::Bool
    end

    enabled = false
    if Sys.ARCH === :x86_64 || Sys.ARCH === :i686
        enabled = ccall(:jl_test_cpu_feature, Bool, (UInt32,), JL_X86_BMI2_FEATURE_ID)
    end

    BMI2_CACHE[] = enabled
    return enabled
end

"""
    pext_u64(val::UInt64, mask::UInt64) -> UInt64

Wrapper for the x86 BMI2 `pext` instruction using LLVM intrinsics.
Extracts bits from `val` at the positions where `mask` is 1, 
compression them to the least significant bits of the result.
"""
@inline function pext_u64(val::UInt64, mask::UInt64)
    Base.llvmcall(
        (
         """
         declare i64 @llvm.x86.bmi.pext.64(i64, i64)
         define i64 @entry(i64 %0, i64 %1) #0 {
             %res = call i64 @llvm.x86.bmi.pext.64(i64 %0, i64 %1)
             ret i64 %res
         }
         attributes #0 = { "target-features"="+bmi2" }
         """, 
         "entry"
        ),
        UInt64, Tuple{UInt64, UInt64}, val, mask
    )
end

# PACKED_TABLE: Precomputed cost changes for 8-bit PEXT chunks.
const PACKED_TABLE = let
    table = Vector{Tuple{Int8, Int8}}(undef, 256)
    for i in 0:255
        min_val = 0
        cur_val = 0
        for j in 0:7
            bit = (i >> j) & 1
            delta = (bit == 1) ? -1 : 1
            cur_val += delta
            if cur_val < min_val
                min_val = cur_val
            end
        end
        table[i+1] = (Int8(min_val), Int8(cur_val))
    end
    table
end


# --- NIBBLE_TABLE Support (Safe Fallback for Zen 1/Legacy) ---

const NIBBLE_TABLE = let
    table = Vector{Tuple{Int8, Int8}}(undef, 256)
    for i in 0:255
        min_val = 0
        cur_val = 0
        pos_mask = i & 15
        neg_mask = i >> 4
        for j in 0:3
            pos_bit = (pos_mask >> j) & 1
            neg_bit = (neg_mask >> j) & 1
            delta = pos_bit - neg_bit
            cur_val += delta
            if cur_val < min_val
                min_val = cur_val
            end
        end
        table[i+1] = (Int8(min_val), Int8(cur_val))
    end
    table
end


"""
    scan_block_minima(vp_mask, vm_mask, start_cost, k, block_start_pos, text_len, matches_out, lane_state, ::Val{USE_PEXT}, all_minima=false)
"""
@inline function scan_block_minima(
    vp_mask::UInt64,
    vm_mask::UInt64,
    start_cost::Int,
    k::Int,
    block_start_pos::Int,
    text_len::Int,
    matches_out::Vector{Tuple{Int, Int}},
    decreasing::Bool, # Replaced Object param
    ::Val{false},
    all_minima::Bool = false
)::Bool # Add return type
    cost = start_cost
    prev_cost = start_cost
    prev_pos = block_start_pos

    if block_start_pos >= text_len
        return decreasing
    end

    if all_minima && cost <= k && prev_pos == 0
        push!(matches_out, (prev_pos, max(0, cost)))
    end

    bit = 1
    while bit <= 64
        pos = block_start_pos + bit
        if pos > text_len
            break
        end

        if bit <= 61 && !all_minima
            nibble_p = (vp_mask >> (bit - 1)) & 15
            nibble_m = (vm_mask >> (bit - 1)) & 15
            byte = nibble_p | (nibble_m << 4)
            min_cost, end_cost = NIBBLE_TABLE[byte + 1]

            can_skip = (cost + min_cost > k) && !(decreasing && prev_cost <= k)
            if can_skip
                cost += end_cost
                prev_cost = cost
                delta4 = Int((nibble_p >> 3) & 1) - Int((nibble_m >> 3) & 1)
                decreasing = (delta4 <= 0)
                prev_pos = pos + 3
                bit += 4
                continue
            end
        end

        cost += Int((vp_mask >> (bit - 1)) & 0x1)
        cost -= Int((vm_mask >> (bit - 1)) & 0x1)

        if all_minima
            if cost <= k
                push!(matches_out, (pos, max(0, cost)))
            end
        else
            costs_are_equal = cost == prev_cost
            costs_are_increasing = cost > prev_cost
            costs_are_decreasing = cost < prev_cost

            if decreasing && costs_are_increasing && prev_cost <= k
                push!(matches_out, (prev_pos, max(0, prev_cost)))
            end

            decreasing = costs_are_decreasing || (decreasing && costs_are_equal)
        end

        prev_cost = cost
        prev_pos = pos
        bit += 1
    end

    if !all_minima && prev_pos == text_len && decreasing && prev_cost <= k
        push!(matches_out, (prev_pos, max(0, prev_cost)))
    end
    
    return decreasing
end

@inline function scan_block_minima(
    vp_mask::UInt64,
    vm_mask::UInt64,
    start_cost::Int,
    k::Int,
    block_start_pos::Int,
    text_len::Int,
    matches_out::Vector{Tuple{Int, Int}},
    decreasing::Bool,
    ::Val{true},
    all_minima::Bool = false
)::Bool
    if all_minima
        return scan_block_minima(
            vp_mask, vm_mask, start_cost, k, block_start_pos, text_len,
            matches_out, decreasing, Val(false), all_minima
        )
    end

    cost = start_cost
    prev_cost = start_cost
    prev_pos = block_start_pos

    if block_start_pos >= text_len
        return decreasing
    end

    bits_to_scan = min(64, text_len - block_start_pos)
    if bits_to_scan <= 0
        return decreasing
    end

    valid_mask = bits_to_scan == 64 ? typemax(UInt64) : ((UInt64(1) << bits_to_scan) - UInt64(1))
    vp_valid = vp_mask & valid_mask
    vm_valid = vm_mask & valid_mask
    active_mask = vp_valid | vm_valid

    if active_mask == 0
        prev_pos += bits_to_scan
        if prev_pos == text_len && decreasing && prev_cost <= k
            push!(matches_out, (prev_pos, max(0, prev_cost)))
        end
        return decreasing
    end

    packed_pos = pext_u64(vp_valid, active_mask)
    packed_neg = pext_u64(vm_valid, active_mask)

    bit_idx = 0
    event_idx = 0

    while bit_idx < bits_to_scan
        shifted = active_mask >> bit_idx
        if shifted == 0
            prev_pos += bits_to_scan - bit_idx
            bit_idx = bits_to_scan
            break
        end

        zeros_until_event = trailing_zeros(shifted)
        if zeros_until_event > 0
            run_len = min(zeros_until_event, bits_to_scan - bit_idx)
            prev_pos += run_len
            bit_idx += run_len
            if bit_idx >= bits_to_scan
                break
            end
        end

        pos = block_start_pos + bit_idx + 1
        pbit = Int((packed_pos >> event_idx) & UInt64(1))
        nbit = Int((packed_neg >> event_idx) & UInt64(1))
        cost += pbit - nbit

        costs_are_equal = cost == prev_cost
        costs_are_increasing = cost > prev_cost
        costs_are_decreasing = cost < prev_cost

        if decreasing && costs_are_increasing && prev_cost <= k
            push!(matches_out, (prev_pos, max(0, prev_cost)))
        end

        decreasing = costs_are_decreasing || (decreasing && costs_are_equal)

        prev_cost = cost
        prev_pos = pos
        bit_idx += 1
        event_idx += 1
    end

    if prev_pos == text_len && decreasing && prev_cost <= k
        push!(matches_out, (prev_pos, max(0, prev_cost)))
    end
    
    return decreasing
end
