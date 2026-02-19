"""
    encode_text_profile(text, bases)
    
Matches Rust's `encode_ref` logic.
Input: A 64-byte text slice.
Output: Vector of UInt64 masks (one per base).


THis was tested:
cd test/verification/julia_test
julia check_encoding.jl
Loading Encoding Test...
Checking mask for Base: A
Checking mask for Base: C
Checking mask for Base: T
Checking mask for Base: G
Test Summary:  | Pass  Total  Time
Encoding Logic |    4      4  0.2s

"""
struct TextBlockProfile
    masks::Vector{UInt64}
end

function encode_text_profile(text::AbstractVector{UInt8}, bases::Vector{UInt8})
    n_bases = length(bases)
    masks = zeros(UInt64, n_bases)
    len = length(text)

    for (j, base) in enumerate(bases)
        base_mask = get_iupac_mask(base)
        result = UInt64(0)
        for i in 1:min(64, len)
            text_char = text[i]
            text_mask = get_iupac_mask(text_char)
            if (base_mask & text_mask) != 0
                result |= (UInt64(1) << (i - 1))
            end
        end
        masks[j] = result
    end
    return masks
end

@inline function encode_text_block(text::AbstractVector{UInt8}, bases::Vector{UInt8})
    return TextBlockProfile(encode_text_profile(text, bases))
end

# --- Pattern Encoding ---

function encode_pattern_sassy(pattern::AbstractVector{UInt8})
    bases = UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')]
    indices = Int[]
    for c in pattern
        c_upper = c & ~UInt8(0x20)
        idx = findfirst(==(c_upper), bases)
        if idx === nothing
            push!(bases, c_upper)
            push!(indices, length(bases))
        else
            push!(indices, idx)
        end
    end
    return (bases, indices)
end
