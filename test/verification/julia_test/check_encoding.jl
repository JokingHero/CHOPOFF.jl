using JSON
using Test
include("../../../src/sassy/constants.jl")
include("../../../src/sassy/encoding.jl")

println("Loading Encoding Test...")
const data = JSON.parsefile("../encoding_test.json")

@testset "Encoding Logic" begin
    text_str = data["text"]
    text_bytes = Vector{UInt8}(text_str)
    
    # Replicate 64-byte padding logic
    padded = fill(UInt8('X'), 64)
    padded[1:length(text_bytes)] .= text_bytes

    # Rust provided the expected masks for [A, C, T, G]
    expected_masks = UInt64.(data["encoded"])
    bases_to_check = [UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')]

    # Run Julia Logic
    actual_masks = encode_text_profile(padded, bases_to_check)

    # Compare
    for (i, base) in enumerate(bases_to_check)
        char_base = Char(base)
        println("Checking mask for Base: $char_base")
        
        if actual_masks[i] != expected_masks[i]
            println("MISMATCH!")
            println("Expected: $(bitstring(expected_masks[i]))")
            println("Actual:   $(bitstring(actual_masks[i]))")
            
            # Find the bit difference
            diff = actual_masks[i] ⊻ expected_masks[i]
            pos = trailing_zeros(diff) + 1
            println("First diff at index: $pos (Character: $(Char(padded[pos])))")
        end
        @test actual_masks[i] == expected_masks[i]
    end
end
