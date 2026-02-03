using JSON
using Test

# Include Sassy implementation
include("../../../src/sassy/constants.jl")
include("../../../src/sassy/encoding.jl")
include("../../../src/sassy/core.jl")

println("Loading Search Test Vectors...")
path = joinpath(@__DIR__, "../test_vectors_search.json")
# Reuse existing vectors
if !isfile(path)
    # Fallback if running from a different dir context
    path = joinpath(@__DIR__, "../rust_gen/test_vectors_search.json")
end
const vectors = JSON.parsefile(path)

println("Running Validation on search_sassy_impl with Val(8) (AVX-512 emulation)...")

@testset "Sassy Search Verification AVX-512" begin
    for (i, case) in enumerate(vectors)
        inp = case["input"]
        pattern_str = inp["pattern"]
        text_str = inp["text"]
        k = inp["k"]
        
        expected_matches = case["output"]

        # Julia Execution
        pattern = Vector{UInt8}(pattern_str)
        text = Vector{UInt8}(text_str)
        (bases, indices) = encode_pattern_sassy(pattern)
        
        # EXPLICITLY CALL WITH Val(8) to test 8-lane Simd
        matches = search_sassy_impl(indices, text, k, bases, Val(8))
        
        # Rust ground truth set (Position, Cost)
        rust_set = Set{Tuple{Int, Int}}()
        for m in expected_matches
            push!(rust_set, (m["text_end"], m["cost"]))
        end
        
        julia_set = Set(matches)
        
        missing = Set{Tuple{Int, Int}}()
        for rm in rust_set
            # check if rm is in julia_set
            if !(rm in julia_set)
                # It might be in julia set with a DIFFERENT cost?
                # Check if there is any match at rm[1] with cost <= rm[2]
                better_or_equal = false
                for jm in matches
                    if jm[1] == rm[1] && jm[2] <= rm[2]
                        better_or_equal = true
                        break
                    end
                end
                if !better_or_equal
                    push!(missing, rm)
                end
            end
        end
        
        if !isempty(missing)
            println("\n[FAIL] Case #$i")
            println("Pattern: $pattern_str")
            println("Text:    $text_str")
            println("K:       $k")
            println("Missing Matches (Rust found, Julia missed or found worse cost): ", missing)
            exit(1)
        end
    end
end
println("AVX-512 Verification Passed!")
