using JSON
using Test
using CHOPOFF.Sassy

println("Loading Search Test Vectors...")
path = joinpath(@__DIR__, "../test_vectors_search.json")
if !isfile(path)
    println("Error: $path not found. Run 'cargo run' in ../rust_gen first.")
    exit(1)
end
const vectors = JSON.parsefile(path)

println("Running Validation on search_sassy_impl...")

@testset "Sassy Search Verification" begin
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
        
        matches = search_sassy_impl(indices, text, k, bases)
        
        # Rust ground truth set (Position, Cost)
        rust_set = Set{Tuple{Int, Int}}()
        for m in expected_matches
             # Rust matches are 0-based exclusive end.
             # Julia matches are 1-based inclusive end.
             # "ACGT" (0..4) -> End 4.
             # "ACGT" (1..4) -> End 4.
             # These should match exactly.
            push!(rust_set, (m["text_end"], m["cost"]))
        end
        
        # For comparison, we filter Julia matches to keep only the BEST cost for each position,
        # because search_sassy_impl currently might return non-optimal costs if not sorted by cost?
        # Actually search_sassy_impl does some dedup, but let's see what it produces.
        
        julia_set = Set(matches)
        
        # Check if we are missing anything Rust found.
        # Implies Rust found a valid match <= k that we missed.
        
        # We allow Julia to find MORE matches (since Rust filters local minima), 
        # but Julia MUST find at least everything Rust found (modulo local minima overlapping?).
        # Actually, if Rust filters non-local minima, it might return FEWER matches.
        # If Julia returns ALL matches <= k (it seems to try to), then RustSet should be a specific subset of JuliaSet.
        
        # However, check for exact matches first.
        
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
            
            # Print nearby julia matches
            println("Julia matches dump:")
            for jm in matches
                if any(rm -> abs(rm[1] - jm[1]) < 5, missing)
                    println("  Found: $jm")
                end
            end
            
            exit(1)
        end
    end
end
println("Verification Passed!")
