using Test
using CHOPOFF
using BioSequences
using DataFrames
using CSV

# Get the test directory for file operations
const TEST_DIR = @__DIR__

@testset "Sassy Algorithm Tests" begin

    # Mock Genome file - using absolute paths
    # For distance=3 or 4, leading/trailing chars are okay
    genome_seq = "AAAA" * 
                 "TTTTTTTTTTTTTTTTTTTT" * "GGG" *   # Match: 20T + GGG (Perfect Cas9 with N=G at PAM)
                 "AAAA"
    # Total: 4 + 20 + 3 + 4 = 31 chars
                 
    # 20T target
    # Guide (Spacer): TTTTTTTTTTTTTTTTTTTT
    # PAM: NGG -> pattern is 20T + NGG = 23 chars
    # The genome substring at positions 5-27 matches pattern exactly (edit distance 0)
    
    # Use absolute paths
    genome_path = joinpath(TEST_DIR, "test_genome_sassy.fa")
    guide_path = joinpath(TEST_DIR, "test_guides_sassy.txt")
    output_path = joinpath(TEST_DIR, "sassy_results.csv")
    
    # helper to clean up
    function cleanup()
        rm(genome_path, force=true)
        rm(genome_path * ".fai", force=true)
        rm(guide_path, force=true)
        rm(output_path, force=true)
    end
    
    # Clean up before test
    cleanup()

    # Store genome to file with proper line breaks
    len = length(genome_seq)
    open(genome_path, "w") do f
        write(f, ">chr1\n")
        write(f, genome_seq)
        write(f, "\n")
    end
    # Create proper .fai index
    open(genome_path * ".fai", "w") do f
        # NAME LENGTH OFFSET LINEBASES LINEWIDTH
        # offset = 6 (length of ">chr1\n")
        write(f, "chr1\t$len\t6\t$len\t$(len + 1)\n")
    end
    
    # Store guide to file
    guide_seq = "TTTTTTTTTTTTTTTTTTTT"
    open(guide_path, "w") do f
        write(f, guide_seq * "\n")
    end
    
    try
        # Test with distance=4 (standard CHOPOFF distance)
        # With 4 leading A's, the match at position 27 has score=4 (4 insertions)
        # So with k=4 we should find it
        motif = Motif("Cas9_Test", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 4, true, 0)
        
        guides = [LongDNA{4}(guide_seq)]
        
        @info "Running Sassy Search with distance=4..."
        CHOPOFF.search_sassy(guides, genome_path, motif, output_path; distance=4)
        
        # Check results
        @test isfile(output_path)
        df = CSV.read(output_path, DataFrame)
        @info "Results:" df
        @test nrow(df) > 0
        
        first_row = df[1, :]
        @test first_row.distance <= 4
        @test first_row.chromosome == "chr1"
        
        @info "Test passed: Match found with distance=$(first_row.distance)"

    finally
        cleanup()
    end

end
