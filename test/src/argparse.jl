using Test

using CRISPRofftargetHunter
using CSV
using DataFrames
using BioSequences


@testset "ArgParse" begin
    @testset "hashDB through command line" begin
        tdir = tempname()
        mkpath(tdir)
        genome = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..", "test", "sample_data", "genome", "semirandom.fa")
        guides_s = "./sample_data/crispritz_results/guides.txt"
        guides = LongDNA{4}.(Set(readlines(guides_s)))

        args = ["build", "--name", "test_hash", "--genome", genome, "--output", tdir, "--motif", "Cas9", 
            "--distance", "1", "--ambig_max", "0", "hashDB"]
        CRISPRofftargetHunter.main(args)
        
        res_file = joinpath(tdir, "hashDB_results.csv")
        args = ["search", tdir, "hashDB", guides_s, res_file, "--right"]
        CRISPRofftargetHunter.main(args)

        # compare the results file with the local results
        hdb_res = search_hashDB(tdir, guides, true)
        res = DataFrame(CSV.File(res_file))
        @test nrow(res) == length(guides)
        @test all(res.guide .== String.(guides))
        @test all(Matrix(res[:, 1:2]) == Matrix(hdb_res[:, 1:2]))
    end
end

