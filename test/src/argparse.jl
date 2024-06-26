using Test

using CHOPOFF
using CSV
using DataFrames
using BioSequences

## SET WD when debugging
# cd("test")

@testset "ArgParse" begin
    @testset "hashDB through command line" begin
        tdir = tempname()
        mkpath(tdir)
        genome = joinpath(dirname(pathof(CHOPOFF)), "..", "test", "sample_data", "genome", "semirandom.fa")
        guides_s = "./sample_data/guides.txt"
        guides = LongDNA{4}.(Set(readlines(guides_s)))

        tdirDB = joinpath(tdir, "hashDB.bin")
        args = ["build", "--name", "test_hash", "--genome", genome, "--output", tdirDB, "--motif", "Cas9", 
            "--distance", "1", "--ambig_max", "0", "hashDB"]
        CHOPOFF.main(args)
        
        res_file = joinpath(tdir, "hashDB_results.csv")
        args = ["estimate", "--database", tdirDB, "--guides", guides_s, "--output", res_file, "--right"]
        CHOPOFF.main(args)

        # compare the results file with the local results
        db = load(tdirDB)
        hdb_res = search_hashDB(db, guides, true)
        res = DataFrame(CSV.File(res_file))
        @test nrow(res) == length(guides)
        @test all(res.guide .== String.(guides))
        @test all(Matrix(res[:, 2:3]) == Matrix(hdb_res[:, 2:3]))
    end
end

