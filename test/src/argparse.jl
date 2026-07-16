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

    @testset "prefixHashScan direct search" begin
        tdir = tempname()
        mkpath(tdir)
        root = joinpath(dirname(pathof(CHOPOFF)), "..")
        genome = joinpath(root, "test", "sample_data", "genome", "semirandom.fa")
        guides_path = joinpath(root, "test", "sample_data", "guides.txt")
        guides = LongDNA{4}.(readlines(guides_path))
        expected = joinpath(tdir, "expected.csv")
        actual = joinpath(tdir, "actual.csv")
        search_prefixHashScan(guides, genome, expected)

        args = [
            "search", "--guides", guides_path, "--output", actual,
            "prefixHashScan", "--genome", genome,
        ]
        parsed = CHOPOFF.parse_commandline(args)
        @test parsed["search"]["database"] === nothing
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(args)
        @test read(actual) == read(expected)

        sassy_args = [
            "search", "--guides", guides_path, "--output", actual,
            "sassy", "--genome", genome, "--motif", "Cas9",
        ]
        @test CHOPOFF.parse_commandline(sassy_args)["search"]["database"] === nothing

        missing_database = [
            "search", "--guides", guides_path, "--output", actual,
            "prefixHashDB",
        ]
        @test_throws ErrorException CHOPOFF.main(missing_database)

        wrong_distance = [
            "search", "--distance", "2", "--guides", guides_path,
            "--output", actual, "prefixHashScan", "--genome", genome,
        ]
        @test_throws ErrorException CHOPOFF.main(wrong_distance)
    end
end
