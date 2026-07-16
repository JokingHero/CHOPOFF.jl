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

        cas12a_guide = "TGCATGCATGCATGCATGCAT"
        cas12a_guides_path = joinpath(tdir, "cas12a_guides.txt")
        write(cas12a_guides_path, cas12a_guide * "\n")
        cas12a_genome = joinpath(tdir, "cas12a.fa")
        cas12a_seq = repeat("A", 40) * "TTTA" * cas12a_guide * repeat("A", 40)
        open(cas12a_genome, "w") do io
            write(io, ">chr1\n", cas12a_seq, "\n")
        end
        open(cas12a_genome * ".fai", "w") do io
            write(io, "chr1\t", string(length(cas12a_seq)), "\t6\t",
                string(length(cas12a_seq)), "\t",
                string(length(cas12a_seq) + 1), "\n")
        end
        cas12a_expected = joinpath(tdir, "cas12a_expected.csv")
        cas12a_actual = joinpath(tdir, "cas12a_actual.csv")
        search_prefixHashScan(
            [LongDNA{4}(cas12a_guide)], cas12a_genome, cas12a_expected;
            motif = "Cas12a")
        cas12a_args = [
            "search", "--guides", cas12a_guides_path,
            "--output", cas12a_actual, "prefixHashScan",
            "--genome", cas12a_genome, "--motif", "Cas12a",
        ]
        parsed_cas12a = CHOPOFF.parse_commandline(cas12a_args)
        @test parsed_cas12a["search"]["prefixHashScan"]["motif"] == "Cas12a"
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(cas12a_args)
        @test read(cas12a_actual) == read(cas12a_expected)

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
