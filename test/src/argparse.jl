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

        twobit_genome = joinpath(
            root, "test", "sample_data", "genome", "semirandom.2bit")
        twobit_actual = joinpath(tdir, "actual_2bit.csv")
        twobit_args = [
            "search", "--guides", guides_path, "--output", twobit_actual,
            "prefixHashScan", "--genome", twobit_genome,
        ]
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(twobit_args)
        @test read(twobit_actual) == read(expected)

        large_guides_path = joinpath(tdir, "large_guides.txt")
        write(
            large_guides_path,
            join(fill(string(first(guides)), 65), "\n") * "\n")
        large_expected = joinpath(tdir, "large_expected.csv")
        large_actual = joinpath(tdir, "large_actual.csv")
        search_prefixHashScan(
            fill(first(guides), 65), genome, large_expected)
        large_args = [
            "search", "--guides", large_guides_path,
            "--output", large_actual, "prefixHashScan", "--genome", genome,
        ]
        @test_logs (:info, r"prefixHashScan guide batching") (:info, r"prefixHashScan execution") CHOPOFF.main(
            large_args)
        @test read(large_actual) == read(large_expected)

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


        cas12a_d1_expected = joinpath(tdir, "cas12a_d1_expected.csv")
        cas12a_d1_actual = joinpath(tdir, "cas12a_d1_actual.csv")
        search_prefixHashScan(
            [LongDNA{4}(cas12a_guide)], cas12a_genome, cas12a_d1_expected;
            motif = "Cas12a", distance = 1)
        cas12a_d1_args = [
            "search", "--distance", "1", "--guides", cas12a_guides_path,
            "--output", cas12a_d1_actual, "prefixHashScan",
            "--genome", cas12a_genome, "--motif", "Cas12a",
        ]
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(
            cas12a_d1_args)
        @test read(cas12a_d1_actual) == read(cas12a_d1_expected)

        generic_expected = joinpath(tdir, "generic_expected.csv")
        generic_actual = joinpath(tdir, "generic_actual.csv")
        search_prefixHashScan(
            guides, genome, generic_expected;
            motif = "Cas9_NGA", distance = 1)
        generic_args = [
            "search", "--distance", "1", "--guides", guides_path,
            "--output", generic_actual, "prefixHashScan",
            "--genome", genome, "--motif", "Cas9_NGA",
        ]
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(
            generic_args)
        @test read(generic_actual) == read(generic_expected)

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

        distance_two_expected = joinpath(tdir, "distance_two_expected.csv")
        distance_two_actual = joinpath(tdir, "distance_two_actual.csv")
        search_prefixHashScan(
            guides, genome, distance_two_expected; distance = 2)
        distance_two = [
            "search", "--distance", "2", "--guides", guides_path,
            "--output", distance_two_actual,
            "prefixHashScan", "--genome", genome,
        ]
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(
            distance_two)
        @test read(distance_two_actual) == read(distance_two_expected)

        distance_four_expected = joinpath(tdir, "distance_four_expected.csv")
        distance_four_actual = joinpath(tdir, "distance_four_actual.csv")
        search_prefixHashScan(
            guides, genome, distance_four_expected; distance = 4)
        distance_four = [
            "search", "--distance", "4", "--guides", guides_path,
            "--output", distance_four_actual,
            "prefixHashScan", "--genome", genome,
        ]
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(
            distance_four)
        @test read(distance_four_actual) == read(distance_four_expected)

        ambiguous_guide = "ACGTACGTACGTACGTACGT"
        ambiguous_guides_path = joinpath(tdir, "ambiguous_guides.txt")
        write(ambiguous_guides_path, ambiguous_guide * "\n")
        ambiguous_genome = joinpath(tdir, "ambiguous.fa")
        ambiguous_seq = repeat("A", 40) *
            ambiguous_guide[1:9] * "N" * ambiguous_guide[11:end] *
            "AGG" * repeat("A", 40)
        open(ambiguous_genome, "w") do io
            write(io, ">chr1\n", ambiguous_seq, "\n")
        end
        open(ambiguous_genome * ".fai", "w") do io
            write(io, "chr1\t", string(length(ambiguous_seq)), "\t6\t",
                string(length(ambiguous_seq)), "\t",
                string(length(ambiguous_seq) + 1), "\n")
        end
        ambiguous_expected = joinpath(tdir, "ambiguous_expected.csv")
        ambiguous_actual = joinpath(tdir, "ambiguous_actual.csv")
        search_prefixHashScan(
            [LongDNA{4}(ambiguous_guide)], ambiguous_genome,
            ambiguous_expected;
            motif = Motif("Cas9"; distance = 0, ambig_max = 1),
            distance = 0)
        ambiguous_args = [
            "search", "--distance", "0", "--guides", ambiguous_guides_path,
            "--output", ambiguous_actual, "prefixHashScan",
            "--genome", ambiguous_genome, "--ambig_max", "1",
        ]
        @test CHOPOFF.parse_commandline(ambiguous_args)["search"]["prefixHashScan"]["ambig_max"] == 1
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(
            ambiguous_args)
        @test read(ambiguous_actual) == read(ambiguous_expected)
        invalid_ambiguous_args = copy(ambiguous_args)
        invalid_ambiguous_args[end] = "4"
        @test_throws ErrorException CHOPOFF.main(invalid_ambiguous_args)

        wrong_distance = [
            "search", "--distance", "5", "--guides", guides_path,
            "--output", actual, "prefixHashScan", "--genome", genome,
        ]
        @test_throws ErrorException CHOPOFF.main(wrong_distance)
    end
end
