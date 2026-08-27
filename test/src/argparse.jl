using Test

using CHOPOFF
using CSV
using DataFrames
using BioSequences

## SET WD when debugging
# cd("test")

function write_argparse_fasta(path::String, sequence::String)
    open(path, "w") do io
        write(io, ">chr1\n", sequence, "\n")
    end
    open(path * ".fai", "w") do io
        write(io, "chr1\t", string(length(sequence)), "\t6\t",
            string(length(sequence)), "\t", string(length(sequence) + 1), "\n")
    end
end

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
        @test parsed["search"]["prefixHashScan"]["simd_backend"] == "auto"
        forced_simd_args = vcat(args, ["--simd_backend", "avx2"])
        parsed_forced_simd = CHOPOFF.parse_commandline(forced_simd_args)
        @test parsed_forced_simd["search"]["prefixHashScan"]["simd_backend"] ==
            "avx2"
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(args)
        @test read(actual) == read(expected)

        counts_expected = joinpath(tdir, "counts_expected.csv")
        counts_actual = joinpath(tdir, "counts_actual.csv")
        search_prefixHashScan(
            guides, genome, counts_expected; output = :counts)
        counts_args = [
            "search", "--guides", guides_path, "--output", counts_actual,
            "prefixHashScan", "--genome", genome,
            "--output_mode", "counts",
        ]
        parsed_counts = CHOPOFF.parse_commandline(counts_args)
        @test parsed_counts["search"]["prefixHashScan"]["output_mode"] ==
            "counts"
        @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(counts_args)
        @test read(counts_actual) == read(counts_expected)

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

        @testset "custom prefixHashScan motifs" begin
            guide20 = "ACGTACGTACGTACGTACGT"
            guide21 = "TGCATGCATGCATGCATGCAT"
            downstream_motif = Motif(
                "downstream_forward", repeat("N", 20) * "XXX",
                repeat("X", 20) * "NGA", true, false, 0, true, 0)
            upstream_motif = Motif(
                "upstream_extend3", "XXXX" * repeat("N", 21),
                "TTTV" * repeat("X", 21), true, false, 0, false, 0)
            internal_motif = Motif(
                "internal", repeat("N", 10) * "XXX" * repeat("N", 10),
                repeat("X", 10) * "AGG" * repeat("X", 10),
                true, true, 0, true, 0)
            pamless_motif = Motif(
                "pamless", repeat("N", 20), repeat("X", 20),
                true, true, 0, true, 0)
            reverse_motif = Motif(
                "downstream_reverse", repeat("N", 20) * "XXX",
                repeat("X", 20) * "NGA", false, true, 0, true, 0)
            forward_site = guide20 * "AGA"
            cases = [
                (
                    motif = downstream_motif,
                    guide = guide20,
                    site = forward_site,
                    fwd_motif = repeat("N", 20) * "XXX",
                    fwd_pam = repeat("X", 20) * "NGA",
                    flags = ["--not_reverse"],
                    strand = "+",
                ),
                (
                    motif = upstream_motif,
                    guide = guide21,
                    site = "TTTA" * guide21,
                    fwd_motif = "XXXX" * repeat("N", 21),
                    fwd_pam = "TTTV" * repeat("X", 21),
                    flags = ["--not_reverse", "--extend3"],
                    strand = "+",
                ),
                (
                    motif = internal_motif,
                    guide = guide20,
                    site = guide20[1:10] * "AGG" * guide20[11:end],
                    fwd_motif = repeat("N", 10) * "XXX" * repeat("N", 10),
                    fwd_pam = repeat("X", 10) * "AGG" * repeat("X", 10),
                    flags = String[],
                    strand = nothing,
                ),
                (
                    motif = pamless_motif,
                    guide = guide20,
                    site = guide20,
                    fwd_motif = repeat("N", 20),
                    fwd_pam = repeat("X", 20),
                    flags = String[],
                    strand = nothing,
                ),
                (
                    motif = reverse_motif,
                    guide = guide20,
                    site = string(reverse_complement(LongDNA{4}(forward_site))),
                    fwd_motif = repeat("N", 20) * "XXX",
                    fwd_pam = repeat("X", 20) * "NGA",
                    flags = ["--not_forward"],
                    strand = "-",
                ),
            ]

            for case in cases
                case_dir = joinpath(tdir, case.motif.alias)
                mkpath(case_dir)
                case_genome = joinpath(case_dir, "genome.fa")
                case_guides = joinpath(case_dir, "guides.txt")
                expected_output = joinpath(case_dir, "expected.csv")
                actual_output = joinpath(case_dir, "actual.csv")
                write_argparse_fasta(
                    case_genome, repeat("C", 40) * case.site * repeat("C", 40))
                write(case_guides, case.guide * "\n")
                search_prefixHashScan(
                    [LongDNA{4}(case.guide)], case_genome, expected_output;
                    motif = case.motif, distance = 0)
                custom_args = [
                    "search", "--distance", "0", "--guides", case_guides,
                    "--output", actual_output, "prefixHashScan",
                    "--genome", case_genome, "--name", case.motif.alias,
                    "--fwd_motif", case.fwd_motif,
                    "--fwd_pam", case.fwd_pam,
                    case.flags...,
                ]
                @test_logs (:info, r"prefixHashScan execution") CHOPOFF.main(
                    custom_args)
                @test read(actual_output) == read(expected_output)
                result = DataFrame(CSV.File(actual_output))
                @test nrow(result) > 0
                if case.strand !== nothing
                    @test all(==(case.strand), result.strand)
                end
            end

            invalid_base = [
                "search", "--distance", "0", "--guides", guides_path,
                "--output", joinpath(tdir, "invalid_custom.csv"),
                "prefixHashScan", "--genome", genome,
            ]
            @test_throws ErrorException CHOPOFF.main(vcat(
                invalid_base,
                ["--motif", "Cas9", "--fwd_motif", repeat("N", 20) * "XXX",
                 "--fwd_pam", repeat("X", 20) * "NGG"]))
            @test_throws ErrorException CHOPOFF.main(vcat(
                invalid_base, ["--fwd_motif", repeat("N", 20) * "XXX"]))
            @test_throws ErrorException CHOPOFF.main(vcat(
                invalid_base,
                ["--fwd_motif", repeat("N", 20) * "XXX",
                 "--fwd_pam", repeat("X", 20) * "NGG",
                 "--not_forward", "--not_reverse"]))
            @test_throws ErrorException CHOPOFF.main(vcat(
                invalid_base, ["--motif", "unknown_motif"]))
            @test_throws ErrorException CHOPOFF.main(vcat(
                invalid_base,
                ["--fwd_motif", "NNNXNNNX", "--fwd_pam", "XXXAXXXG"]))
        end

        sassy_args = [
            "search", "--guides", guides_path, "--output", actual,
            "sassy", "--genome", genome, "--motif", "Cas9",
            "--backend", "avx2_safe",
        ]
        parsed_sassy = CHOPOFF.parse_commandline(sassy_args)
        @test parsed_sassy["search"]["database"] === nothing
        @test parsed_sassy["search"]["sassy"]["backend"] == "avx2_safe"

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
