using Test
using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "helpers", "prefix_parity.jl"))
using .PrefixParity

function write_parity_fixture(path::String, rows)
    CSV.write(path, DataFrame(rows))
    return path
end

function empty_parity_fixture(path::String)
    CSV.write(path, DataFrame(
        guide = String[], alignment_guide = String[],
        alignment_reference = String[], distance = Int[],
        chromosome = String[], start = Int[], strand = String[]))
    return path
end

@testset "prefixHashDB parity helpers" begin
    tdir = mktempdir()
    row1 = (
        guide = "AAAAAAAAAAAAAAAAAAAA",
        alignment_guide = "AAAAAAAAAAAAAAAAAAAA",
        alignment_reference = "AAAAAAAAAAAAAAAAAAAA",
        distance = 0,
        chromosome = "chr1",
        start = 10,
        strand = "+",
    )
    row2 = (
        guide = "CCCCCCCCCCCCCCCCCCCC",
        alignment_guide = "CCCCCCCCCCCCCCCCCCCC",
        alignment_reference = "CCCCCCCCCCCCCCCCCCCT",
        distance = 1,
        chromosome = "chr2",
        start = 20,
        strand = "-",
    )

    @testset "row order is ignored" begin
        scan = write_parity_fixture(
            joinpath(tdir, "ordered_scan.csv"), [row1, row2])
        prefix = write_parity_fixture(
            joinpath(tdir, "ordered_prefix.csv"), [row2, row1])
        scan_only, prefix_only =
            write_exact_parity_diffs(scan, prefix, tdir)
        @test scan_only == 0
        @test prefix_only == 0
    end

    @testset "detail differences fail exact parity" begin
        changed = merge(row1, (alignment_reference = "TAAAAAAAAAAAAAAAAAAA",))
        scan = write_parity_fixture(
            joinpath(tdir, "detail_scan.csv"), [changed])
        prefix = write_parity_fixture(
            joinpath(tdir, "detail_prefix.csv"), [row1])
        scan_only, prefix_only =
            write_exact_parity_diffs(scan, prefix, tdir)
        @test scan_only == 1
        @test prefix_only == 1

        scan_only, prefix_only = write_named_parity_diffs(
            scan, prefix, tdir, "sassy", "prefix")
        @test scan_only == 0
        @test prefix_only == 0
    end

    @testset "duplicate multiplicity is exact" begin
        scan = write_parity_fixture(
            joinpath(tdir, "duplicate_scan.csv"), [row1, row1])
        prefix = write_parity_fixture(
            joinpath(tdir, "duplicate_prefix.csv"), [row1])
        scan_only, prefix_only =
            write_exact_parity_diffs(scan, prefix, tdir)
        @test scan_only == 1
        @test prefix_only == 0
        diff = DataFrame(CSV.File(
            joinpath(tdir, "parity_scan_only.csv")))
        @test nrow(diff) == 1
        @test Symbol.(names(diff)) == DETAIL_PARITY_COLS
    end

    @testset "empty outputs retain concrete schema" begin
        empty_rows = DataFrame(
            guide = String[],
            alignment_guide = String[],
            alignment_reference = String[],
            distance = Int[],
            chromosome = String[],
            start = Int[],
            strand = String[],
        )
        scan = joinpath(tdir, "empty_scan.csv")
        prefix = joinpath(tdir, "empty_prefix.csv")
        CSV.write(scan, empty_rows)
        CSV.write(prefix, empty_rows)
        scan_only, prefix_only =
            write_exact_parity_diffs(scan, prefix, tdir)
        @test scan_only == 0
        @test prefix_only == 0
        @test nrow(DataFrame(CSV.File(
            joinpath(tdir, "parity_scan_only.csv")))) == 0
    end


    @testset "semantic ambiguity classifications" begin
        ambiguous = merge(row1, (
            alignment_reference = "NAAAAAAAAAAAAAAAAAAA", start = 30))
        changed_traceback = merge(ambiguous, (
            alignment_guide = "A-AAAAAAAAAAAAAAAAAAA",
            alignment_reference = "NA-AAAAAAAAAAAAAAAAAA"))
        inspect = row -> (
            ambiguity_count = row.start == 30 ? 1 : 0,
            valid = true,
            optimal = true,
        )

        scan = write_parity_fixture(joinpath(tdir, "semantic_scan.csv"), [ambiguous])
        prefix = write_parity_fixture(joinpath(tdir, "semantic_prefix.csv"), [ambiguous])
        result = classify_semantic_parity(scan, prefix, 1, inspect)
        @test result.passed
        @test result.exact

        prefix = write_parity_fixture(
            joinpath(tdir, "semantic_prefix_duplicate.csv"),
            [ambiguous, ambiguous])
        result = classify_semantic_parity(scan, prefix, 1, inspect)
        @test result.passed
        @test result.prefix_only == 1
        @test only(result.differences.reason) == "legacy_duplicate"

        scan_duplicate = write_parity_fixture(
            joinpath(tdir, "semantic_scan_duplicate.csv"),
            [ambiguous, ambiguous])
        prefix_single = write_parity_fixture(
            joinpath(tdir, "semantic_prefix_single.csv"), [ambiguous])
        result = classify_semantic_parity(
            scan_duplicate, prefix_single, 1, inspect)
        @test !result.passed
        @test only(result.differences.reason) == "scan_duplicate"

        prefix = empty_parity_fixture(joinpath(tdir, "semantic_prefix_empty.csv"))
        result = classify_semantic_parity(scan, prefix, 1, inspect)
        @test result.passed
        @test only(result.differences.reason) == "legacy_false_negative"

        result = classify_semantic_parity(prefix, scan, 1, inspect)
        @test !result.passed
        @test only(result.differences.reason) == "scan_miss"

        result = classify_semantic_parity(prefix, scan, 0, inspect)
        @test result.passed
        @test only(result.differences.reason) ==
            "baseline_above_ambiguity_limit"

        prefix = write_parity_fixture(
            joinpath(tdir, "semantic_prefix_tie.csv"), [changed_traceback])
        result = classify_semantic_parity(scan, prefix, 1, inspect)
        @test result.passed
        @test result.scan_only == 1
        @test result.prefix_only == 1
        @test all(==("optimal_traceback_tie"), result.differences.reason)

        invalid_inspect = row -> (
            ambiguity_count = 1, valid = false, optimal = false)
        prefix = empty_parity_fixture(joinpath(tdir, "invalid_prefix.csv"))
        result = classify_semantic_parity(scan, prefix, 1, invalid_inspect)
        @test !result.passed
        @test only(result.differences.reason) == "invalid_scan_row"

        result = classify_semantic_parity(prefix, scan, 1, invalid_inspect)
        @test result.passed
        @test only(result.differences.reason) == "legacy_false_positive"

        unambiguous = write_parity_fixture(
            joinpath(tdir, "unambiguous_scan.csv"), [row1])
        result = classify_semantic_parity(
            unambiguous, prefix, 1,
            row -> (ambiguity_count = 0, valid = true, optimal = true))
        @test !result.passed
        @test only(result.differences.reason) == "unambiguous_difference"
    end

    @testset "semantic detail counts" begin
        detail = write_parity_fixture(
            joinpath(tdir, "semantic_counts_detail.csv"),
            [row1, row1, row2])
        counts = expected_counts_from_detail(
            detail,
            [row1.guide, row2.guide, "GGGGGGGGGGGGGGGGGGGG", row1.guide],
            1)
        @test counts.guide == [
            row1.guide, row2.guide, "GGGGGGGGGGGGGGGGGGGG"]
        @test counts.D0 == [1, 0, 0]
        @test counts.D1 == [0, 1, 0]
        @test all(counts.complete)
    end
end
