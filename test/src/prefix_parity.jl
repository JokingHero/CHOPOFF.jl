using Test
using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "helpers", "prefix_parity.jl"))
using .PrefixParity

function write_parity_fixture(path::String, rows)
    CSV.write(path, DataFrame(rows))
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
end
