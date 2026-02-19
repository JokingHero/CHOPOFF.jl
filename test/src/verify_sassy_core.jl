using Test


using CHOPOFF
using BioSequences
using CSV
using DataFrames

## SET WD when debugging
cd(dirname(dirname(@__FILE__)))

const REQUIRED_PARITY_COLS = [
    :guide,
    :distance,
    :chromosome,
    :start,
    :strand,
]

const FULL_PARITY_COLS = [
    :guide,
    :alignment_guide,
    :alignment_reference,
    :distance,
    :chromosome,
    :start,
    :strand,
]

function print_mismatch_context(
    prefix::String,
    failed::DataFrame,
    source_df::DataFrame,
    target_df::DataFrame,
)
    nrow(failed) == 0 && return
    row = failed[1, :]
    println("$prefix First mismatch:")
    println(DataFrame(row))

    rel_src = filter(
        x -> x.guide == row.guide && x.chromosome == row.chromosome && x.strand == row.strand,
        source_df,
    )
    rel_tgt = filter(
        x -> x.guide == row.guide && x.chromosome == row.chromosome && x.strand == row.strand,
        target_df,
    )
    if nrow(rel_src) > 0
        sort!(rel_src, [:start, :distance, :alignment_guide, :alignment_reference])
    end
    if nrow(rel_tgt) > 0
        sort!(rel_tgt, [:start, :distance, :alignment_guide, :alignment_reference])
    end

    println("$prefix Related source rows:")
    println(first(rel_src, min(nrow(rel_src), 12)))
    println("$prefix Related target rows:")
    println(first(rel_tgt, min(nrow(rel_tgt), 12)))
end

function verify_lineardb_vs_sassy(
    guides::Vector{LongDNA{4}},
    genome::String,
    motif::Motif,
    tdir::String,
    label::String,
)
    @testset "linearDB vs sassy ($label)" begin
        sassy_ldb_path = joinpath(tdir, "sassy_linearDB_$label")
        mkpath(sassy_ldb_path)
        build_linearDB("semirandom_sassy_$label", genome, motif, sassy_ldb_path, 7)
        sassy_detail_path = joinpath(sassy_ldb_path, "detail.csv")

        for d in 1:3
            sassy_path = joinpath(tdir, "sassy_$(label)_d$d.csv")
            search_sassy(guides, genome, motif, sassy_path; distance = d)
            sassy_df = DataFrame(CSV.File(sassy_path))

            search_linearDB(sassy_ldb_path, guides, sassy_detail_path; distance = d)
            ldb_df = DataFrame(CSV.File(sassy_detail_path))

            failed = antijoin(sassy_df, ldb_df, on = REQUIRED_PARITY_COLS)
            if nrow(failed) > 0
                println("Sassy matches not in LinearDB ($label, d=$d): n=$(nrow(failed))")
                println(first(failed, min(40, nrow(failed))))
                print_mismatch_context("Sassy -> LinearDB ($label, d=$d).", failed, sassy_df, ldb_df)
            end
            @test nrow(failed) == 0

            failed2 = antijoin(ldb_df, sassy_df, on = REQUIRED_PARITY_COLS)
            if nrow(failed2) > 0
                println("LinearDB matches not in Sassy ($label, d=$d): n=$(nrow(failed2))")
                println(first(failed2, min(40, nrow(failed2))))
                print_mismatch_context("LinearDB -> Sassy ($label, d=$d).", failed2, ldb_df, sassy_df)
            end
            @test nrow(failed2) == 0

            # Informational only: alignment path strings may differ even when
            # distance/start/strand/chromosome parity is exact.
            path_delta_sassy = antijoin(sassy_df, ldb_df, on = FULL_PARITY_COLS)
            path_delta_ldb = antijoin(ldb_df, sassy_df, on = FULL_PARITY_COLS)
            if nrow(path_delta_sassy) > 0 || nrow(path_delta_ldb) > 0
                println(
                    "Alignment-path deltas ($label, d=$d): ",
                    "sassy-only=", nrow(path_delta_sassy),
                    ", linearDB-only=", nrow(path_delta_ldb),
                    " (non-failing diagnostic)"
                )
            end
        end
    end
end

@testset "databases" begin
    data_root = joinpath(dirname(pathof(CHOPOFF)), "..", "test", "sample_data")
    tdir = tempname()
    mkpath(tdir)

    # Cas9 (extends5 = true)
    cas9_genome = joinpath(data_root, "genome", "semirandom.fa")
    cas9_guides = sort(collect(Set(readlines(joinpath(data_root, "guides.txt")))))
    cas9_guides = LongDNA{4}.(cas9_guides)
    verify_lineardb_vs_sassy(cas9_guides, cas9_genome, Motif("Cas9"), tdir, "Cas9")

    # Cas12a (extends5 = false)
    cas12a_genome = joinpath(data_root, "genome", "semirandom.2bit")
    cas12a_guides = LongDNA{4}.([
        "TCGATTGTTTGGCTCTCTAAA",
        "GCAGGGGGACGCAAGTACGAA",
        "GGGCCGAAACGCGACACCGCC",
    ])
    verify_lineardb_vs_sassy(cas12a_guides, cas12a_genome, Motif("Cas12a"), tdir, "Cas12a")
end
