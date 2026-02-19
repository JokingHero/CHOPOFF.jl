#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames

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

function get_arg(name::String, default::String)
    flag = "--$name"
    for (i, a) in enumerate(ARGS)
        if a == flag
            if i == length(ARGS)
                error("Missing value for $flag")
            end
            return ARGS[i + 1]
        elseif startswith(a, flag * "=")
            return split(a, "="; limit = 2)[2]
        end
    end
    return default
end

function load_fixture(motif_name::String)
    data_root = joinpath(dirname(pathof(CHOPOFF)), "..", "test", "sample_data")
    if motif_name == "Cas9"
        genome = joinpath(data_root, "genome", "semirandom.fa")
        guides = sort(collect(Set(readlines(joinpath(data_root, "guides.txt")))))
        return genome, LongDNA{4}.(guides)
    elseif motif_name == "Cas12a"
        genome = joinpath(data_root, "genome", "semirandom.2bit")
        guides = LongDNA{4}.([
            "TCGATTGTTTGGCTCTCTAAA",
            "GCAGGGGGACGCAAGTACGAA",
            "GGGCCGAAACGCGACACCGCC",
        ])
        return genome, guides
    else
        error("Unsupported motif '$motif_name'. Use Cas9 or Cas12a.")
    end
end

function print_related(title::String, df::DataFrame, row, max_rows::Int)
    rel = filter(
        x -> x.guide == row.guide && x.chromosome == row.chromosome && x.strand == row.strand,
        df,
    )
    if nrow(rel) > 0
        sort!(rel, [:start, :distance, :alignment_guide, :alignment_reference])
    end
    println(title, " (", nrow(rel), " rows)")
    println(first(rel, min(nrow(rel), max_rows)))
end

function main()
    motif_name = get_arg("motif", "Cas9")
    distance = parse(Int, get_arg("distance", "1"))
    max_rows = parse(Int, get_arg("max-rows", "20"))
    if !(1 <= distance <= 3)
        error("This script expects --distance in 1:3 to mirror verify_sassy_core.")
    end

    motif = Motif(motif_name)
    genome, guides = load_fixture(motif_name)

    tdir = mktempdir()
    ldb_path = joinpath(tdir, "sassy_linearDB")
    mkpath(ldb_path)
    build_linearDB("semirandom_sassy_debug_$motif_name", genome, motif, ldb_path, 7)
    detail_path = joinpath(ldb_path, "detail.csv")

    sassy_path = joinpath(tdir, "sassy_debug_$motif_name.csv")
    search_sassy(guides, genome, motif, sassy_path; distance = distance)
    sassy_df = DataFrame(CSV.File(sassy_path))

    search_linearDB(ldb_path, guides, detail_path; distance = distance)
    ldb_df = DataFrame(CSV.File(detail_path))

    failed_sassy = antijoin(sassy_df, ldb_df, on = REQUIRED_PARITY_COLS)
    failed_linear = antijoin(ldb_df, sassy_df, on = REQUIRED_PARITY_COLS)
    path_sassy = antijoin(sassy_df, ldb_df, on = FULL_PARITY_COLS)
    path_linear = antijoin(ldb_df, sassy_df, on = FULL_PARITY_COLS)

    println("motif=$motif_name distance=$distance")
    println("sassy rows: ", nrow(sassy_df))
    println("linearDB rows: ", nrow(ldb_df))
    println("sassy-only mismatches (required cols): ", nrow(failed_sassy))
    println("linearDB-only mismatches (required cols): ", nrow(failed_linear))
    println("sassy-only mismatches (full incl. path): ", nrow(path_sassy))
    println("linearDB-only mismatches (full incl. path): ", nrow(path_linear))

    if nrow(failed_sassy) == 0 && nrow(failed_linear) == 0
        println("No mismatches found.")
        return
    end

    if nrow(failed_linear) > 0
        row = failed_linear[1, :]
        println("\nFirst mismatch source: linearDB -> sassy")
        println(DataFrame(row))
        print_related("Related linearDB rows", ldb_df, row, max_rows)
        print_related("Related sassy rows", sassy_df, row, max_rows)
    else
        row = failed_sassy[1, :]
        println("\nFirst mismatch source: sassy -> linearDB")
        println(DataFrame(row))
        print_related("Related sassy rows", sassy_df, row, max_rows)
        print_related("Related linearDB rows", ldb_df, row, max_rows)
    end
end

main()
