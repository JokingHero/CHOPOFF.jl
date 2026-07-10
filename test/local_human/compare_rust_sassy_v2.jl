#!/usr/bin/env julia

using CSV
using DataFrames

include("run_human_sassy.jl")

function env_required(name::String)
    val = strip(get(ENV, name, ""))
    isempty(val) && error("Missing required environment variable $name")
    return String(val)
end

function core_for_guides(path::String, guides::Set{String})
    core = normalize_core_parity(path)
    isempty(guides) && return core
    return core[in.(core.guide, Ref(guides)), :]
end

function main()
    rust_path = env_required("CHOPOFF_RUST_SASSY_CSV")
    prefix_path = String(strip(get(ENV, "CHOPOFF_PREFIXHASH_CSV", "")))
    validated_prefix_path = String(strip(get(ENV, "CHOPOFF_VALIDATED_PREFIXHASH_CSV", "")))
    output_dir = String(strip(get(ENV, "CHOPOFF_RUST_PARITY_OUT", dirname(String(rust_path)))))
    genome = String(strip(get(ENV, "CHOPOFF_HUMAN_GENOME", DEFAULT_GENOME)))
    motif_name = String(strip(get(ENV, "CHOPOFF_HUMAN_MOTIF", "Cas9")))
    distance = parse_int_env("CHOPOFF_HUMAN_DISTANCE", 3)
    mkpath(output_dir)

    motif = Motif(motif_name; distance = distance)
    if isempty(validated_prefix_path)
        isempty(prefix_path) && error("Set CHOPOFF_VALIDATED_PREFIXHASH_CSV or CHOPOFF_PREFIXHASH_CSV")
        validated_prefix_path, validated_rows, rejected_rows =
            write_validated_prefix_oracle(prefix_path, genome, motif, distance, output_dir)
    else
        validated_rows = count_result_rows(validated_prefix_path)
        rejected_rows = missing
    end

    rust_core_all = normalize_core_parity(rust_path)
    rust_guides = Set(String.(rust_core_all.guide))
    prefix_core = core_for_guides(validated_prefix_path, rust_guides)
    rust_core = core_for_guides(rust_path, rust_guides)

    rust_only = antijoin(rust_core, prefix_core, on = CORE_PARITY_COLS)
    prefix_only = antijoin(prefix_core, rust_core, on = CORE_PARITY_COLS)
    sort!(rust_only, CORE_PARITY_COLS)
    sort!(prefix_only, CORE_PARITY_COLS)

    rust_only_path = joinpath(output_dir, "parity_rust_only.csv")
    prefix_only_path = joinpath(output_dir, "parity_validated_prefix_only.csv")
    CSV.write(rust_only_path, rust_only)
    CSV.write(prefix_only_path, prefix_only)

    summary = DataFrame([(
        rust_csv = rust_path,
        validated_prefix_csv = validated_prefix_path,
        rust_rows = nrow(rust_core),
        validated_prefix_rows = nrow(prefix_core),
        raw_validated_prefix_rows = validated_rows,
        rejected_prefix_rows = rejected_rows,
        rust_only = nrow(rust_only),
        validated_prefix_only = nrow(prefix_only),
        output_dir = output_dir,
    )])
    CSV.write(joinpath(output_dir, "rust_vs_validated_prefix_summary.csv"), summary)
    println(summary)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
