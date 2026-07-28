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
    prefix_path = env_required("CHOPOFF_PREFIXHASH_CSV")
    output_dir = String(strip(get(ENV, "CHOPOFF_RUST_PARITY_OUT", dirname(String(rust_path)))))
    mkpath(output_dir)

    rust_core_all = normalize_core_parity(rust_path)
    rust_guides = Set(String.(rust_core_all.guide))
    prefix_core = core_for_guides(prefix_path, rust_guides)
    rust_core = core_for_guides(rust_path, rust_guides)

    join_columns = vcat(CORE_PARITY_COLS, :occurrence)
    rust_only = antijoin(rust_core, prefix_core; on = join_columns)
    prefix_only = antijoin(prefix_core, rust_core; on = join_columns)
    sort!(rust_only, CORE_PARITY_COLS)
    sort!(prefix_only, CORE_PARITY_COLS)
    select!(rust_only, Not(:occurrence))
    select!(prefix_only, Not(:occurrence))

    rust_only_path = joinpath(output_dir, "parity_rust_only.csv")
    prefix_only_path = joinpath(output_dir, "parity_prefix_only.csv")
    CSV.write(rust_only_path, rust_only)
    CSV.write(prefix_only_path, prefix_only)

    summary = DataFrame([(
        rust_csv = rust_path,
        prefix_csv = prefix_path,
        rust_rows = nrow(rust_core),
        prefix_rows = nrow(prefix_core),
        rust_only = nrow(rust_only),
        prefix_only = nrow(prefix_only),
        output_dir = output_dir,
    )])
    CSV.write(joinpath(output_dir, "rust_vs_prefix_summary.csv"), summary)
    println(summary)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
