module PrefixParity

using CSV
using DataFrames

export CORE_PARITY_COLS,
    DETAIL_PARITY_COLS,
    normalize_core_parity,
    normalize_detail_parity,
    write_exact_parity_diffs,
    write_named_parity_diffs,
    write_parity_diffs

const CORE_PARITY_COLS =
    [:guide, :distance, :chromosome, :start, :strand]
const DETAIL_PARITY_COLS = [
    :guide,
    :alignment_guide,
    :alignment_reference,
    :distance,
    :chromosome,
    :start,
    :strand,
]

function normalize_parity(path::String, columns::Vector{Symbol})
    df = DataFrame(CSV.File(path))
    available = Set(Symbol.(names(df)))
    for column in columns
        column in available ||
            error("Missing required parity column '$column' in $path")
    end
    normalized = select(df, columns)
    for column in columns
        if column in (:distance, :start)
            normalized[!, column] =
                Int[Int(value) for value in normalized[!, column]]
        else
            normalized[!, column] =
                String[String(value) for value in normalized[!, column]]
        end
    end
    sort!(normalized, columns)

    occurrence = Vector{Int}(undef, nrow(normalized))
    current = 0
    for row_idx in 1:nrow(normalized)
        if row_idx == 1 || any(
                normalized[row_idx, column] !=
                    normalized[row_idx - 1, column]
                for column in columns)
            current = 1
        else
            current += 1
        end
        occurrence[row_idx] = current
    end
    normalized.occurrence = occurrence
    return normalized
end

normalize_core_parity(path::String) =
    normalize_parity(path, CORE_PARITY_COLS)

normalize_detail_parity(path::String) =
    normalize_parity(path, DETAIL_PARITY_COLS)

function write_parity_diffs(
    lhs_path::String,
    rhs_path::String,
    output_dir::String,
    lhs_name::String,
    rhs_name::String,
    columns::Vector{Symbol})

    lhs = normalize_parity(lhs_path, columns)
    rhs = normalize_parity(rhs_path, columns)
    join_columns = vcat(columns, :occurrence)
    lhs_only = antijoin(lhs, rhs; on = join_columns)
    rhs_only = antijoin(rhs, lhs; on = join_columns)
    select!(lhs_only, Not(:occurrence))
    select!(rhs_only, Not(:occurrence))
    CSV.write(
        joinpath(output_dir, "parity_$(lhs_name)_only.csv"), lhs_only)
    CSV.write(
        joinpath(output_dir, "parity_$(rhs_name)_only.csv"), rhs_only)
    return nrow(lhs_only), nrow(rhs_only)
end

function write_exact_parity_diffs(
    scan_path::String,
    prefix_path::String,
    output_dir::String)

    return write_parity_diffs(
        scan_path,
        prefix_path,
        output_dir,
        "scan",
        "prefix",
        DETAIL_PARITY_COLS)
end

function write_named_parity_diffs(
    lhs_path::String,
    rhs_path::String,
    output_dir::String,
    lhs_name::String,
    rhs_name::String)

    return write_parity_diffs(
        lhs_path,
        rhs_path,
        output_dir,
        lhs_name,
        rhs_name,
        CORE_PARITY_COLS)
end

function write_parity_diffs(
    sassy_path::String,
    prefix_path::String,
    output_dir::String)

    return write_named_parity_diffs(
        sassy_path, prefix_path, output_dir, "sassy", "prefix")
end

end
