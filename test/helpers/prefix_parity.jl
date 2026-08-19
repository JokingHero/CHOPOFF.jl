module PrefixParity

using CSV
using DataFrames

export CORE_PARITY_COLS,
    DETAIL_PARITY_COLS,
    PARITY_REASON_COLS,
    classify_semantic_parity,
    expected_counts_from_detail,
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
const PARITY_REASON_COLS = vcat(
    [:side], DETAIL_PARITY_COLS,
    [:ambiguity_count, :valid, :optimal, :reason, :accepted])

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

function detail_value_key(row)
    return Tuple(getproperty(row, column) for column in DETAIL_PARITY_COLS)
end

function core_value_key(row)
    return Tuple(getproperty(row, column) for column in CORE_PARITY_COLS)
end

function semantic_difference_row(side::String, row, inspected, reason, accepted)
    values = NamedTuple{Tuple(DETAIL_PARITY_COLS)}(
        Tuple(getproperty(row, column) for column in DETAIL_PARITY_COLS))
    return merge((side = side,), values, (
        ambiguity_count = Int(inspected.ambiguity_count),
        valid = Bool(inspected.valid),
        optimal = Bool(inspected.optimal),
        reason = String(reason),
        accepted = Bool(accepted),
    ))
end

"""
    classify_semantic_parity(scan_path, prefix_path, requested_ambig, inspect)

Compare detail-row multisets. `inspect(row)` must return
`(ambiguity_count, valid, optimal)`. Exact unambiguous parity remains strict;
only validated ambiguity-specific differences are accepted.
"""
function classify_semantic_parity(
    scan_path::String,
    prefix_path::String,
    requested_ambig::Int,
    inspect::Function;
    output_path::Union{Nothing, String} = nothing)

    scan = normalize_detail_parity(scan_path)
    prefix = normalize_detail_parity(prefix_path)
    join_columns = vcat(DETAIL_PARITY_COLS, :occurrence)
    scan_only = antijoin(scan, prefix; on = join_columns)
    prefix_only = antijoin(prefix, scan; on = join_columns)
    scan_values = Set(detail_value_key(row) for row in eachrow(scan))
    prefix_values = Set(detail_value_key(row) for row in eachrow(prefix))

    scan_inspected = [inspect(row) for row in eachrow(scan_only)]
    prefix_inspected = [inspect(row) for row in eachrow(prefix_only)]
    scan_reason = fill("", nrow(scan_only))
    prefix_reason = fill("", nrow(prefix_only))
    scan_accepted = falses(nrow(scan_only))
    prefix_accepted = falses(nrow(prefix_only))

    # Pair different, equally optimal traceback strings at the same semantic hit.
    prefix_by_core = Dict{Tuple, Vector{Int}}()
    for (idx, row) in enumerate(eachrow(prefix_only))
        push!(get!(prefix_by_core, core_value_key(row), Int[]), idx)
    end
    used_prefix = falses(nrow(prefix_only))
    for (scan_idx, row) in enumerate(eachrow(scan_only))
        inspected = scan_inspected[scan_idx]
        inspected.ambiguity_count > 0 || continue
        inspected.valid && inspected.optimal || continue
        candidates = get(prefix_by_core, core_value_key(row), Int[])
        prefix_idx = findfirst(idx -> !used_prefix[idx] &&
            prefix_inspected[idx].valid && prefix_inspected[idx].optimal,
            candidates)
        prefix_idx === nothing && continue
        matched = candidates[prefix_idx]
        used_prefix[matched] = true
        scan_reason[scan_idx] = "optimal_traceback_tie"
        prefix_reason[matched] = "optimal_traceback_tie"
        scan_accepted[scan_idx] = true
        prefix_accepted[matched] = true
    end

    for (idx, row) in enumerate(eachrow(scan_only))
        isempty(scan_reason[idx]) || continue
        inspected = scan_inspected[idx]
        if detail_value_key(row) in prefix_values
            scan_reason[idx] = "scan_duplicate"
        elseif inspected.ambiguity_count == 0
            scan_reason[idx] = "unambiguous_difference"
        elseif !inspected.valid || !inspected.optimal ||
                inspected.ambiguity_count > requested_ambig
            scan_reason[idx] = "invalid_scan_row"
        else
            scan_reason[idx] = "legacy_false_negative"
            scan_accepted[idx] = true
        end
    end

    for (idx, row) in enumerate(eachrow(prefix_only))
        isempty(prefix_reason[idx]) || continue
        inspected = prefix_inspected[idx]
        if detail_value_key(row) in scan_values && inspected.ambiguity_count > 0
            prefix_reason[idx] = "legacy_duplicate"
            prefix_accepted[idx] = true
        elseif inspected.ambiguity_count == 0
            prefix_reason[idx] = "unambiguous_difference"
        elseif inspected.ambiguity_count > requested_ambig
            prefix_reason[idx] = "baseline_above_ambiguity_limit"
            prefix_accepted[idx] = true
        elseif !inspected.valid || !inspected.optimal
            prefix_reason[idx] = "legacy_false_positive"
            prefix_accepted[idx] = true
        else
            prefix_reason[idx] = "scan_miss"
        end
    end

    rows = NamedTuple[]
    for (idx, row) in enumerate(eachrow(scan_only))
        push!(rows, semantic_difference_row(
            "scan", row, scan_inspected[idx], scan_reason[idx], scan_accepted[idx]))
    end
    for (idx, row) in enumerate(eachrow(prefix_only))
        push!(rows, semantic_difference_row(
            "prefix", row, prefix_inspected[idx], prefix_reason[idx], prefix_accepted[idx]))
    end
    differences = isempty(rows) ? DataFrame(
        [column => (column in (:distance, :start, :ambiguity_count) ? Int[] :
            column in (:valid, :optimal, :accepted) ? Bool[] : String[])
         for column in PARITY_REASON_COLS]) : DataFrame(rows)
    output_path === nothing || CSV.write(output_path, differences)
    passed = all(differences.accepted)
    return (
        passed = passed,
        exact = isempty(rows),
        scan_only = nrow(scan_only),
        prefix_only = nrow(prefix_only),
        accepted = count(differences.accepted),
        rejected = count(!, differences.accepted),
        differences = differences,
    )
end

function expected_counts_from_detail(
    detail_path::String, guides::Vector{String}, distance::Int)

    ordered_guides = unique(guides)
    counts = zeros(Int, length(ordered_guides), distance + 1)
    row_by_guide = Dict(guide => idx for (idx, guide) in enumerate(ordered_guides))
    detail = normalize_core_parity(detail_path)
    unique!(detail, CORE_PARITY_COLS)
    for row in eachrow(detail)
        row.distance <= distance || continue
        counts[row_by_guide[row.guide], row.distance + 1] += 1
    end
    out = DataFrame(guide = ordered_guides)
    for dist in 0:distance
        out[!, Symbol("D$dist")] = counts[:, dist + 1]
    end
    out.complete = fill(true, length(ordered_guides))
    return out
end

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
