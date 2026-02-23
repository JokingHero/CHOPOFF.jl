#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Dates
using Printf
using Random
using Statistics

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const SAMPLE_DIR = joinpath(ROOT_DIR, "test", "sample_data")
const CORE_PARITY_COLS = [:guide, :distance, :chromosome, :start, :strand]
const TRACEBACK_DNA_ALPHABET = ('A', 'C', 'G', 'T')
const CAS12A_GUIDES = LongDNA{4}.([
    "TCGATTGTTTGGCTCTCTAAA",
    "GCAGGGGGACGCAAGTACGAA",
    "GGGCCGAAACGCGACACCGCC",
])

function parse_int_env(name::String, default::Int)
    raw = strip(get(ENV, name, string(default)))
    try
        return parse(Int, raw)
    catch
        error("Invalid integer in environment variable $name: '$raw'")
    end
end

function parse_float_env(name::String, default::Float64)
    raw = strip(get(ENV, name, string(default)))
    try
        return parse(Float64, raw)
    catch
        error("Invalid float in environment variable $name: '$raw'")
    end
end

function parse_bool_env(name::String, default::Bool)
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "y")
        return true
    elseif raw in ("0", "false", "no", "n")
        return false
    end
    error("Invalid boolean in environment variable $name: '$raw'")
end

function parse_distances_env()
    raw = strip(get(ENV, "SASSY_TRACEBACK_DISTANCES", "3"))
    vals = Int[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        d = parse(Int, t)
        (1 <= d <= 8) || error("Distance must be in 1:8, got $d")
        push!(vals, d)
    end
    isempty(vals) && error("SASSY_TRACEBACK_DISTANCES resolved to an empty set.")
    return unique(vals)
end

function parse_strict_modes_env()
    raw = lowercase(strip(get(ENV, "SASSY_TRACEBACK_STRICT_MODES", "strict")))
    if raw == "strict"
        return [true]
    elseif raw == "nonstrict"
        return [false]
    elseif raw == "both"
        return [true, false]
    end
    error("Invalid SASSY_TRACEBACK_STRICT_MODES='$raw'. Use strict, nonstrict, or both.")
end

function require_file(path::String)
    isfile(path) || error("Required file not found: $path")
    return path
end

function load_cases(cas9_limit::Int)
    cas9_genome = require_file(joinpath(SAMPLE_DIR, "genome", "semirandom.fa"))
    require_file(cas9_genome * ".fai")
    cas9_guides_path = require_file(joinpath(SAMPLE_DIR, "guides.txt"))
    cas9_guides_all = LongDNA{4}.(readlines(cas9_guides_path))
    n_cas9 = min(length(cas9_guides_all), cas9_limit)
    n_cas9 > 0 || error("No Cas9 guides available after applying SASSY_TRACEBACK_CAS9_GUIDES.")
    cas9_guides = cas9_guides_all[1:n_cas9]

    cas12a_genome = require_file(joinpath(SAMPLE_DIR, "genome", "semirandom.2bit"))

    return [
        (
            label = "Cas9",
            motif_name = "Cas9",
            genome = cas9_genome,
            guides = cas9_guides,
        ),
        (
            label = "Cas12a",
            motif_name = "Cas12a",
            genome = cas12a_genome,
            guides = CAS12A_GUIDES,
        ),
    ]
end

function normalize_core_parity(df::DataFrame)
    cols = Set(Symbol.(names(df)))
    for c in CORE_PARITY_COLS
        c in cols || error("Missing required parity column '$c'.")
    end
    core = select(df, CORE_PARITY_COLS)
    core.guide = String.(core.guide)
    core.distance = Int.(core.distance)
    core.chromosome = String.(core.chromosome)
    core.start = Int.(core.start)
    core.strand = String.(core.strand)
    sort!(core, CORE_PARITY_COLS)
    return core
end

function check_core_parity(lhs_path::String, rhs_path::String)
    lhs = normalize_core_parity(DataFrame(CSV.File(lhs_path)))
    rhs = normalize_core_parity(DataFrame(CSV.File(rhs_path)))
    if lhs == rhs
        return true, DataFrame(), DataFrame()
    end
    lhs_only = antijoin(lhs, rhs, on = CORE_PARITY_COLS)
    rhs_only = antijoin(rhs, lhs, on = CORE_PARITY_COLS)
    return false, lhs_only, rhs_only
end

@inline function custom_traceback(
    guide_for_aln::LongDNA{4},
    ref_for_aln::LongDNA{4},
    k::Int,
)
    return CHOPOFF.Sassy.traceback_custom(guide_for_aln, ref_for_aln, k)
end

@inline function random_dna(rng::AbstractRNG, len::Int)
    return String([rand(rng, TRACEBACK_DNA_ALPHABET) for _ in 1:len])
end

function mutate_for_micro(seq::String, k::Int, rng::AbstractRNG)
    chars = collect(seq)
    n_sub = rand(rng, 0:k)
    for _ in 1:n_sub
        idx = rand(rng, 1:length(chars))
        old = chars[idx]
        new_base = old
        while new_base == old
            new_base = rand(rng, TRACEBACK_DNA_ALPHABET)
        end
        chars[idx] = new_base
    end
    ext_len = rand(rng, 0:k)
    for _ in 1:ext_len
        push!(chars, rand(rng, TRACEBACK_DNA_ALPHABET))
    end
    return String(chars)
end

function build_micro_inputs(cases, reps_per_k::Int)
    rng = MersenneTwister(20260222)
    inputs = Vector{Tuple{LongDNA{4}, LongDNA{4}, Int}}()
    for case in cases
        for guide in case.guides
            g = String(guide)
            for k in 1:3
                for _ in 1:reps_per_k
                    ref = mutate_for_micro(g, k, rng)
                    push!(inputs, (LongDNA{4}(g), LongDNA{4}(ref), k))
                end
            end
        end
    end
    return inputs
end

function run_micro_once(inputs, backend::Symbol)
    checksum = 0
    elapsed = @elapsed begin
        for (guide, ref, k) in inputs
            aln = if backend === :align
                CHOPOFF.align(guide, ref, k)
            else
                custom_traceback(guide, ref, k)
            end
            checksum += aln.dist
        end
    end
    return elapsed, checksum
end

function verify_micro_distance_parity(inputs)
    mismatch_count = 0
    for (guide, ref, k) in inputs
        align_dist = CHOPOFF.align(guide, ref, k).dist
        custom_dist = custom_traceback(guide, ref, k).dist
        mismatch_count += (align_dist == custom_dist) ? 0 : 1
    end
    return mismatch_count
end

function summarize_times(times::Vector{Float64})
    if isempty(times)
        return (median = NaN, mean = NaN, std = NaN, min = NaN, max = NaN)
    end
    return (
        median = median(times),
        mean = mean(times),
        std = std(times),
        min = minimum(times),
        max = maximum(times),
    )
end

function run_end_to_end_once(
    case,
    motif::Motif,
    distance::Int,
    strict_pam::Bool,
    backend::Symbol,
    output_path::String,
)
    es = fill(1_000_000, distance + 1)
    elapsed = @elapsed CHOPOFF.search_sassy(
        case.guides,
        case.genome,
        motif,
        output_path;
        distance = distance,
        strict_pam = strict_pam,
        early_stopping = es,
        traceback_backend = backend,
    )
    return elapsed
end

function run_end_to_end_case(
    case,
    distance::Int,
    strict_pam::Bool,
    runs::Int,
    tmp_root::String,
)
    case_tag = "$(case.label)_d$(distance)_strict$(strict_pam)"
    case_dir = joinpath(tmp_root, case_tag)
    mkpath(case_dir)

    motif = Motif(case.motif_name; distance = distance)
    align_warm = joinpath(case_dir, "warm_align.csv")
    custom_warm = joinpath(case_dir, "warm_custom.csv")
    run_end_to_end_once(case, motif, distance, strict_pam, :align, align_warm)
    run_end_to_end_once(case, motif, distance, strict_pam, :custom, custom_warm)

    align_times = Float64[]
    custom_times = Float64[]
    parity_flags = Bool[]

    for run_id in 1:runs
        align_out = joinpath(case_dir, "align_run_$(run_id).csv")
        custom_out = joinpath(case_dir, "custom_run_$(run_id).csv")
        align_elapsed = NaN
        custom_elapsed = NaN

        if isodd(run_id)
            align_elapsed = run_end_to_end_once(case, motif, distance, strict_pam, :align, align_out)
            custom_elapsed = run_end_to_end_once(case, motif, distance, strict_pam, :custom, custom_out)
        else
            custom_elapsed = run_end_to_end_once(case, motif, distance, strict_pam, :custom, custom_out)
            align_elapsed = run_end_to_end_once(case, motif, distance, strict_pam, :align, align_out)
        end

        parity_ok, lhs_only, rhs_only = check_core_parity(align_out, custom_out)
        push!(align_times, align_elapsed)
        push!(custom_times, custom_elapsed)
        push!(parity_flags, parity_ok)

        @printf(
            "  run %d/%d | align=%.6fs | custom=%.6fs | parity=%s\n",
            run_id,
            runs,
            align_elapsed,
            custom_elapsed,
            parity_ok ? "PASS" : "FAIL",
        )
        if !parity_ok
            println("    align-only rows: ", nrow(lhs_only))
            if nrow(lhs_only) > 0
                println(first(lhs_only, min(nrow(lhs_only), 5)))
            end
            println("    custom-only rows: ", nrow(rhs_only))
            if nrow(rhs_only) > 0
                println(first(rhs_only, min(nrow(rhs_only), 5)))
            end
        end
    end

    valid_idx = findall(parity_flags)
    align_valid = align_times[valid_idx]
    custom_valid = custom_times[valid_idx]
    align_stats = summarize_times(align_valid)
    custom_stats = summarize_times(custom_valid)
    speedup = (
        isfinite(align_stats.median) && isfinite(custom_stats.median) && custom_stats.median > 0
    ) ? (align_stats.median / custom_stats.median) : NaN

    @printf(
        "  summary %s | valid=%d/%d | median(s): align=%.6f custom=%.6f | speedup(align/custom)=%.3f\n",
        case_tag,
        length(valid_idx),
        runs,
        align_stats.median,
        custom_stats.median,
        speedup,
    )

    return (
        case_tag = case_tag,
        align_times = align_times,
        custom_times = custom_times,
        parity_flags = parity_flags,
        align_valid = align_valid,
        custom_valid = custom_valid,
        align_stats = align_stats,
        custom_stats = custom_stats,
        speedup = speedup,
    )
end

function main()
    runs = parse_int_env("SASSY_TRACEBACK_RUNS", 11)
    cas9_limit = parse_int_env("SASSY_TRACEBACK_CAS9_GUIDES", 6)
    micro_reps_per_k = parse_int_env("SASSY_TRACEBACK_MICRO_REPS_PER_K", 25)
    strict_modes = parse_strict_modes_env()
    distances = parse_distances_env()
    min_speedup = parse_float_env("SASSY_TRACEBACK_MIN_SPEEDUP", 0.10)
    min_e2e_speedup = parse_float_env("SASSY_TRACEBACK_MIN_E2E_SPEEDUP", 0.0)
    enforce = parse_bool_env("SASSY_TRACEBACK_ENFORCE", false)
    keep_tmp = parse_bool_env("SASSY_TRACEBACK_KEEP_TMP", false)
    out_csv = strip(get(ENV, "SASSY_TRACEBACK_OUT", ""))

    cases = load_cases(cas9_limit)
    tmp_root = mktempdir()

    println("Sassy traceback benchmark")
    println("Runs per section: $runs")
    println("Distances: $(join(distances, ','))")
    println("Strict modes: ", strict_modes == [true, false] ? "both" : (strict_modes[1] ? "strict" : "nonstrict"))
    println("Cases: ", join([c.label for c in cases], ", "))
    println("Enforce speed gate: ", enforce)
    println("Minimum required traceback speedup (micro): ", min_speedup * 100, "%")
    println("Minimum required end-to-end speedup: ", min_e2e_speedup * 100, "%")
    println("Temporary output: $tmp_root")

    summary_rows = NamedTuple[]

    println("\n== Microbenchmark (traceback only) ==")
    micro_inputs = build_micro_inputs(cases, micro_reps_per_k)
    println("Micro inputs: ", length(micro_inputs))
    mismatch_count = verify_micro_distance_parity(micro_inputs)
    println("Distance parity mismatches (align vs custom): ", mismatch_count)
    mismatch_count == 0 || error("Micro distance parity failed with $mismatch_count mismatches.")

    run_micro_once(micro_inputs, :align)
    run_micro_once(micro_inputs, :custom)

    micro_align_times = Float64[]
    micro_custom_times = Float64[]
    for run_id in 1:runs
        align_elapsed = NaN
        custom_elapsed = NaN
        if isodd(run_id)
            align_elapsed, _ = run_micro_once(micro_inputs, :align)
            custom_elapsed, _ = run_micro_once(micro_inputs, :custom)
        else
            custom_elapsed, _ = run_micro_once(micro_inputs, :custom)
            align_elapsed, _ = run_micro_once(micro_inputs, :align)
        end
        push!(micro_align_times, align_elapsed)
        push!(micro_custom_times, custom_elapsed)
        @printf("  run %d/%d | align=%.6fs | custom=%.6fs\n", run_id, runs, align_elapsed, custom_elapsed)
    end
    micro_align_stats = summarize_times(micro_align_times)
    micro_custom_stats = summarize_times(micro_custom_times)
    micro_speedup = (
        isfinite(micro_align_stats.median) && isfinite(micro_custom_stats.median) &&
        micro_custom_stats.median > 0
    ) ? (micro_align_stats.median / micro_custom_stats.median) : NaN
    @printf(
        "Micro median(s): align=%.6f custom=%.6f | speedup(align/custom)=%.3f\n",
        micro_align_stats.median,
        micro_custom_stats.median,
        micro_speedup,
    )
    push!(
        summary_rows,
        (
            timestamp = Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
            section = "micro",
            case = "all",
            distance = -1,
            strict_pam = missing,
            runs_total = runs,
            runs_valid = runs,
            align_median_s = micro_align_stats.median,
            custom_median_s = micro_custom_stats.median,
            speedup_align_over_custom = micro_speedup,
            parity_ok = true,
        ),
    )

    println("\n== End-to-end benchmark (search_sassy) ==")
    e2e_results = NamedTuple[]
    all_parity_ok = true
    pooled_align_valid = Float64[]
    pooled_custom_valid = Float64[]

    for case in cases
        for strict_pam in strict_modes
            for distance in distances
                println("\nCase=$(case.label) distance=$distance strict_pam=$strict_pam")
                result = run_end_to_end_case(case, distance, strict_pam, runs, tmp_root)
                push!(e2e_results, result)
                all_parity_ok &= all(result.parity_flags)
                append!(pooled_align_valid, result.align_valid)
                append!(pooled_custom_valid, result.custom_valid)

                push!(
                    summary_rows,
                    (
                        timestamp = Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
                        section = "end_to_end",
                        case = result.case_tag,
                        distance = distance,
                        strict_pam = strict_pam,
                        runs_total = runs,
                        runs_valid = length(result.align_valid),
                        align_median_s = result.align_stats.median,
                        custom_median_s = result.custom_stats.median,
                        speedup_align_over_custom = result.speedup,
                        parity_ok = all(result.parity_flags),
                    ),
                )
            end
        end
    end

    align_pooled_stats = summarize_times(pooled_align_valid)
    custom_pooled_stats = summarize_times(pooled_custom_valid)
    pooled_speedup = (
        isfinite(align_pooled_stats.median) && isfinite(custom_pooled_stats.median) &&
        custom_pooled_stats.median > 0
    ) ? (align_pooled_stats.median / custom_pooled_stats.median) : NaN
    @printf(
        "\nPooled end-to-end median(s): align=%.6f custom=%.6f | speedup(align/custom)=%.3f | parity_all=%s\n",
        align_pooled_stats.median,
        custom_pooled_stats.median,
        pooled_speedup,
        all_parity_ok ? "PASS" : "FAIL",
    )

    if !isempty(out_csv)
        out_dir = dirname(out_csv)
        mkpath(out_dir)
        CSV.write(out_csv, DataFrame(summary_rows))
        println("Wrote summary CSV: $out_csv")
    end

    required_micro_speedup = 1.0 + min_speedup
    required_e2e_speedup = 1.0 + min_e2e_speedup
    if enforce
        all_parity_ok || error("Speed gate failed: end-to-end parity mismatches were detected.")
        isfinite(micro_speedup) || error("Speed gate failed: micro speedup is not finite.")
        isfinite(pooled_speedup) || error("Speed gate failed: pooled speedup is not finite.")
        if micro_speedup < required_micro_speedup
            error(
                "Speed gate failed (micro): expected align/custom speedup >= " *
                "$(round(required_micro_speedup; digits = 3)), observed " *
                "$(round(micro_speedup; digits = 3)).",
            )
        end
        if pooled_speedup < required_e2e_speedup
            error(
                "Speed gate failed (end-to-end): expected align/custom speedup >= " *
                "$(round(required_e2e_speedup; digits = 3)), observed " *
                "$(round(pooled_speedup; digits = 3)).",
            )
        end
    end

    if !keep_tmp
        try
            rm(tmp_root; recursive = true, force = true)
        catch
        end
    end
end

main()
