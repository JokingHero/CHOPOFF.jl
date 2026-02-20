#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Dates
using Printf
using Statistics

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const SAMPLE_DIR = joinpath(ROOT_DIR, "test", "sample_data")
const CORE_PARITY_COLS = [:guide, :distance, :chromosome, :start, :strand]
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

function parse_bool_env(name::String, default::Bool)
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "y")
        return true
    elseif raw in ("0", "false", "no", "n")
        return false
    end
    error("Invalid boolean in environment variable $name: '$raw'")
end

function resolve_modes()
    mode = lowercase(strip(get(ENV, "CHOPOFF_BENCH_MODE", "both")))
    if mode == "both"
        return ["single", "env"]
    elseif mode in ("single", "env")
        return [mode]
    end
    error("Invalid CHOPOFF_BENCH_MODE='$mode'. Allowed values: both, single, env.")
end

function maybe_cpu_name()
    try
        return string(Sys.CPU_NAME)
    catch
        return "unknown"
    end
end

function require_file(path::String)
    isfile(path) || error("Required file not found: $path")
    return path
end

function load_benchmark_cases(distance::Int)
    cas9_genome = require_file(joinpath(SAMPLE_DIR, "genome", "semirandom.fa"))
    require_file(cas9_genome * ".fai")
    cas9_guides_path = require_file(joinpath(SAMPLE_DIR, "guides.txt"))
    cas9_guides = LongDNA{4}.(sort(collect(Set(readlines(cas9_guides_path)))))

    cas12a_genome = require_file(joinpath(SAMPLE_DIR, "genome", "semirandom.2bit"))
    cas12a_guides = copy(CAS12A_GUIDES)

    return [
        (
            label = "Cas9",
            motif = Motif("Cas9"; distance = distance),
            genome = cas9_genome,
            guides = cas9_guides,
        ),
        (
            label = "Cas12a",
            motif = Motif("Cas12a"; distance = distance),
            genome = cas12a_genome,
            guides = cas12a_guides,
        ),
    ]
end

function read_result(path::String)
    isfile(path) || error("Missing result file: $path")
    first_line = open(path, "r") do io
        eof(io) ? "" : replace(chomp(readline(io)), "\ufeff" => "")
    end

    if isempty(first_line)
        return DataFrame(
            guide = String[],
            alignment_guide = String[],
            alignment_reference = String[],
            distance = Int[],
            chromosome = String[],
            start = Int[],
            strand = String[],
        )
    end

    has_header = startswith(lowercase(first_line), "guide,")
    df = has_header ? DataFrame(CSV.File(path)) : DataFrame(CSV.File(path; header = false))

    if !has_header
        expected = [
            :guide,
            :alignment_guide,
            :alignment_reference,
            :distance,
            :chromosome,
            :start,
            :strand,
        ]
        if ncol(df) == 0
            return DataFrame(
                guide = String[],
                alignment_guide = String[],
                alignment_reference = String[],
                distance = Int[],
                chromosome = String[],
                start = Int[],
                strand = String[],
            )
        elseif ncol(df) < length(expected)
            error("Unexpected number of columns ($(ncol(df))) in headerless results: $path")
        end
        rename!(df, Dict(names(df)[i] => expected[i] for i in 1:length(expected)))
    end

    return df
end

function normalize_core_parity(df::DataFrame)
    cols = Set(Symbol.(names(df)))
    for c in CORE_PARITY_COLS
        c in cols || error("Missing required parity column '$c' in results.")
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
    lhs = normalize_core_parity(read_result(lhs_path))
    rhs = normalize_core_parity(read_result(rhs_path))
    if lhs == rhs
        return true, DataFrame(), DataFrame()
    end
    lhs_only = antijoin(lhs, rhs, on = CORE_PARITY_COLS)
    rhs_only = antijoin(rhs, lhs, on = CORE_PARITY_COLS)
    return false, lhs_only, rhs_only
end

function summarize_times(times::Vector{Float64})
    if isempty(times)
        return (
            median = NaN,
            mean = NaN,
            std = NaN,
            min = NaN,
            max = NaN,
        )
    end
    return (
        median = median(times),
        mean = mean(times),
        std = std(times),
        min = minimum(times),
        max = maximum(times),
    )
end

function run_case(case, distance::Int, runs::Int, mode::AbstractString, tmp_root::String)
    println("\n=== Case: $(case.label) | distance=$distance | guides=$(length(case.guides)) ===")

    case_dir = joinpath(tmp_root, "$(case.label)_$mode")
    db_dir = joinpath(case_dir, "prefixHashDB")
    mkpath(db_dir)

    build_prefixHashDB("bench_$(case.label)_d$distance", case.genome, case.motif, db_dir)
    es = fill(1_000_000, distance + 1)

    warm_sassy = joinpath(case_dir, "warmup_sassy.csv")
    warm_prefix = joinpath(case_dir, "warmup_prefixhash.csv")
    search_sassy(case.guides, case.genome, case.motif, warm_sassy; distance = distance, early_stopping = es)
    search_prefixHashDB(db_dir, case.guides, warm_prefix; distance = distance, early_stopping = es)
    warm_ok, warm_lhs_only, warm_rhs_only = check_core_parity(warm_sassy, warm_prefix)
    println("Warmup parity: ", warm_ok ? "PASS" : "FAIL")
    if !warm_ok
        println("Warmup Sassy-only tuples: ", nrow(warm_lhs_only))
        println("Warmup PrefixHashDB-only tuples: ", nrow(warm_rhs_only))
    end

    sassy_times = Float64[]
    prefix_times = Float64[]
    parity_flags = Bool[]

    for run_id in 1:runs
        sassy_out = joinpath(case_dir, "sassy_run_$(run_id).csv")
        prefix_out = joinpath(case_dir, "prefixhash_run_$(run_id).csv")

        sassy_elapsed = @elapsed search_sassy(
            case.guides,
            case.genome,
            case.motif,
            sassy_out;
            distance = distance,
            early_stopping = es,
        )
        prefix_elapsed = @elapsed search_prefixHashDB(
            db_dir,
            case.guides,
            prefix_out;
            distance = distance,
            early_stopping = es,
        )

        parity_ok, lhs_only, rhs_only = check_core_parity(sassy_out, prefix_out)
        push!(sassy_times, sassy_elapsed)
        push!(prefix_times, prefix_elapsed)
        push!(parity_flags, parity_ok)

        @printf(
            "run %d/%d | sassy=%.6fs | prefixHashDB=%.6fs | parity=%s\n",
            run_id,
            runs,
            sassy_elapsed,
            prefix_elapsed,
            parity_ok ? "PASS" : "FAIL",
        )
        if !parity_ok
            println("  Sassy-only tuples: ", nrow(lhs_only))
            if nrow(lhs_only) > 0
                println(first(lhs_only, min(nrow(lhs_only), 5)))
            end
            println("  PrefixHashDB-only tuples: ", nrow(rhs_only))
            if nrow(rhs_only) > 0
                println(first(rhs_only, min(nrow(rhs_only), 5)))
            end
        end
    end

    valid_idx = findall(parity_flags)
    sassy_valid = sassy_times[valid_idx]
    prefix_valid = prefix_times[valid_idx]
    sassy_stats = summarize_times(sassy_valid)
    prefix_stats = summarize_times(prefix_valid)
    parity_ok = all(parity_flags)

    ratio_prefix_over_sassy = (
        isfinite(prefix_stats.median) && isfinite(sassy_stats.median) && sassy_stats.median > 0
    ) ? (prefix_stats.median / sassy_stats.median) : NaN
    ratio_sassy_over_prefix = (
        isfinite(prefix_stats.median) && isfinite(sassy_stats.median) && prefix_stats.median > 0
    ) ? (sassy_stats.median / prefix_stats.median) : NaN

    @printf(
        "Summary %s | valid_runs=%d/%d | median(s): sassy=%.6f prefixHashDB=%.6f | ratio(prefix/sassy)=%.3f | parity_all=%s\n",
        case.label,
        length(valid_idx),
        runs,
        sassy_stats.median,
        prefix_stats.median,
        ratio_prefix_over_sassy,
        parity_ok ? "PASS" : "FAIL",
    )

    stamp = Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS")
    rows = [
        (
            timestamp = stamp,
            motif = case.label,
            distance = distance,
            thread_mode = mode,
            threads = Threads.nthreads(),
            algorithm = "sassy",
            runs_total = runs,
            runs_valid = length(valid_idx),
            median_s = sassy_stats.median,
            mean_s = sassy_stats.mean,
            std_s = sassy_stats.std,
            min_s = sassy_stats.min,
            max_s = sassy_stats.max,
            parity_ok = parity_ok,
            ratio_vs_other = ratio_sassy_over_prefix,
        ),
        (
            timestamp = stamp,
            motif = case.label,
            distance = distance,
            thread_mode = mode,
            threads = Threads.nthreads(),
            algorithm = "prefixHashDB",
            runs_total = runs,
            runs_valid = length(valid_idx),
            median_s = prefix_stats.median,
            mean_s = prefix_stats.mean,
            std_s = prefix_stats.std,
            min_s = prefix_stats.min,
            max_s = prefix_stats.max,
            parity_ok = parity_ok,
            ratio_vs_other = ratio_prefix_over_sassy,
        ),
    ]
    return rows
end

function run_child()
    runs = parse_int_env("CHOPOFF_BENCH_RUNS", 7)
    distance = parse_int_env("CHOPOFF_BENCH_DISTANCE", 3)
    keep_tmp = parse_bool_env("CHOPOFF_BENCH_KEEP_TMP", false)
    mode = strip(get(ENV, "CHOPOFF_BENCH_CHILD_MODE", "env"))
    out_path = strip(get(ENV, "CHOPOFF_BENCH_OUT", ""))
    append_out = parse_bool_env("CHOPOFF_BENCH_APPEND", false)

    println("CHOPOFF benchmark child mode")
    println("time_utc: ", Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println("mode: $mode")
    println("threads: ", Threads.nthreads())
    println("cpu: ", maybe_cpu_name())
    println("runs: $runs")
    println("distance: $distance")

    cases = load_benchmark_cases(distance)
    tmp_root = mktempdir(prefix = "chopoff_speed_")
    rows = NamedTuple[]
    for case in cases
        append!(rows, run_case(case, distance, runs, mode, tmp_root))
    end

    summary = DataFrame(rows)
    if out_path != ""
        mkpath(dirname(out_path))
        if append_out && isfile(out_path)
            CSV.write(out_path, summary; append = true, writeheader = false)
        else
            CSV.write(out_path, summary)
        end
        println("\nWrote summary CSV: $out_path")
    end

    if keep_tmp
        println("Keeping temporary benchmark directory: $tmp_root")
    else
        rm(tmp_root; recursive = true, force = true)
    end
end

function run_parent()
    modes = resolve_modes()
    user_out = strip(get(ENV, "CHOPOFF_BENCH_OUT", ""))
    summary_path = user_out == "" ? joinpath(mktempdir(prefix = "chopoff_bench_summary_"), "summary.csv") : abspath(user_out)
    if isfile(summary_path)
        rm(summary_path; force = true)
    end

    println("CHOPOFF speed benchmark parent mode")
    println("project_root: $ROOT_DIR")
    println("modes: ", join(modes, ", "))
    println("summary_path: $summary_path")

    child_cmd = `$(Base.julia_cmd()) --project=$ROOT_DIR $(@__FILE__)`
    for (idx, mode) in enumerate(modes)
        env_copy = copy(ENV)
        env_copy["CHOPOFF_BENCH_CHILD"] = "1"
        env_copy["CHOPOFF_BENCH_CHILD_MODE"] = mode
        env_copy["CHOPOFF_BENCH_OUT"] = summary_path
        env_copy["CHOPOFF_BENCH_APPEND"] = idx == 1 ? "0" : "1"
        if mode == "single"
            env_copy["JULIA_NUM_THREADS"] = "1"
        end

        println("\nLaunching mode '$mode' ...")
        run(setenv(child_cmd, env_copy))
    end

    println("\nCombined summary")
    summary = DataFrame(CSV.File(summary_path))
    show(summary; allrows = true, allcols = true, truncate = 0)
    println()
    println("Summary CSV: $summary_path")
end

if get(ENV, "CHOPOFF_BENCH_CHILD", "0") == "1"
    run_child()
else
    run_parent()
end
