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

    selected = Set(split(lowercase(strip(get(ENV, "CHOPOFF_BENCH_CASES", "cas9,cas12a"))), ','))
    cases = NamedTuple[]
    if "cas9" in selected
        push!(cases, (
            label = "Cas9",
            motif = Motif("Cas9"; distance = distance),
            genome = cas9_genome,
            guides = cas9_guides,
        ))
    end
    if "cas12a" in selected
        push!(cases, (
            label = "Cas12a",
            motif = Motif("Cas12a"; distance = distance),
            genome = cas12a_genome,
            guides = cas12a_guides,
        ))
    end
    isempty(cases) && error("No benchmark cases selected by CHOPOFF_BENCH_CASES.")
    return cases
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


function scan_stat_fields(stats)
    if stats === nothing
        return (
            scan_motif_candidates = missing,
            scan_ambiguous_prefixes = missing,
            scan_prefix_hits = missing,
            scan_guide_pairs = missing,
            scan_alignment_calls = missing,
            scan_emitted_rows = missing,
            scan_path_rows = missing,
            scan_query_hashes = missing,
            scan_bruteforce_guide_pairs = missing,
            scan_path_source = missing,
            scan_query_variant = missing,
            scan_query_build_s = missing,
            scan_path_load_s = missing,
            scan_query_hash_s = missing,
            scan_query_format_s = missing,
            scan_query_fold_s = missing,
            scan_query_dedup_s = missing,
            scan_query_insert_s = missing,
            scan_query_lookup_s = missing,
            scan_chrom_load_s = missing,
            scan_findguides_s = missing,
            scan_candidate_prefix_s = missing,
            scan_candidate_hash_s = missing,
            scan_candidate_materialize_s = missing,
            scan_align_s = missing,
            scan_emit_s = missing,
            scan_scan_s = missing,
            scan_verify_s = missing,
        )
    end
    return (
        scan_motif_candidates = stats.motif_candidates,
        scan_ambiguous_prefixes = stats.ambiguous_prefixes,
        scan_prefix_hits = stats.prefix_hits,
        scan_guide_pairs = stats.guide_pairs,
        scan_alignment_calls = stats.alignment_calls,
        scan_emitted_rows = stats.emitted_rows,
        scan_path_rows = stats.path_rows,
        scan_query_hashes = stats.query_hashes,
        scan_bruteforce_guide_pairs = stats.bruteforce_guide_pairs,
        scan_path_source = string(stats.path_source),
        scan_query_variant = string(stats.query_variant),
        scan_query_build_s = Float64(stats.query_build_ns) / 1e9,
        scan_path_load_s = Float64(stats.path_load_ns) / 1e9,
        scan_query_hash_s = Float64(stats.query_hash_ns) / 1e9,
        scan_query_format_s = Float64(stats.query_format_ns) / 1e9,
        scan_query_fold_s = Float64(stats.query_fold_ns) / 1e9,
        scan_query_dedup_s = Float64(stats.query_dedup_ns) / 1e9,
        scan_query_insert_s = Float64(stats.query_insert_ns) / 1e9,
        scan_query_lookup_s = Float64(stats.query_lookup_ns) / 1e9,
        scan_chrom_load_s = Float64(stats.chrom_load_ns) / 1e9,
        scan_findguides_s = Float64(stats.findguides_ns) / 1e9,
        scan_candidate_prefix_s = Float64(stats.candidate_prefix_ns) / 1e9,
        scan_candidate_hash_s = Float64(stats.candidate_hash_ns) / 1e9,
        scan_candidate_materialize_s = Float64(stats.candidate_materialize_ns) / 1e9,
        scan_align_s = Float64(stats.align_ns) / 1e9,
        scan_emit_s = Float64(stats.emit_ns) / 1e9,
        scan_scan_s = Float64(stats.scan_ns) / 1e9,
        scan_verify_s = Float64(stats.verify_ns) / 1e9,
    )
end

function run_case(case, distance::Int, runs::Int, mode::AbstractString, tmp_root::String, scan_enabled::Bool)
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
    if scan_enabled
        warm_scan = joinpath(case_dir, "warmup_prefixhashscan.csv")
        CHOPOFF.search_prefixHashScan(case.guides, case.genome, case.motif, warm_scan; distance = distance, early_stopping = es)
        scan_warm_ok, scan_warm_lhs, scan_warm_rhs = check_core_parity(warm_scan, warm_prefix)
        println("Warmup prefixHashScan parity: ", scan_warm_ok ? "PASS" : "FAIL")
        if !scan_warm_ok
            println("Warmup prefixHashScan-only tuples: ", nrow(scan_warm_lhs))
            println("Warmup PrefixHashDB-only tuples: ", nrow(scan_warm_rhs))
        end
    end
    warm_ok, warm_lhs_only, warm_rhs_only = check_core_parity(warm_sassy, warm_prefix)
    println("Warmup parity: ", warm_ok ? "PASS" : "FAIL")
    if !warm_ok
        println("Warmup Sassy-only tuples: ", nrow(warm_lhs_only))
        println("Warmup PrefixHashDB-only tuples: ", nrow(warm_rhs_only))
    end

    sassy_times = Float64[]
    prefix_times = Float64[]
    scan_times = Float64[]
    scan_stats_runs = CHOPOFF.PrefixHashScanStats[]
    parity_flags = Bool[]
    scan_parity_flags = Bool[]

    for run_id in 1:runs
        sassy_out = joinpath(case_dir, "sassy_run_$(run_id).csv")
        prefix_out = joinpath(case_dir, "prefixhash_run_$(run_id).csv")
        scan_out = joinpath(case_dir, "prefixhashscan_run_$(run_id).csv")

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
        scan_elapsed = NaN
        scan_parity_ok = true
        if scan_enabled
            scan_stats = CHOPOFF.PrefixHashScanStats()
            scan_elapsed = @elapsed CHOPOFF.search_prefixHashScan(
                case.guides,
                case.genome,
                case.motif,
                scan_out;
                distance = distance,
                early_stopping = es,
                stats = scan_stats,
            )
            scan_parity_ok, scan_lhs_only, scan_rhs_only = check_core_parity(scan_out, prefix_out)
            push!(scan_times, scan_elapsed)
            push!(scan_stats_runs, scan_stats)
            push!(scan_parity_flags, scan_parity_ok)
        end

        parity_ok, lhs_only, rhs_only = check_core_parity(sassy_out, prefix_out)
        push!(sassy_times, sassy_elapsed)
        push!(prefix_times, prefix_elapsed)
        push!(parity_flags, parity_ok)

        if scan_enabled
            @printf(
                "run %d/%d | sassy=%.6fs | prefixHashDB=%.6fs | prefixHashScan=%.6fs | parity=%s | scan_parity=%s\n",
                run_id,
                runs,
                sassy_elapsed,
                prefix_elapsed,
                scan_elapsed,
                parity_ok ? "PASS" : "FAIL",
                scan_parity_ok ? "PASS" : "FAIL",
            )
        else
            @printf(
                "run %d/%d | sassy=%.6fs | prefixHashDB=%.6fs | parity=%s\n",
                run_id,
                runs,
                sassy_elapsed,
                prefix_elapsed,
                parity_ok ? "PASS" : "FAIL",
            )
        end
        if scan_enabled
            @printf(
                "  scan stats | source=%s | variant=%s | path_rows=%d | query_hashes=%d | brute_pairs=%d | candidates=%d | prefix_hits=%d | guide_pairs=%d | align_calls=%d | rows=%d | ambiguous_prefixes=%d | query=%.3fs | load=%.3fs | hash=%.3fs | format=%.3fs | fold=%.3fs | dedup=%.3fs | insert=%.3fs | lookup=%.3fs | chrom=%.3fs | find=%.3fs | prefix=%.3fs | cand_hash=%.3fs | materialize=%.3fs | align=%.3fs | emit=%.3fs | scan=%.3fs | verify=%.3fs\n",
                string(scan_stats.path_source),
                string(scan_stats.query_variant),
                scan_stats.path_rows,
                scan_stats.query_hashes,
                scan_stats.bruteforce_guide_pairs,
                scan_stats.motif_candidates,
                scan_stats.prefix_hits,
                scan_stats.guide_pairs,
                scan_stats.alignment_calls,
                scan_stats.emitted_rows,
                scan_stats.ambiguous_prefixes,
                Float64(scan_stats.query_build_ns) / 1e9,
                Float64(scan_stats.path_load_ns) / 1e9,
                Float64(scan_stats.query_hash_ns) / 1e9,
                Float64(scan_stats.query_format_ns) / 1e9,
                Float64(scan_stats.query_fold_ns) / 1e9,
                Float64(scan_stats.query_dedup_ns) / 1e9,
                Float64(scan_stats.query_insert_ns) / 1e9,
                Float64(scan_stats.query_lookup_ns) / 1e9,
                Float64(scan_stats.chrom_load_ns) / 1e9,
                Float64(scan_stats.findguides_ns) / 1e9,
                Float64(scan_stats.candidate_prefix_ns) / 1e9,
                Float64(scan_stats.candidate_hash_ns) / 1e9,
                Float64(scan_stats.candidate_materialize_ns) / 1e9,
                Float64(scan_stats.align_ns) / 1e9,
                Float64(scan_stats.emit_ns) / 1e9,
                Float64(scan_stats.scan_ns) / 1e9,
                Float64(scan_stats.verify_ns) / 1e9,
            )
        end

        if scan_enabled && !scan_parity_ok
            println("  prefixHashScan-only tuples: ", nrow(scan_lhs_only))
            if nrow(scan_lhs_only) > 0
                println(first(scan_lhs_only, min(nrow(scan_lhs_only), 5)))
            end
            println("  PrefixHashDB-only tuples vs scan: ", nrow(scan_rhs_only))
            if nrow(scan_rhs_only) > 0
                println(first(scan_rhs_only, min(nrow(scan_rhs_only), 5)))
            end
        end

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
    scan_valid_idx = findall(scan_parity_flags)
    scan_stats = summarize_times(scan_times[scan_valid_idx])
    scan_summary_stats = isempty(scan_valid_idx) ? nothing : scan_stats_runs[last(scan_valid_idx)]
    parity_ok = all(parity_flags)
    scan_parity_ok = !scan_enabled || all(scan_parity_flags)

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
    rows = Any[
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
            scan_stat_fields(nothing)...,
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
            scan_stat_fields(nothing)...,
        ),
    ]

    if scan_enabled
        push!(rows, (
            timestamp = stamp,
            motif = case.label,
            distance = distance,
            thread_mode = mode,
            threads = Threads.nthreads(),
            algorithm = "prefixHashScan",
            runs_total = runs,
            runs_valid = length(scan_valid_idx),
            median_s = scan_stats.median,
            mean_s = scan_stats.mean,
            std_s = scan_stats.std,
            min_s = scan_stats.min,
            max_s = scan_stats.max,
            parity_ok = scan_parity_ok,
            ratio_vs_other = isfinite(scan_stats.median) && isfinite(prefix_stats.median) && prefix_stats.median > 0 ? scan_stats.median / prefix_stats.median : NaN,
            scan_stat_fields(scan_summary_stats)...,
        ))
    end
    return rows
end

function run_child()
    runs = parse_int_env("CHOPOFF_BENCH_RUNS", 7)
    distance = parse_int_env("CHOPOFF_BENCH_DISTANCE", 3)
    keep_tmp = parse_bool_env("CHOPOFF_BENCH_KEEP_TMP", false)
    scan_enabled = parse_bool_env("CHOPOFF_BENCH_SCAN", false)
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
    println("scan_enabled: $scan_enabled")
    println("cases: ", get(ENV, "CHOPOFF_BENCH_CASES", "cas9,cas12a"))

    cases = load_benchmark_cases(distance)
    tmp_root = mktempdir(prefix = "chopoff_speed_")
    rows = NamedTuple[]
    for case in cases
        append!(rows, run_case(case, distance, runs, mode, tmp_root, scan_enabled))
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
