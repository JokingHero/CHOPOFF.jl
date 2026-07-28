#!/usr/bin/env julia

include(joinpath(@__DIR__, "run_human_sassy.jl"))
using Statistics

const DEFAULT_CAS12A_GUIDES =
    joinpath(@__DIR__, "data", "cas12a_guides_for_tests.txt")

function parse_csv_values(name::String, default::String)
    raw = strip(get(ENV, name, default))
    values = String[
        strip(value) for value in split(raw, ',') if !isempty(strip(value))]
    isempty(values) && error("$name must contain at least one value.")
    return values
end

function sweep_motifs()
    motifs = parse_csv_values(
        "CHOPOFF_HUMAN_SWEEP_MOTIFS", "Cas9,Cas12a")
    all(in(("Cas9", "Cas12a")), motifs) ||
        error("CHOPOFF_HUMAN_SWEEP_MOTIFS supports only Cas9 and Cas12a.")
    return unique(motifs)
end

function sweep_distances()
    distances = parse.(Int, parse_csv_values(
        "CHOPOFF_HUMAN_SWEEP_DISTANCES", "0,1,2,3,4"))
    all(in(0:4), distances) ||
        error("CHOPOFF_HUMAN_SWEEP_DISTANCES must be within 0:4.")
    return unique(distances)
end

function sweep_output_dir()
    configured = strip(get(ENV, "CHOPOFF_HUMAN_SWEEP_OUT", ""))
    if !isempty(configured)
        return abspath(configured)
    end
    stamp = Dates.format(now(Dates.UTC), dateformat"yyyymmdd_HHMMSS")
    return joinpath(@__DIR__, "outputs", "prefix_sweep_" * stamp)
end

function sweep_index_parent()
    configured = strip(get(
        ENV, "CHOPOFF_HUMAN_SWEEP_INDEX_PARENT",
        joinpath(@__DIR__, "indexes")))
    return abspath(configured)
end

function genome_key(genome::String)
    return replace(basename(genome), r"[^A-Za-z0-9_.-]" => "_")
end

function sweep_index_dir(
    index_parent::String, genome::String, motif_name::String)

    return joinpath(
        index_parent, "$(genome_key(genome))_$(motif_name)_d4")
end

function directory_size(path::String)
    total = 0
    for (root, _, files) in walkdir(path), file in files
        total += stat(joinpath(root, file)).size
    end
    return total
end

function guide_path(motif_name::String)
    env_name = motif_name == "Cas9" ?
        "CHOPOFF_HUMAN_SWEEP_CAS9_GUIDES" :
        "CHOPOFF_HUMAN_SWEEP_CAS12A_GUIDES"
    default = motif_name == "Cas9" ? DEFAULT_GUIDES : DEFAULT_CAS12A_GUIDES
    return abspath(strip(get(ENV, env_name, default)))
end

function load_sweep_guides(motif_name::String, distance::Int)
    guides = load_guides(
        guide_path(motif_name), Motif(motif_name; distance = distance))
    limit = parse_int_env("CHOPOFF_HUMAN_SWEEP_GUIDE_LIMIT", 0)
    limit >= 0 || error("CHOPOFF_HUMAN_SWEEP_GUIDE_LIMIT must be nonnegative.")
    limit == 0 && return guides
    return guides[1:min(limit, length(guides))]
end

function print_progress(message)
    println(message)
    flush(stdout)
end

function write_search_tables(
    output_dir::String,
    timing_rows,
    summary_rows,
    parity_rows,
    stats_rows)

    CSV.write(joinpath(output_dir, "timings.csv"), DataFrame(timing_rows))
    CSV.write(joinpath(output_dir, "summary.csv"), DataFrame(summary_rows))
    CSV.write(joinpath(output_dir, "parity.csv"), DataFrame(parity_rows))
    CSV.write(
        joinpath(output_dir, "prefixhashscan_stats.csv"),
        DataFrame(stats_rows))
    return nothing
end

function run_build_stage()
    Threads.nthreads() == 1 ||
        error("The prefixHashDB build stage requires JULIA_NUM_THREADS=1.")
    output_dir = sweep_output_dir()
    mkpath(output_dir)
    genome = validate_genome(abspath(strip(get(
        ENV, "CHOPOFF_HUMAN_GENOME", DEFAULT_GENOME))))
    index_parent = sweep_index_parent()
    mkpath(index_parent)
    rebuild = parse_bool_env("CHOPOFF_HUMAN_SWEEP_REBUILD", false)
    rows = NamedTuple[]

    for motif_name in sweep_motifs()
        index_dir = sweep_index_dir(index_parent, genome, motif_name)
        db_file = joinpath(index_dir, "prefixHashDB.bin")
        reused = !rebuild && isfile(db_file)
        elapsed = 0.0
        if reused
            print_progress("Reusing d4 prefixHashDB: $index_dir")
        else
            print_progress("Building $motif_name d4 prefixHashDB: $index_dir")
            mkpath(index_dir)
            motif = Motif(motif_name; distance = 4)
            elapsed = @elapsed build_prefixHashDB(
                "human_$(motif_name)_d4", genome, motif, index_dir)
            print_progress(
                "Finished $motif_name d4 prefixHashDB in $(round(elapsed; digits=3)) s")
        end
        push!(rows, (
            timestamp = Dates.format(
                now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
            genome = genome,
            motif = motif_name,
            build_distance = 4,
            threads = Threads.nthreads(),
            reused = reused,
            elapsed_s = elapsed,
            index_bytes = directory_size(index_dir),
            index_dir = index_dir,
        ))
        CSV.write(joinpath(output_dir, "builds.csv"), DataFrame(rows))
    end
    return nothing
end

function run_prefix_db(
    index_dir::String,
    guides::Vector{LongDNA{4}},
    output::String,
    distance::Int,
    early_stopping::Vector{Int})

    return search_prefixHashDB(
        index_dir, guides, output;
        distance = distance, early_stopping = early_stopping)
end

function run_prefix_scan(
    guides::Vector{LongDNA{4}},
    genome::String,
    motif::Motif,
    output::String,
    distance::Int,
    early_stopping::Vector{Int};
    stats::Union{Nothing, CHOPOFF.PrefixHashScanStats} = nothing)

    return CHOPOFF.search_prefixHashScan(
        guides, genome, motif, output;
        distance = distance,
        early_stopping = early_stopping,
        scan_threads = Threads.nthreads(),
        stats = stats)
end

function run_search_stage()
    expected_threads = parse_int_env("CHOPOFF_HUMAN_SWEEP_THREADS", 24)
    Threads.nthreads() == expected_threads ||
        error("Search stage expected $expected_threads Julia threads, got $(Threads.nthreads()).")
    runs = parse_int_env("CHOPOFF_HUMAN_SWEEP_RUNS", 5)
    runs >= 1 || error("CHOPOFF_HUMAN_SWEEP_RUNS must be positive.")
    output_dir = sweep_output_dir()
    mkpath(output_dir)
    cases_dir = joinpath(output_dir, "cases")
    mkpath(cases_dir)
    genome = validate_genome(abspath(strip(get(
        ENV, "CHOPOFF_HUMAN_GENOME", DEFAULT_GENOME))))
    index_parent = sweep_index_parent()

    timing_rows = NamedTuple[]
    summary_rows = NamedTuple[]
    parity_rows = NamedTuple[]
    stats_rows = NamedTuple[]

    for motif_name in sweep_motifs()
        index_dir = sweep_index_dir(index_parent, genome, motif_name)
        isfile(joinpath(index_dir, "prefixHashDB.bin")) ||
            error("Missing d4 prefixHashDB: $index_dir")
        guides = load_sweep_guides(motif_name, 4)
        for distance in sweep_distances()
            motif = Motif(motif_name; distance = distance)
            early_stopping = fill(typemax(Int), distance + 1)
            case_dir = joinpath(cases_dir, "$(motif_name)_d$(distance)")
            mkpath(case_dir)
            prefix_output = joinpath(case_dir, "prefixhash.csv")
            scan_output = joinpath(case_dir, "prefixhashscan.csv")

            print_progress(
                "Warmup motif=$motif_name distance=$distance guides=$(length(guides))")
            run_prefix_db(
                index_dir, guides, prefix_output, distance, early_stopping)
            run_prefix_scan(
                guides, genome, motif, scan_output, distance, early_stopping)

            elapsed = Dict(
                "prefixHashDB" => Float64[],
                "prefixHashScan" => Float64[],
            )
            for run_idx in 1:runs
                order = isodd(run_idx) ?
                    ("prefixHashDB", "prefixHashScan") :
                    ("prefixHashScan", "prefixHashDB")
                for algorithm in order
                    GC.gc()
                    output = algorithm == "prefixHashDB" ?
                        prefix_output : scan_output
                    seconds = if algorithm == "prefixHashDB"
                        @elapsed run_prefix_db(
                            index_dir, guides, output,
                            distance, early_stopping)
                    else
                        @elapsed run_prefix_scan(
                            guides, genome, motif, output,
                            distance, early_stopping)
                    end
                    push!(elapsed[algorithm], seconds)
                    rows = count_result_rows(output)
                    bytes = filesize_bytes(output)
                    push!(timing_rows, (
                        timestamp = Dates.format(
                            now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
                        genome = genome,
                        motif = motif_name,
                        distance = distance,
                        algorithm = algorithm,
                        run = run_idx,
                        execution_order = findfirst(==(algorithm), order),
                        threads = Threads.nthreads(),
                        guide_count = length(guides),
                        elapsed_s = seconds,
                        rows = rows,
                        bytes = bytes,
                        output = output,
                        index_dir = index_dir,
                    ))
                    print_progress(
                        "timed motif=$motif_name distance=$distance " *
                        "algorithm=$algorithm run=$run_idx/$runs " *
                        "elapsed=$(round(seconds; digits=3))s rows=$rows")
                end
            end

            scan_only, prefix_only = write_exact_parity_diffs(
                scan_output, prefix_output, case_dir)
            parity_ok = scan_only == 0 && prefix_only == 0
            push!(parity_rows, (
                motif = motif_name,
                distance = distance,
                guide_count = length(guides),
                prefix_rows = count_result_rows(prefix_output),
                scan_rows = count_result_rows(scan_output),
                scan_only = scan_only,
                prefix_only = prefix_only,
                parity = parity_ok,
                case_dir = case_dir,
            ))

            for algorithm in ("prefixHashDB", "prefixHashScan")
                values = elapsed[algorithm]
                push!(summary_rows, (
                    motif = motif_name,
                    distance = distance,
                    algorithm = algorithm,
                    threads = Threads.nthreads(),
                    guide_count = length(guides),
                    runs = length(values),
                    median_s = median(values),
                    mean_s = mean(values),
                    min_s = minimum(values),
                    max_s = maximum(values),
                    parity = parity_ok,
                    case_dir = case_dir,
                    index_dir = index_dir,
                ))
            end

            scan_stats = CHOPOFF.PrefixHashScanStats()
            diagnostic_elapsed = @elapsed run_prefix_scan(
                guides, genome, motif, scan_output, distance, early_stopping;
                stats = scan_stats)
            push!(stats_rows, (
                motif = motif_name,
                distance = distance,
                threads = Threads.nthreads(),
                guide_count = length(guides),
                diagnostic_elapsed_s = diagnostic_elapsed,
                path_source = string(scan_stats.path_source),
                scan_backend = string(scan_stats.scan_backend),
                path_rows = scan_stats.path_rows,
                query_hashes = scan_stats.query_hashes,
                motif_candidates = scan_stats.motif_candidates,
                prefix_hits = scan_stats.prefix_hits,
                guide_pairs = scan_stats.guide_pairs,
                alignment_calls = scan_stats.alignment_calls,
                distance_calls = scan_stats.distance_calls,
                traceback_calls = scan_stats.traceback_calls,
                emitted_rows = scan_stats.emitted_rows,
                metadata_s = Float64(scan_stats.metadata_ns) / 1e9,
                query_build_s = Float64(scan_stats.query_build_ns) / 1e9,
                path_load_s = Float64(scan_stats.path_load_ns) / 1e9,
                scan_s = Float64(scan_stats.scan_ns) / 1e9,
                verify_s = Float64(scan_stats.verify_ns) / 1e9,
                align_s = Float64(scan_stats.align_ns) / 1e9,
                output = scan_output,
            ))
            print_progress(
                "parity motif=$motif_name distance=$distance " *
                "status=$(parity_ok ? "PASS" : "FAIL") " *
                "scan_only=$scan_only prefix_only=$prefix_only")
            write_search_tables(
                output_dir, timing_rows, summary_rows, parity_rows, stats_rows)
        end
    end
    print_progress("Sweep complete. Summary: $(joinpath(output_dir, "summary.csv"))")
    return nothing
end

function required_sweep_output_dir()
    configured = strip(get(ENV, "CHOPOFF_HUMAN_SWEEP_OUT", ""))
    isempty(configured) &&
        error("CHOPOFF_HUMAN_SWEEP_OUT is required for parity stage.")
    output_dir = abspath(configured)
    isdir(output_dir) || error("Sweep output directory not found: $output_dir")
    return output_dir
end

function sweep_parity_row(
    output_dir::String,
    motif_name::String,
    distance::Int,
    guide_count::Int)

    case_dir = joinpath(output_dir, "cases", "$(motif_name)_d$(distance)")
    prefix_output = require_file(joinpath(case_dir, "prefixhash.csv"))
    scan_output = require_file(joinpath(case_dir, "prefixhashscan.csv"))
    scan_only, prefix_only = write_exact_parity_diffs(
        scan_output, prefix_output, case_dir)
    return (
        motif = motif_name,
        distance = distance,
        guide_count = guide_count,
        prefix_rows = count_result_rows(prefix_output),
        scan_rows = count_result_rows(scan_output),
        scan_only = scan_only,
        prefix_only = prefix_only,
        parity = scan_only == 0 && prefix_only == 0,
        case_dir = case_dir,
    )
end

function canonical_parity_table(path::String)
    isfile(path) || return DataFrame()
    source = DataFrame(CSV.File(path))
    prefix_column = :prefix_rows in propertynames(source) ?
        :prefix_rows : :prefix_raw_rows
    return DataFrame(
        motif = String.(source.motif),
        distance = Int.(source.distance),
        guide_count = Int.(source.guide_count),
        prefix_rows = Int.(source[!, prefix_column]),
        scan_rows = Int.(source.scan_rows),
        scan_only = Int.(source.scan_only),
        prefix_only = Int.(source.prefix_only),
        parity = Bool.(source.parity),
        case_dir = String.(source.case_dir),
    )
end

function replace_parity_row(table::DataFrame, row)
    nrow(table) == 0 && return DataFrame([row])
    keep = .!((table.motif .== row.motif) .&
        (table.distance .== row.distance))
    return vcat(table[keep, :], DataFrame([row]))
end

function run_parity_stage()
    output_dir = required_sweep_output_dir()
    summary_path = require_file(joinpath(output_dir, "summary.csv"))
    parity_path = joinpath(output_dir, "parity.csv")
    summary = DataFrame(CSV.File(summary_path))
    parity = canonical_parity_table(parity_path)

    for motif_name in sweep_motifs(), distance in sweep_distances()
        summary_mask =
            (String.(summary.motif) .== motif_name) .&
            (Int.(summary.distance) .== distance)
        any(summary_mask) || error(
            "Missing summary rows for motif=$motif_name distance=$distance")
        guide_count = Int(summary[findfirst(summary_mask), :guide_count])
        row = sweep_parity_row(
            output_dir, motif_name, distance, guide_count)
        parity = replace_parity_row(parity, row)
        sort!(parity, [:motif, :distance])
        summary[summary_mask, :parity] .= row.parity
        CSV.write(parity_path, parity)
        CSV.write(summary_path, summary)
        print_progress(
            "exact parity motif=$motif_name distance=$distance " *
            "status=$(row.parity ? "PASS" : "FAIL") " *
            "scan_only=$(row.scan_only) prefix_only=$(row.prefix_only)")
    end
    print_progress("Parity repair complete: $parity_path")
    return nothing
end

function run_parent()
    output_dir = sweep_output_dir()
    mkpath(output_dir)
    script = abspath(@__FILE__)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    search_threads = parse_int_env("CHOPOFF_HUMAN_SWEEP_THREADS", 24)

    for (stage, threads) in (("build", 1), ("search", search_threads))
        env = copy(ENV)
        env["CHOPOFF_HUMAN_SWEEP_STAGE"] = stage
        env["CHOPOFF_HUMAN_SWEEP_OUT"] = output_dir
        cmd = `$julia --project=$ROOT_DIR --threads=$threads $script`
        print_progress("Launching stage=$stage threads=$threads")
        run(setenv(cmd, env))
    end
    print_progress("Overnight sweep finished: $output_dir")
    return nothing
end

function sweep_main()
    stage = lowercase(strip(get(ENV, "CHOPOFF_HUMAN_SWEEP_STAGE", "parent")))
    if stage == "parent"
        run_parent()
    elseif stage == "build"
        run_build_stage()
    elseif stage == "search"
        run_search_stage()
    elseif stage == "parity"
        run_parity_stage()
    else
        error("CHOPOFF_HUMAN_SWEEP_STAGE must be parent, build, search, or parity.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    sweep_main()
end
