#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Dates
using Printf

include(joinpath(@__DIR__, "..", "helpers", "prefix_parity.jl"))
using .PrefixParity

const ROOT_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_GENOME = "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
const DEFAULT_GUIDES = joinpath(@__DIR__, "data", "guides_for_tests.txt")
const SAMPLE_GENOME = joinpath(ROOT_DIR, "test", "sample_data", "genome", "semirandom.fa")
const SAMPLE_GUIDES = joinpath(ROOT_DIR, "test", "sample_data", "guides.txt")

function parse_bool_env(name::String, default::Bool)
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    raw in ("1", "true", "yes", "y") && return true
    raw in ("0", "false", "no", "n") && return false
    error("Invalid boolean in environment variable $name: '$raw'")
end

function parse_int_env(name::String, default::Int)
    raw = strip(get(ENV, name, string(default)))
    try
        return parse(Int, raw)
    catch
        error("Invalid integer in environment variable $name: '$raw'")
    end
end

function parse_early_stopping(distance::Int)
    raw = strip(get(ENV, "CHOPOFF_HUMAN_EARLY_STOPPING", ""))
    isempty(raw) && return fill(1_000_000, distance + 1)
    vals = parse.(Int, split(raw, ','))
    length(vals) == distance + 1 ||
        error("CHOPOFF_HUMAN_EARLY_STOPPING must have distance+1 comma-separated values")
    return vals
end

function require_file(path::String)
    isfile(path) || error("Required file not found: $path")
    return path
end

function validate_genome(path::String)
    require_file(path)
    if !endswith(lowercase(path), ".2bit")
        require_file(path * ".fai")
    end
    return path
end

function load_guides(path::String, motif::Motif)
    require_file(path)
    raw = String[]
    for line in eachline(path)
        s = uppercase(strip(line))
        isempty(s) && continue
        startswith(s, "#") && continue
        push!(raw, s)
    end
    isempty(raw) && error("No guides found in $path")

    expected_len = length_noPAM(motif)
    bad = findall(!=(expected_len), length.(raw))
    if !isempty(bad)
        first_bad = bad[1]
        error(
            "Guide length mismatch for motif $(motif.alias): expected $expected_len, " *
            "got $(length(raw[first_bad])) at guide $first_bad ($(raw[first_bad]))",
        )
    end
    return LongDNA{4}.(raw)
end

function count_result_rows(path::String)
    isfile(path) || return 0
    rows = 0
    for _ in eachline(path)
        rows += 1
    end
    return max(rows - 1, 0)
end

filesize_bytes(path::String) = isfile(path) ? stat(path).size : 0

function count_ambiguous_reference_rows(path::String)
    isfile(path) || return 0
    rows = 0
    for row in CSV.Rows(path)
        ref = getproperty(row, :alignment_reference)
        if any(c -> !(c in ('A', 'C', 'G', 'T', '-', 'a', 'c', 'g', 't')), ref)
            rows += 1
        end
    end
    return rows
end

function warmup_compile(distance::Int, backend::Symbol)
    isfile(SAMPLE_GENOME) || return
    isfile(SAMPLE_GENOME * ".fai") || return
    isfile(SAMPLE_GUIDES) || return

    motif = Motif("Cas9"; distance = min(distance, 3))
    guides = LongDNA{4}.(readlines(SAMPLE_GUIDES)[1:1])
    tmp = mktempdir(prefix = "chopoff_human_warmup_")
    try
        search_sassy(
            guides,
            SAMPLE_GENOME,
            motif,
            joinpath(tmp, "warmup.csv");
            distance = min(distance, 3),
            early_stopping = fill(10, min(distance, 3) + 1),
            backend = backend,
        )
    finally
        rm(tmp; recursive = true, force = true)
    end
end

function main()
    genome = validate_genome(String(strip(get(ENV, "CHOPOFF_HUMAN_GENOME", DEFAULT_GENOME))))
    guides_path = String(strip(get(ENV, "CHOPOFF_HUMAN_GUIDES", DEFAULT_GUIDES)))
    motif_name = String(strip(get(ENV, "CHOPOFF_HUMAN_MOTIF", "Cas9")))
    distance = parse_int_env("CHOPOFF_HUMAN_DISTANCE", 3)
    compare_prefix = parse_bool_env("CHOPOFF_HUMAN_COMPARE_PREFIX", true)
    compare_scan = parse_bool_env("CHOPOFF_HUMAN_COMPARE_SCAN", true)
    scan_query_variant = Symbol(strip(get(ENV, "CHOPOFF_HUMAN_SCAN_QUERY_VARIANT", "auto")))
    scan_backend = Symbol(strip(get(ENV, "CHOPOFF_HUMAN_SCAN_BACKEND", "auto")))
    scan_bucket_bases = parse_int_env("CHOPOFF_HUMAN_SCAN_BUCKET_BASES", 11)
    scan_threads = parse_int_env("CHOPOFF_HUMAN_SCAN_THREADS", Threads.nthreads())
    scan_prefilter_bits = parse_int_env("CHOPOFF_HUMAN_SCAN_PREFILTER_BITS", 26)
    scan_verify_variant = Symbol(strip(get(ENV, "CHOPOFF_HUMAN_SCAN_VERIFY_VARIANT", "auto")))
    rebuild_prefix = parse_bool_env("CHOPOFF_HUMAN_REBUILD_PREFIX", false)
    keep_outputs = parse_bool_env("CHOPOFF_HUMAN_KEEP_OUTPUTS", true)
    backend = Symbol(strip(get(ENV, "CHOPOFF_HUMAN_BACKEND", "auto")))

    motif = Motif(motif_name; distance = distance)
    guides = load_guides(guides_path, motif)
    early_stopping = parse_early_stopping(distance)

    stamp = Dates.format(now(Dates.UTC), dateformat"yyyymmdd_HHMMSS")
    output_dir = joinpath(@__DIR__, "outputs", stamp)
    genome_key = replace(basename(genome), r"[^A-Za-z0-9_.-]" => "_")
    index_dir = joinpath(@__DIR__, "indexes", "$(genome_key)_$(motif_name)_d$(distance)")
    mkpath(output_dir)

    println("Human SASSY benchmark")
    println("time_utc: ", Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println("genome: ", genome)
    println("guides: ", guides_path)
    println("guide_count: ", length(guides))
    println("motif: ", motif_name)
    println("distance: ", distance)
    println("threads: ", Threads.nthreads())
    println("compare_prefix: ", compare_prefix)
    println("compare_scan: ", compare_scan)
    println("scan_query_variant: ", scan_query_variant)
    println("scan_backend: ", scan_backend)
    println("scan_bucket_bases: ", scan_bucket_bases)
    println("scan_threads: ", scan_threads)
    println("scan_prefilter_bits: ", scan_prefilter_bits)
    println("scan_verify_variant: ", scan_verify_variant)
    println("sassy_backend: ", backend)
    println("sassy_backend_resolved: ", CHOPOFF.Sassy.resolve_sassy_backend(backend))
    println("output_dir: ", output_dir)

    warmup_compile(distance, backend)

    sassy_path = joinpath(output_dir, "sassy.csv")
    sassy_elapsed = @elapsed search_sassy(
        guides,
        genome,
        motif,
        sassy_path;
        distance = distance,
        early_stopping = early_stopping,
        backend = backend,
    )
    sassy_rows = count_result_rows(sassy_path)
    sassy_bytes = filesize_bytes(sassy_path)
    sassy_ambig_ref_rows = count_ambiguous_reference_rows(sassy_path)

    prefix_elapsed = missing
    prefix_rows = missing
    prefix_bytes = missing
    prefix_ambig_ref_rows = missing
    parity_sassy_only = missing
    parity_prefix_only = missing
    prefix_path = ""
    scan_elapsed = missing
    scan_rows = missing
    scan_bytes = missing
    scan_ambig_ref_rows = missing
    scan_path_source = missing
    scan_query_variant_used = missing
    scan_backend_used = missing
    scan_alignment_calls = missing
    scan_distance_calls = missing
    scan_traceback_calls = missing
    scan_metadata_s = missing
    scan_query_build_s = missing
    scan_path_load_s = missing
    scan_record_io_s = missing
    scan_sequence_convert_s = missing
    scan_scan_s = missing
    scan_findguides_s = missing
    scan_candidate_hash_s = missing
    scan_align_s = missing
    parity_scan_only = missing
    parity_prefix_only_vs_scan = missing

    if compare_prefix
        mkpath(index_dir)
        db_file = joinpath(index_dir, "prefixHashDB.bin")
        if rebuild_prefix || !isfile(db_file)
            build_prefixHashDB("human_$(motif_name)_d$(distance)", genome, motif, index_dir)
        else
            println("Reusing prefixHashDB index: ", index_dir)
        end

        prefix_path = joinpath(output_dir, "prefixhash.csv")
        prefix_elapsed = @elapsed search_prefixHashDB(
            index_dir,
            guides,
            prefix_path;
            distance = distance,
            early_stopping = early_stopping,
        )
        prefix_rows = count_result_rows(prefix_path)
        prefix_bytes = filesize_bytes(prefix_path)
        prefix_ambig_ref_rows = count_ambiguous_reference_rows(prefix_path)
        parity_sassy_only, parity_prefix_only =
            write_parity_diffs(sassy_path, prefix_path, output_dir)
    end

    if compare_scan
        scan_path = joinpath(output_dir, "prefixhashscan.csv")
        scan_stats = CHOPOFF.PrefixHashScanStats()
        scan_elapsed = @elapsed CHOPOFF.search_prefixHashScan(
            guides,
            genome,
            motif,
            scan_path;
            distance = distance,
            early_stopping = early_stopping,
            query_variant = scan_query_variant,
            scan_backend = scan_backend,
            bucket_bases = scan_bucket_bases,
            scan_threads = scan_threads,
            prefilter_bits = scan_prefilter_bits,
            verify_variant = scan_verify_variant,
            stats = scan_stats,
        )
        scan_rows = count_result_rows(scan_path)
        scan_bytes = filesize_bytes(scan_path)
        scan_ambig_ref_rows = count_ambiguous_reference_rows(scan_path)
        scan_path_source = scan_stats.path_source
        scan_query_variant_used = scan_stats.query_variant
        scan_backend_used = scan_stats.scan_backend
        scan_alignment_calls = scan_stats.alignment_calls
        scan_distance_calls = scan_stats.distance_calls
        scan_traceback_calls = scan_stats.traceback_calls
        scan_metadata_s = scan_stats.metadata_ns / 1e9
        scan_query_build_s = scan_stats.query_build_ns / 1e9
        scan_path_load_s = scan_stats.path_load_ns / 1e9
        scan_record_io_s = scan_stats.record_io_ns / 1e9
        scan_sequence_convert_s = scan_stats.sequence_convert_ns / 1e9
        scan_scan_s = scan_stats.scan_ns / 1e9
        scan_findguides_s = scan_stats.findguides_ns / 1e9
        scan_candidate_hash_s = scan_stats.candidate_hash_ns / 1e9
        scan_align_s = scan_stats.align_ns / 1e9
        scan_summary = DataFrame([scan_stats])
        CSV.write(joinpath(output_dir, "prefixhashscan_stats.csv"), scan_summary)
        if !isempty(prefix_path)
            parity_scan_only, parity_prefix_only_vs_scan =
                write_exact_parity_diffs(scan_path, prefix_path, output_dir)
        end
    end

    summary = DataFrame([(
        timestamp = Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
        genome = genome,
        guides = guides_path,
        guide_count = length(guides),
        motif = motif_name,
        distance = distance,
        threads = Threads.nthreads(),
        sassy_backend = backend,
        sassy_backend_resolved = CHOPOFF.Sassy.resolve_sassy_backend(backend),
        early_stopping = join(early_stopping, ","),
        sassy_elapsed_s = sassy_elapsed,
        sassy_rows = sassy_rows,
        sassy_bytes = sassy_bytes,
        sassy_ambig_ref_rows = sassy_ambig_ref_rows,
        compare_prefix = compare_prefix,
        compare_scan = compare_scan,
        scan_query_variant = scan_query_variant,
        scan_backend = scan_backend,
        scan_bucket_bases = scan_bucket_bases,
        scan_threads = scan_threads,
        scan_prefilter_bits = scan_prefilter_bits,
        scan_verify_variant = scan_verify_variant,
        prefix_elapsed_s = prefix_elapsed,
        prefix_rows = prefix_rows,
        prefix_bytes = prefix_bytes,
        prefix_ambig_ref_rows = prefix_ambig_ref_rows,
        parity_sassy_only = parity_sassy_only,
        parity_prefix_only = parity_prefix_only,
        scan_elapsed_s = scan_elapsed,
        scan_rows = scan_rows,
        scan_bytes = scan_bytes,
        scan_ambig_ref_rows = scan_ambig_ref_rows,
        scan_path_source = scan_path_source,
        scan_query_variant_used = scan_query_variant_used,
        scan_backend_used = scan_backend_used,
        scan_alignment_calls = scan_alignment_calls,
        scan_distance_calls = scan_distance_calls,
        scan_traceback_calls = scan_traceback_calls,
        scan_metadata_s = scan_metadata_s,
        scan_query_build_s = scan_query_build_s,
        scan_path_load_s = scan_path_load_s,
        scan_record_io_s = scan_record_io_s,
        scan_sequence_convert_s = scan_sequence_convert_s,
        scan_scan_s = scan_scan_s,
        scan_findguides_s = scan_findguides_s,
        scan_candidate_hash_s = scan_candidate_hash_s,
        scan_align_s = scan_align_s,
        parity_scan_only = parity_scan_only,
        parity_prefix_only_vs_scan = parity_prefix_only_vs_scan,
        output_dir = output_dir,
        index_dir = compare_prefix ? index_dir : "",
    )])

    summary_path = joinpath(output_dir, "summary.csv")
    CSV.write(summary_path, summary)

    @printf("SASSY elapsed: %.3fs | rows=%d | bytes=%d
", sassy_elapsed, sassy_rows, sassy_bytes)
    @printf("SASSY ambiguous-reference rows: %d
", sassy_ambig_ref_rows)
    if compare_scan
        @printf(
            "prefixHashScan elapsed: %.3fs | rows=%d | bytes=%d | ambiguous_reference_rows=%d | path_source=%s | query_variant=%s | backend=%s | verify=%s | verification_calls=%d | distance_calls=%d | tracebacks=%d
",
            scan_elapsed,
            scan_rows,
            scan_bytes,
            scan_ambig_ref_rows,
            string(scan_path_source),
            string(scan_query_variant_used),
            string(scan_backend_used),
            string(scan_verify_variant),
            scan_alignment_calls,
            scan_distance_calls,
            scan_traceback_calls,
        )
        @printf(
            "prefixHashScan timers: metadata=%.3fs | query_build=%.3fs | path_load=%.3fs | record_io=%.3fs | sequence_convert=%.3fs | scan=%.3fs | findguides=%.3fs | candidate_hash=%.3fs | align=%.3fs
",
            scan_metadata_s,
            scan_query_build_s,
            scan_path_load_s,
            scan_record_io_s,
            scan_sequence_convert_s,
            scan_scan_s,
            scan_findguides_s,
            scan_candidate_hash_s,
            scan_align_s,
        )
        if !ismissing(parity_scan_only)
            @printf(
                "Exact parity prefixHashScan vs prefixHashDB: scan_only=%s | prefix_only=%s
",
                string(parity_scan_only),
                string(parity_prefix_only_vs_scan),
            )
        end
    end
    if compare_prefix
        @printf(
            "prefixHashDB elapsed: %.3fs | rows=%d | bytes=%d | ambiguous_reference_rows=%d
",
            prefix_elapsed,
            prefix_rows,
            prefix_bytes,
            prefix_ambig_ref_rows,
        )
        @printf(
            "Core parity SASSY vs prefixHashDB: sassy_only=%s | prefix_only=%s
",
            string(parity_sassy_only),
            string(parity_prefix_only),
        )
    end
    println("summary: ", summary_path)

    if !keep_outputs
        rm(output_dir; recursive = true, force = true)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
