#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Printf
using Statistics

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_GENOME = "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
const DEFAULT_GUIDES = joinpath(ROOT_DIR, "test", "local_human", "data", "guides_for_tests.txt")

struct TuningConfig
    chunk_bases::Int
    prefilter_bits::Int
    bucket_bases::Int
end

function parse_int_env(name::String, default::Int)
    return parse(Int, strip(get(ENV, name, string(default))))
end

function parse_int_list(name::String, default::String)
    return parse.(Int, filter(!isempty, strip.(split(get(ENV, name, default), ","))))
end

function parse_symbol_list(name::String, default::String)
    return Symbol.(filter(!isempty, strip.(split(get(ENV, name, default), ","))))
end

function load_guides(path::String)
    return LongDNA{4}.(filter(!isempty, strip.(readlines(path))))
end

function config_label(config::TuningConfig)
    return "c$(config.chunk_bases)_p$(config.prefilter_bits)_b$(config.bucket_bases)"
end

function parse_configs(raw::String)
    configs = TuningConfig[]
    for item in filter(!isempty, strip.(split(raw, ',')))
        values = parse.(Int, split(item, ':'))
        length(values) == 3 || error("Expected chunk:prefilter:bucket, got '$item'.")
        push!(configs, TuningConfig(values...))
    end
    return configs
end

function result_signature(result)
    results, ranges = result
    values = Tuple[]
    for chrom_idx in eachindex(ranges), strand in (:plus, :minus)
        for result_idx in ranges[chrom_idx]
            for hit in getfield(something(results[result_idx]), strand)
                push!(values, (
                    chrom_idx, strand, hit.guide_idx, hit.pos, hit.dist,
                    hit.is_antisense, hit.aln_guide, hit.aln_ref))
            end
        end
    end
    return values
end

function build_prepared(
    config::TuningConfig, guides, motif;
    query_build_backend::Symbol = :serial, query_threads::Int = 1)

    query = nothing
    guides_ = nothing
    elapsed = @elapsed query, guides_ = CHOPOFF.build_prefix_hash_scan_compact_query(
        guides, motif, 3, 16, nothing;
        bucket_bases = config.bucket_bases,
        prefilter_bits = config.prefilter_bits,
        query_build_backend = query_build_backend,
        query_threads = query_threads,
    )
    profiles = CHOPOFF.build_prefix_hash_scan_myers_profiles(guides_)
    return query, guides_, profiles, elapsed
end

function scan_prepared(
    genome, reference_lengths, dbi, prepared, config, scan_threads;
    lookup_variant::Symbol = :auto)

    query, guides_, profiles = prepared
    use_bucketed = lookup_variant == :bucketed ||
        (lookup_variant == :auto && config.prefilter_bits != 0 &&
            config.bucket_bases == 11 &&
            config.chunk_bases + 25 <= typemax(UInt32))
    mode = use_bucketed ? :bucketed_reuse : :buffered_reuse
    return CHOPOFF.stream_prefix_hash_scan(
        genome, reference_lengths, query, dbi, guides_, profiles, 3,
        config.chunk_bases, scan_threads, Val(mode), nothing, Val(:chunk))
end

function search_kwargs(
    config::TuningConfig, scan_threads::Int;
    query_build_backend::Symbol = :auto,
    lookup_variant::Symbol = :auto)

    return (
        distance = 3,
        early_stopping = fill(1_000_000, 4),
        query_variant = :bitmask64,
        scan_backend = :streaming_fasta_simd,
        bucket_bases = config.bucket_bases,
        scan_threads = scan_threads,
        stream_chunk_bases = config.chunk_bases,
        prefilter_bits = config.prefilter_bits,
        query_build_backend = query_build_backend,
        lookup_variant = lookup_variant,
        verify_variant = :myers_raw,
    )
end

function stage_configs(stage::Symbol)
    fixed_chunk = parse_int_env("CHOPOFF_TUNING_FIXED_CHUNK", 2 * 1024 * 1024)
    fixed_prefilter = parse_int_env("CHOPOFF_TUNING_FIXED_PREFILTER", 26)
    fixed_bucket = parse_int_env("CHOPOFF_TUNING_FIXED_BUCKET", 11)
    if stage == :chunk
        return [TuningConfig(x, fixed_prefilter, fixed_bucket) for x in
            parse_int_list("CHOPOFF_TUNING_CHUNKS", "2097152,4194304,8388608,16777216")]
    elseif stage == :prefilter
        return [TuningConfig(fixed_chunk, x, fixed_bucket) for x in
            parse_int_list("CHOPOFF_TUNING_PREFILTERS", "22,24,26")]
    elseif stage == :bucket
        return [TuningConfig(fixed_chunk, fixed_prefilter, x) for x in
            parse_int_list("CHOPOFF_TUNING_BUCKETS", "9,10,11,12")]
    elseif stage in (:query, :lookup)
        return [TuningConfig(fixed_chunk, fixed_prefilter, fixed_bucket)]
    elseif stage == :final
        return parse_configs(get(ENV, "CHOPOFF_TUNING_CONFIGS", "2097152:26:11"))
    end
    error("CHOPOFF_TUNING_STAGE must be chunk, prefilter, bucket, query, lookup, or final.")
end

function run_query_stage(
    genome, reference_lengths, dbi, guides, motif, config, output_dir,
    runs::Int, warmups::Int, allocation_runs::Int, scan_threads::Int)

    backends = parse_symbol_list(
        "CHOPOFF_TUNING_QUERY_BACKENDS", "serial,parallel")
    all(in((:serial, :parallel)), backends) ||
        error("CHOPOFF_TUNING_QUERY_BACKENDS must contain serial or parallel.")
    :serial in backends || pushfirst!(backends, :serial)
    for backend in backends, _ in 1:warmups
        query, guides_, profiles, _ = build_prepared(
            config, guides, motif; query_build_backend = backend,
            query_threads = scan_threads)
        scan_prepared(
            genome, reference_lengths, dbi, (query, guides_, profiles),
            config, scan_threads)
    end

    baseline_built = build_prepared(
        config, guides, motif; query_build_backend = :serial,
        query_threads = scan_threads)
    baseline_prepared = baseline_built[1:3]
    baseline_signature = result_signature(scan_prepared(
        genome, reference_lengths, dbi, baseline_prepared, config, scan_threads))
    baseline_out = joinpath(output_dir, "query_serial_threads$(scan_threads).csv")
    CHOPOFF.search_prefixHashScan(
        guides, genome, motif, baseline_out;
        search_kwargs(config, scan_threads; query_build_backend = :serial)...)
    baseline_bytes = read(baseline_out)
    baseline_stats = CHOPOFF.PrefixHashScanStats()
    CHOPOFF.search_prefixHashScan(
        guides, genome, motif, baseline_out;
        search_kwargs(config, scan_threads; query_build_backend = :serial)...,
        stats = baseline_stats)
    semantic_fields = (
        :path_rows, :query_hashes, :motif_candidates, :prefix_hits,
        :guide_pairs, :alignment_calls, :distance_calls, :traceback_calls,
        :emitted_rows)

    rows = NamedTuple[]
    for round in 1:runs
        order = circshift(backends, -(round - 1) % length(backends))
        for backend in order
            GC.gc()
            query = nothing
            guides_ = nothing
            profiles = nothing
            query_elapsed = @elapsed query, guides_, profiles, _ = build_prepared(
                config, guides, motif; query_build_backend = backend,
                query_threads = scan_threads)
            prepared = (query, guides_, profiles)
            scan_result = nothing
            GC.gc()
            scan_elapsed = @elapsed scan_result = scan_prepared(
                genome, reference_lengths, dbi, prepared, config, scan_threads)
            signature_parity = result_signature(scan_result) == baseline_signature
            out = joinpath(output_dir,
                "query_$(backend)_threads$(scan_threads).csv")
            GC.gc()
            full_elapsed = @elapsed CHOPOFF.search_prefixHashScan(
                guides, genome, motif, out;
                search_kwargs(config, scan_threads;
                    query_build_backend = backend)...)
            byte_parity = read(out) == baseline_bytes
            push!(rows, (
                stage = "query", threads = scan_threads, round = round,
                guides = length(guides), query_build_backend = string(backend),
                query_build_s = query_elapsed, prepared_scan_s = scan_elapsed,
                end_to_end_s = full_elapsed, signature_parity = signature_parity,
                byte_parity = byte_parity, rows = length(baseline_signature),
            ))
            @printf(
                "stage=query threads=%d guides=%d round=%d backend=%s query=%.6f scan=%.6f full=%.6f parity=%s\n",
                scan_threads, length(guides), round, string(backend),
                query_elapsed, scan_elapsed, full_elapsed,
                signature_parity && byte_parity ? "PASS" : "FAIL")
        end
    end

    allocation_rows = NamedTuple[]
    diagnostic_rows = NamedTuple[]
    for backend in backends
        out = joinpath(output_dir,
            "query_allocation_$(backend)_threads$(scan_threads).csv")
        for run in 1:allocation_runs
            GC.gc()
            allocated = @allocated CHOPOFF.search_prefixHashScan(
                guides, genome, motif, out;
                search_kwargs(config, scan_threads;
                    query_build_backend = backend)...)
            push!(allocation_rows, (
                query_build_backend = string(backend), run = run,
                allocated_bytes = allocated, byte_parity = read(out) == baseline_bytes))
        end
        stats = CHOPOFF.PrefixHashScanStats()
        CHOPOFF.search_prefixHashScan(
            guides, genome, motif, out;
            search_kwargs(config, scan_threads; query_build_backend = backend)...,
            stats = stats)
        push!(diagnostic_rows, (
            query_build_backend = string(backend),
            semantic_parity = all(field ->
                getfield(stats, field) == getfield(baseline_stats, field),
                semantic_fields),
            query_hashes = stats.query_hashes, guide_pairs = stats.guide_pairs,
            emitted_rows = stats.emitted_rows,
            query_build_s_stats = Float64(stats.query_build_ns) / 1e9,
            scan_s_stats = Float64(stats.scan_ns) / 1e9,
        ))
    end

    keys = [:query_build_backend]
    summary = combine(groupby(DataFrame(rows), keys),
        :query_build_s => median => :query_build_median_s,
        :prepared_scan_s => median => :prepared_scan_median_s,
        :end_to_end_s => median => :end_to_end_median_s,
        :signature_parity => all => :signature_parity,
        :byte_parity => all => :byte_parity)
    allocations = combine(groupby(DataFrame(allocation_rows), keys),
        :allocated_bytes => median => :allocated_median_bytes,
        :byte_parity => all => :allocation_byte_parity)
    summary = leftjoin(summary, allocations; on = keys)
    summary = leftjoin(summary, DataFrame(diagnostic_rows); on = keys)
    insertcols!(summary, 1, :stage => fill("query", nrow(summary)))
    insertcols!(summary, 2, :threads => fill(scan_threads, nrow(summary)))
    insertcols!(summary, 3, :guides => fill(length(guides), nrow(summary)))
    prefix = "query_guides$(length(guides))_threads$(scan_threads)"
    CSV.write(joinpath(output_dir, prefix * "_runs.csv"), DataFrame(rows))
    CSV.write(joinpath(output_dir, prefix * "_allocations.csv"),
        DataFrame(allocation_rows))
    CSV.write(joinpath(output_dir, prefix * "_diagnostics.csv"),
        DataFrame(diagnostic_rows))
    summary_path = joinpath(output_dir, prefix * "_summary.csv")
    CSV.write(summary_path, summary)
    println(summary)
    println("summary: $summary_path")
    return nothing
end

function run_lookup_stage(
    genome, reference_lengths, dbi, guides, motif, config, output_dir,
    runs::Int, warmups::Int, allocation_runs::Int, scan_threads::Int)

    variants = parse_symbol_list(
        "CHOPOFF_TUNING_LOOKUP_VARIANTS", "inline,bucketed")
    all(in((:inline, :bucketed)), variants) ||
        error("CHOPOFF_TUNING_LOOKUP_VARIANTS must contain inline or bucketed.")
    :inline in variants || pushfirst!(variants, :inline)

    built = build_prepared(
        config, guides, motif; query_build_backend = :parallel,
        query_threads = scan_threads)
    prepared = built[1:3]
    for variant in variants, _ in 1:warmups
        scan_prepared(
            genome, reference_lengths, dbi, prepared, config, scan_threads;
            lookup_variant = variant)
    end

    baseline_signature = result_signature(scan_prepared(
        genome, reference_lengths, dbi, prepared, config, scan_threads;
        lookup_variant = :inline))
    baseline_out = joinpath(
        output_dir, "lookup_inline_threads$(scan_threads).csv")
    CHOPOFF.search_prefixHashScan(
        guides, genome, motif, baseline_out;
        search_kwargs(config, scan_threads; lookup_variant = :inline)...)
    baseline_bytes = read(baseline_out)
    baseline_stats = CHOPOFF.PrefixHashScanStats()
    CHOPOFF.search_prefixHashScan(
        guides, genome, motif, baseline_out;
        search_kwargs(config, scan_threads; lookup_variant = :inline)...,
        stats = baseline_stats)
    semantic_fields = (
        :path_rows, :query_hashes, :motif_candidates, :prefix_hits,
        :guide_pairs, :alignment_calls, :distance_calls, :traceback_calls,
        :emitted_rows)

    rows = NamedTuple[]
    for round in 1:runs
        order = circshift(variants, -(round - 1) % length(variants))
        for variant in order
            GC.gc()
            scan_result = nothing
            scan_elapsed = @elapsed scan_result = scan_prepared(
                genome, reference_lengths, dbi, prepared, config, scan_threads;
                lookup_variant = variant)
            signature_parity =
                result_signature(scan_result) == baseline_signature
            out = joinpath(
                output_dir, "lookup_$(variant)_threads$(scan_threads).csv")
            GC.gc()
            full_elapsed = @elapsed CHOPOFF.search_prefixHashScan(
                guides, genome, motif, out;
                search_kwargs(config, scan_threads;
                    lookup_variant = variant)...)
            byte_parity = read(out) == baseline_bytes
            push!(rows, (
                stage = "lookup", threads = scan_threads, round = round,
                guides = length(guides), lookup_variant = string(variant),
                prepared_scan_s = scan_elapsed,
                end_to_end_s = full_elapsed,
                signature_parity = signature_parity,
                byte_parity = byte_parity, rows = length(baseline_signature),
            ))
            @printf(
                "stage=lookup threads=%d round=%d variant=%s scan=%.6f full=%.6f parity=%s\n",
                scan_threads, round, string(variant), scan_elapsed,
                full_elapsed,
                signature_parity && byte_parity ? "PASS" : "FAIL")
        end
    end

    allocation_rows = NamedTuple[]
    diagnostic_rows = NamedTuple[]
    for variant in variants
        out = joinpath(
            output_dir, "lookup_allocation_$(variant)_threads$(scan_threads).csv")
        for run in 1:allocation_runs
            GC.gc()
            allocated = @allocated CHOPOFF.search_prefixHashScan(
                guides, genome, motif, out;
                search_kwargs(config, scan_threads;
                    lookup_variant = variant)...)
            push!(allocation_rows, (
                lookup_variant = string(variant), run = run,
                allocated_bytes = allocated,
                byte_parity = read(out) == baseline_bytes))
        end
        stats = CHOPOFF.PrefixHashScanStats()
        CHOPOFF.search_prefixHashScan(
            guides, genome, motif, out;
            search_kwargs(config, scan_threads; lookup_variant = variant)...,
            stats = stats)
        push!(diagnostic_rows, (
            lookup_variant = string(variant),
            semantic_parity = all(field ->
                getfield(stats, field) == getfield(baseline_stats, field),
                semantic_fields),
            motif_candidates = stats.motif_candidates,
            prefix_hits = stats.prefix_hits, guide_pairs = stats.guide_pairs,
            emitted_rows = stats.emitted_rows,
            query_build_s_stats = Float64(stats.query_build_ns) / 1e9,
            scan_s_stats = Float64(stats.scan_ns) / 1e9,
        ))
    end

    keys = [:lookup_variant]
    summary = combine(groupby(DataFrame(rows), keys),
        :prepared_scan_s => median => :prepared_scan_median_s,
        :end_to_end_s => median => :end_to_end_median_s,
        :signature_parity => all => :signature_parity,
        :byte_parity => all => :byte_parity)
    allocations = combine(groupby(DataFrame(allocation_rows), keys),
        :allocated_bytes => median => :allocated_median_bytes,
        :byte_parity => all => :allocation_byte_parity)
    summary = leftjoin(summary, allocations; on = keys)
    summary = leftjoin(summary, DataFrame(diagnostic_rows); on = keys)
    insertcols!(summary, 1, :stage => fill("lookup", nrow(summary)))
    insertcols!(summary, 2, :threads => fill(scan_threads, nrow(summary)))
    insertcols!(summary, 3, :guides => fill(length(guides), nrow(summary)))
    prefix = "lookup_guides$(length(guides))_threads$(scan_threads)"
    CSV.write(joinpath(output_dir, prefix * "_runs.csv"), DataFrame(rows))
    CSV.write(joinpath(output_dir, prefix * "_allocations.csv"),
        DataFrame(allocation_rows))
    CSV.write(joinpath(output_dir, prefix * "_diagnostics.csv"),
        DataFrame(diagnostic_rows))
    summary_path = joinpath(output_dir, prefix * "_summary.csv")
    CSV.write(summary_path, summary)
    println(summary)
    println("summary: $summary_path")
    return nothing
end

function main()
    genome = abspath(get(ENV, "CHOPOFF_TUNING_GENOME", DEFAULT_GENOME))
    guides_path = abspath(get(ENV, "CHOPOFF_TUNING_GUIDES", DEFAULT_GUIDES))
    output_dir = abspath(get(ENV, "CHOPOFF_TUNING_OUT", "/tmp/prefix_hash_scan_tuning"))
    stage = Symbol(lowercase(strip(get(ENV, "CHOPOFF_TUNING_STAGE", "chunk"))))
    runs = parse_int_env(
        "CHOPOFF_TUNING_RUNS", stage in (:final, :query, :lookup) ? 15 : 7)
    warmups = parse_int_env("CHOPOFF_TUNING_WARMUPS", 2)
    allocation_runs = parse_int_env(
        "CHOPOFF_TUNING_ALLOCATION_RUNS",
        stage == :lookup ? 5 : (stage in (:final, :query) ? 3 : 1))
    scan_threads = parse_int_env("CHOPOFF_TUNING_THREADS", Threads.nthreads())
    configs = stage_configs(stage)
    baseline = TuningConfig(2 * 1024 * 1024, 26, 11)
    stage == :final && baseline ∉ configs && pushfirst!(configs, baseline)
    isempty(configs) && error("No tuning configurations selected.")
    isfile(genome) || error("Genome not found: $genome")
    isfile(genome * ".fai") || error("FASTA index not found: $(genome).fai")
    isfile(guides_path) || error("Guides not found: $guides_path")
    mkpath(output_dir)

    guides = load_guides(guides_path)
    guide_limit = parse_int_env("CHOPOFF_TUNING_GUIDE_LIMIT", 0)
    guide_limit > 0 && (guides = guides[1:min(guide_limit, length(guides))])
    motif = Motif("Cas9"; distance = 3)
    dbi, reference_lengths = CHOPOFF.prefix_hash_scan_dbinfo(genome, motif)
    if stage == :query
        run_query_stage(
            genome, reference_lengths, dbi, guides, motif, only(configs),
            output_dir, runs, warmups, allocation_runs, scan_threads)
        return
    elseif stage == :lookup
        run_lookup_stage(
            genome, reference_lengths, dbi, guides, motif, only(configs),
            output_dir, runs, warmups, allocation_runs, scan_threads)
        return
    end
    prepared = Dict{TuningConfig, Tuple}()
    for config in configs
        for _ in 1:warmups
            query, guides_, profiles, _ = build_prepared(config, guides, motif)
            prepared[config] = (query, guides_, profiles)
            scan_prepared(genome, reference_lengths, dbi, prepared[config], config, scan_threads)
        end
    end

    baseline_prepared = if haskey(prepared, baseline)
        prepared[baseline]
    else
        build_prepared(baseline, guides, motif)[1:3]
    end
    baseline_signature = result_signature(scan_prepared(
        genome, reference_lengths, dbi, baseline_prepared, baseline, scan_threads))
    baseline_out = joinpath(output_dir, "baseline_threads$(scan_threads).csv")
    CHOPOFF.search_prefixHashScan(
        guides, genome, motif, baseline_out; search_kwargs(baseline, scan_threads)...)
    baseline_stats = CHOPOFF.PrefixHashScanStats()
    CHOPOFF.search_prefixHashScan(
        guides, genome, motif, baseline_out;
        search_kwargs(baseline, scan_threads)..., stats = baseline_stats)
    semantic_fields = (
        :motif_candidates, :prefix_hits, :guide_pairs, :alignment_calls,
        :distance_calls, :traceback_calls, :emitted_rows)

    baseline_bytes = read(baseline_out)

    rows = NamedTuple[]
    for round in 1:runs
        order = circshift(configs, -(round - 1) % length(configs))
        for config in order
            GC.gc()
            query = nothing
            guides_ = nothing
            profiles = nothing
            query_elapsed = @elapsed query, guides_, profiles, _ =
                build_prepared(config, guides, motif)
            prepared[config] = (query, guides_, profiles)

            scan_result = nothing
            GC.gc()
            scan_elapsed = @elapsed scan_result = scan_prepared(
                genome, reference_lengths, dbi, prepared[config], config, scan_threads)
            signature_parity = result_signature(scan_result) == baseline_signature

            out = joinpath(output_dir, "$(config_label(config))_threads$(scan_threads).csv")
            GC.gc()
            end_to_end = @elapsed CHOPOFF.search_prefixHashScan(
                guides, genome, motif, out; search_kwargs(config, scan_threads)...)
            byte_parity = read(out) == baseline_bytes
            push!(rows, (
                stage = string(stage), threads = scan_threads, round = round,
                chunk_bases = config.chunk_bases,
                prefilter_bits = config.prefilter_bits,
                bucket_bases = config.bucket_bases,
                query_build_s = query_elapsed,
                prepared_scan_s = scan_elapsed,
                end_to_end_s = end_to_end,
                signature_parity = signature_parity,
                byte_parity = byte_parity,
                rows = length(baseline_signature),
            ))
            @printf(
                "stage=%s threads=%d round=%d config=%s query=%.6f scan=%.6f full=%.6f parity=%s\n",
                string(stage), scan_threads, round, config_label(config),
                query_elapsed, scan_elapsed, end_to_end,
                signature_parity && byte_parity ? "PASS" : "FAIL")
        end
    end

    allocation_rows = NamedTuple[]
    for config in configs, run in 1:allocation_runs
        out = joinpath(output_dir, "allocation_$(config_label(config))_threads$(scan_threads).csv")
        GC.gc()
        allocated = @allocated CHOPOFF.search_prefixHashScan(
            guides, genome, motif, out; search_kwargs(config, scan_threads)...)
        push!(allocation_rows, (
            stage = string(stage), threads = scan_threads, run = run,
            chunk_bases = config.chunk_bases,
            prefilter_bits = config.prefilter_bits,
            bucket_bases = config.bucket_bases,
            allocated_bytes = allocated,
            byte_parity = read(out) == baseline_bytes,
        ))
    end
    diagnostic_rows = NamedTuple[]
    for config in configs
        out = joinpath(output_dir, "diagnostic_$(config_label(config))_threads$(scan_threads).csv")
        stats = CHOPOFF.PrefixHashScanStats()
        CHOPOFF.search_prefixHashScan(
            guides, genome, motif, out;
            search_kwargs(config, scan_threads)..., stats = stats)
        semantic_parity = all(field ->
            getfield(stats, field) == getfield(baseline_stats, field), semantic_fields)
        push!(diagnostic_rows, (
            chunk_bases = config.chunk_bases,
            prefilter_bits = config.prefilter_bits,
            bucket_bases = config.bucket_bases,
            semantic_parity = semantic_parity,
            motif_candidates = stats.motif_candidates,
            prefix_hits = stats.prefix_hits,
            guide_pairs = stats.guide_pairs,
            alignment_calls = stats.alignment_calls,
            distance_calls = stats.distance_calls,
            traceback_calls = stats.traceback_calls,
            emitted_rows = stats.emitted_rows,
            query_build_s_stats = Float64(stats.query_build_ns) / 1e9,
            scan_s_stats = Float64(stats.scan_ns) / 1e9,
        ))
    end


    runs_path = joinpath(output_dir, "$(stage)_threads$(scan_threads)_runs.csv")
    allocations_path = joinpath(output_dir, "$(stage)_threads$(scan_threads)_allocations.csv")
    diagnostics_path = joinpath(output_dir, "$(stage)_threads$(scan_threads)_diagnostics.csv")
    CSV.write(runs_path, DataFrame(rows))
    CSV.write(allocations_path, DataFrame(allocation_rows))
    CSV.write(diagnostics_path, DataFrame(diagnostic_rows))

    summary = combine(
        groupby(DataFrame(rows), [:chunk_bases, :prefilter_bits, :bucket_bases]),
        :query_build_s => median => :query_build_median_s,
        :prepared_scan_s => median => :prepared_scan_median_s,
        :end_to_end_s => median => :end_to_end_median_s,
        :signature_parity => all => :signature_parity,
        :byte_parity => all => :byte_parity,
    )
    allocation_summary = combine(
        groupby(DataFrame(allocation_rows), [:chunk_bases, :prefilter_bits, :bucket_bases]),
        :allocated_bytes => median => :allocated_median_bytes,
        :byte_parity => all => :allocation_byte_parity,
    )
    summary = leftjoin(summary, allocation_summary;
        on = [:chunk_bases, :prefilter_bits, :bucket_bases])
    summary = leftjoin(summary, DataFrame(diagnostic_rows);
        on = [:chunk_bases, :prefilter_bits, :bucket_bases])
    insertcols!(summary, 1, :stage => fill(string(stage), nrow(summary)))
    insertcols!(summary, 2, :threads => fill(scan_threads, nrow(summary)))
    summary_path = joinpath(output_dir, "$(stage)_threads$(scan_threads)_summary.csv")
    CSV.write(summary_path, summary)
    println(summary)
    println("runs: $runs_path")
    println("summary: $summary_path")
end

main()
