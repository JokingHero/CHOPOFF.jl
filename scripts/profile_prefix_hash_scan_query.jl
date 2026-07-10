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

function parse_variants()
    raw = strip(get(ENV, "CHOPOFF_QUERY_PROFILE_VARIANTS", "baseline,columnwise,bitmask64"))
    variants = Symbol[]
    for item in split(raw, ',')
        name = Symbol(strip(item))
        name in (:baseline, :columnwise, :bitmask64, :bruteforce) || error("Unknown query variant: $name")
        push!(variants, name)
    end
    isempty(variants) && error("No query variants selected.")
    return variants
end

function require_file(path::String)
    isfile(path) || error("Required file not found: $path")
    return path
end

function selected_case_labels()
    raw = lowercase(strip(get(ENV, "CHOPOFF_QUERY_PROFILE_CASES", "cas9,cas12a")))
    labels = Set{String}()
    for item in split(raw, ',')
        label = strip(item)
        label in ("cas9", "cas12a") || error("Unknown query profile case: $label")
        push!(labels, label)
    end
    isempty(labels) && error("No query profile cases selected.")
    return labels
end

function load_cases(distance::Int)
    labels = selected_case_labels()
    cases = NamedTuple[]
    if "cas9" in labels
        cas9_guides_path = require_file(joinpath(SAMPLE_DIR, "guides.txt"))
        cas9_guides = LongDNA{4}.(sort(collect(Set(readlines(cas9_guides_path)))))
        push!(cases, (label = "Cas9", motif = Motif("Cas9"; distance = distance), guides = cas9_guides))
    end
    if "cas12a" in labels
        push!(cases, (label = "Cas12a", motif = Motif("Cas12a"; distance = distance), guides = copy(CAS12A_GUIDES)))
    end
    return cases
end

function hash_type_for(hash_len::Int)
    return CHOPOFF.smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
end

function summarize(values::Vector{Float64})
    return (median = median(values), mean = mean(values), min = minimum(values), max = maximum(values))
end

function normalize_query_map(query::Dict)
    out = Dict{eltype(keys(query)), Vector{Int}}()
    for (k, v) in query
        out[k] = sort(v)
    end
    return out
end

function normalize_query_map(query::CHOPOFF.PrefixHashScanBitmaskQuery)
    out = Dict{eltype(keys(query.masks)), Vector{Int}}()
    for (k, mask0) in query.masks
        mask = mask0
        guides = Int[]
        while mask != 0
            push!(guides, trailing_zeros(mask) + 1)
            mask &= mask - 1
        end
        out[k] = guides
    end
    return out
end

function profile_case(case, distance::Int, hash_len::Int, runs::Int, variants::Vector{Symbol})
    println("\n=== Query profile: $(case.label) | distance=$distance | guides=$(length(case.guides)) | hash_len=$hash_len ===")
    load_stats = CHOPOFF.PrefixHashScanStats()
    paths, source = CHOPOFF.load_prefix_hash_scan_paths(case.motif, distance, hash_len, load_stats)
    guides_ = CHOPOFF.oriented_prefix_hash_scan_guides(case.guides, case.motif)
    hash_type = hash_type_for(hash_len)
    @printf("paths=%d source=%s load=%.3fs\n", size(paths, 1), string(source), Float64(load_stats.path_load_ns) / 1e9)

    baseline_query = nothing
    rows = NamedTuple[]
    for variant in variants
        if variant == :bruteforce
            warm_stats = CHOPOFF.PrefixHashScanStats()
            warm_query = nothing
        else
            warm_stats = CHOPOFF.PrefixHashScanStats()
            warm_query = CHOPOFF.build_prefix_hash_scan_map_from_paths(
                paths,
                guides_,
                hash_type,
                warm_stats;
                query_variant = variant,
            )
            if variant == :baseline
                baseline_query = normalize_query_map(warm_query)
            elseif baseline_query !== nothing
                normalize_query_map(warm_query) == baseline_query || error("Variant $variant does not match baseline query map for $(case.label).")
            end
        end

        elapsed = Float64[]
        allocated = Float64[]
        last_stats = nothing
        for run_id in 1:runs
            stats = CHOPOFF.PrefixHashScanStats()
            GC.gc()
            if variant == :bruteforce
                stats.query_variant = :bruteforce
                push!(elapsed, 0.0)
                push!(allocated, 0.0)
            else
                push!(elapsed, @elapsed CHOPOFF.build_prefix_hash_scan_map_from_paths(
                    paths,
                    guides_,
                    hash_type,
                    stats;
                    query_variant = variant,
                ))
                GC.gc()
                push!(allocated, Float64(@allocated CHOPOFF.build_prefix_hash_scan_map_from_paths(
                    paths,
                    guides_,
                    hash_type,
                    nothing;
                    query_variant = variant,
                )))
            end
            last_stats = stats
        end
        estats = summarize(elapsed)
        astats = summarize(allocated)
        @printf(
            "%s | median=%.6fs mean=%.6fs min=%.6fs max=%.6fs alloc_median=%.0f | hashes=%d format=%.3fs fold=%.3fs dedup=%.3fs insert=%.3fs\n",
            string(variant),
            estats.median,
            estats.mean,
            estats.min,
            estats.max,
            astats.median,
            last_stats.query_hashes,
            Float64(last_stats.query_format_ns) / 1e9,
            Float64(last_stats.query_fold_ns) / 1e9,
            Float64(last_stats.query_dedup_ns) / 1e9,
            Float64(last_stats.query_insert_ns) / 1e9,
        )
        push!(rows, (
            timestamp = Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
            motif = case.label,
            distance = distance,
            hash_len = hash_len,
            guides = length(case.guides),
            paths = size(paths, 1),
            source = string(source),
            variant = string(variant),
            runs = runs,
            median_s = estats.median,
            mean_s = estats.mean,
            min_s = estats.min,
            max_s = estats.max,
            alloc_median = astats.median,
            query_hashes = last_stats.query_hashes,
            format_s = Float64(last_stats.query_format_ns) / 1e9,
            fold_s = Float64(last_stats.query_fold_ns) / 1e9,
            dedup_s = Float64(last_stats.query_dedup_ns) / 1e9,
            insert_s = Float64(last_stats.query_insert_ns) / 1e9,
        ))
    end
    return rows
end

function main()
    distance = parse_int_env("CHOPOFF_QUERY_PROFILE_DISTANCE", 3)
    hash_len = parse_int_env("CHOPOFF_QUERY_PROFILE_HASH_LEN", min(20 - distance, 16))
    runs = parse_int_env("CHOPOFF_QUERY_PROFILE_RUNS", 7)
    variants = parse_variants()
    out_path = strip(get(ENV, "CHOPOFF_QUERY_PROFILE_OUT", ""))

    println("prefixHashScan query profile")
    println("time_utc: ", Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println("threads: ", Threads.nthreads())
    println("runs: $runs")
    println("cases: ", get(ENV, "CHOPOFF_QUERY_PROFILE_CASES", "cas9,cas12a"))
    println("variants: ", join(string.(variants), ", "))

    rows = NamedTuple[]
    for case in load_cases(distance)
        append!(rows, profile_case(case, distance, hash_len, runs, variants))
    end

    if out_path != ""
        mkpath(dirname(out_path))
        CSV.write(out_path, DataFrame(rows))
        println("\nWrote query profile CSV: $out_path")
    end
end

main()
