#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Printf
using Statistics

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const SAMPLE_DIR = joinpath(ROOT_DIR, "test", "sample_data")
const CORE_COLS = [:guide, :distance, :chromosome, :start, :strand]

function parse_int_list(name::String, default::String)
    raw = strip(get(ENV, name, default))
    return [parse(Int, strip(x)) for x in split(raw, ',') if strip(x) != ""]
end

function parse_symbol_list(name::String, default::String)
    raw = strip(get(ENV, name, default))
    return [Symbol(strip(x)) for x in split(raw, ',') if strip(x) != ""]
end

function read_core(path::String)
    df = DataFrame(CSV.File(path))
    core = select(df, CORE_COLS)
    core.guide = String.(core.guide)
    core.distance = Int.(core.distance)
    core.chromosome = String.(core.chromosome)
    core.start = Int.(core.start)
    core.strand = String.(core.strand)
    sort!(core, CORE_COLS)
    return core
end

function main()
    genome = joinpath(SAMPLE_DIR, "genome", "semirandom.fa")
    guides_all = LongDNA{4}.(sort(collect(Set(readlines(joinpath(SAMPLE_DIR, "guides.txt"))))))
    guide_counts = parse_int_list("CHOPOFF_SCAN_EXPERIMENT_GUIDES", "1,8,20")
    distances = parse_int_list("CHOPOFF_SCAN_EXPERIMENT_DISTANCES", "3")
    variants = parse_symbol_list("CHOPOFF_SCAN_EXPERIMENT_VARIANTS", "bitmask64")
    backends = parse_symbol_list(
        "CHOPOFF_SCAN_EXPERIMENT_BACKENDS",
        "legacy,fused_dict,fused_directory",
    )
    bucket_bases = parse_int_list("CHOPOFF_SCAN_EXPERIMENT_BUCKET_BASES", "9,10,11")
    verify_variants = parse_symbol_list(
        "CHOPOFF_SCAN_EXPERIMENT_VERIFY",
        "align,distance_first",
    )
    runs = only(parse_int_list("CHOPOFF_SCAN_EXPERIMENT_RUNS", "3"))
    scan_threads = only(parse_int_list(
        "CHOPOFF_SCAN_EXPERIMENT_THREADS",
        string(Threads.nthreads()),
    ))
    tdir = mktempdir(prefix = "chopoff_scan_experiment_")
    rows = NamedTuple[]

    for distance in distances
        motif = Motif("Cas9"; distance = distance)
        dbdir = joinpath(tdir, "prefixhash_d$distance")
        build_prefixHashDB("scan_experiment_d$distance", genome, motif, dbdir)
        for n in guide_counts
            guides = guides_all[1:min(n, length(guides_all))]
            prefix_out = joinpath(tdir, "prefix_d$(distance)_n$(length(guides)).csv")
            es = fill(1_000_000, distance + 1)
            search_prefixHashDB(dbdir, guides, prefix_out; distance = distance, early_stopping = es)
            prefix_core = read_core(prefix_out)
            for variant in variants
                supported_backends = distance == 3 && variant == :bitmask64 ?
                    backends : [:legacy]
                for backend in supported_backends
                    buckets = backend == :fused_directory ? bucket_bases : [first(bucket_bases)]
                    verifies = backend == :legacy ? [:align] : verify_variants
                    for bucket in buckets, verify in verifies
                        label = "$(variant)_$(backend)_b$(bucket)_$(verify)_d$(distance)_n$(length(guides))"
                        out = joinpath(tdir, "scan_" * label * ".csv")
                        kwargs = (
                            distance = distance,
                            early_stopping = es,
                            query_variant = variant,
                            scan_backend = backend,
                            bucket_bases = bucket,
                            scan_threads = scan_threads,
                            verify_variant = verify,
                        )
                        CHOPOFF.search_prefixHashScan(guides, genome, motif, out; kwargs...)

                        elapsed = Float64[]
                        for _ in 1:runs
                            push!(elapsed, @elapsed CHOPOFF.search_prefixHashScan(
                                guides, genome, motif, out; kwargs...))
                        end

                        stats = CHOPOFF.PrefixHashScanStats()
                        CHOPOFF.search_prefixHashScan(
                            guides, genome, motif, out; kwargs..., stats = stats)
                        parity = read_core(out) == prefix_core
                        elapsed_median = median(elapsed)
                        @printf(
                            "distance=%d guides=%d backend=%s bucket=%d verify=%s median=%.6fs parity=%s pairs=%d query=%.3fs scan=%.3fs rows=%d\n",
                            distance, length(guides), string(backend), bucket,
                            string(verify), elapsed_median, parity ? "PASS" : "FAIL",
                            stats.guide_pairs, Float64(stats.query_build_ns) / 1e9,
                            Float64(stats.scan_ns) / 1e9, stats.emitted_rows,
                        )
                        push!(rows, (
                            distance = distance,
                            guides = length(guides),
                            variant = string(variant),
                            backend = string(backend),
                            bucket_bases = bucket,
                            verify = string(verify),
                            scan_threads = scan_threads,
                            runs = runs,
                            elapsed_median_s = elapsed_median,
                            elapsed_min_s = minimum(elapsed),
                            elapsed_max_s = maximum(elapsed),
                            parity = parity,
                            motif_candidates = stats.motif_candidates,
                            prefix_hits = stats.prefix_hits,
                            guide_pairs = stats.guide_pairs,
                            verification_calls = stats.alignment_calls,
                            distance_calls = stats.distance_calls,
                            traceback_calls = stats.traceback_calls,
                            query_build_s = Float64(stats.query_build_ns) / 1e9,
                            chromosome_load_s = Float64(stats.chrom_load_ns) / 1e9,
                            materialize_s = Float64(stats.candidate_materialize_ns) / 1e9,
                            verify_s = Float64(stats.verify_ns) / 1e9,
                            scan_s = Float64(stats.scan_ns) / 1e9,
                            rows = stats.emitted_rows,
                        ))
                    end
                end
            end
        end
    end
    out_path = strip(get(ENV, "CHOPOFF_SCAN_EXPERIMENT_OUT", ""))
    if out_path != ""
        CSV.write(out_path, DataFrame(rows))
        println("wrote $out_path")
    end
end

main()
