#!/usr/bin/env julia
# PAMless versus Cas9 prefixHashScan on the same guides (GRCh38). Manual; not in runtests.

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Dates
using Statistics

const ROOT_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_GENOME = "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
const DEFAULT_GUIDES = joinpath(ROOT_DIR, "test", "local_human", "data", "guides_for_tests.txt")
const UNLIMITED = 1_000_000_000

env_ints(name, default) = parse.(Int, split(get(ENV, name, default), ","))
env_int(name, default) = parse(Int, get(ENV, name, string(default)))

genome = abspath(get(ENV, "CHOPOFF_PAMLESS_GENOME", DEFAULT_GENOME))
guides_path = abspath(get(ENV, "CHOPOFF_PAMLESS_GUIDES", DEFAULT_GUIDES))
out_dir = abspath(get(ENV, "CHOPOFF_PAMLESS_OUT", joinpath(ROOT_DIR, "test", "local_human",
    "outputs", "pamless_" * Dates.format(now(UTC), "yyyymmdd"))))
threads = env_int("CHOPOFF_PAMLESS_THREADS", Threads.nthreads())
count_distances = env_ints("CHOPOFF_PAMLESS_DISTANCES", "0,1,2,3,4")
detail_distances = env_ints("CHOPOFF_PAMLESS_DETAIL_DISTANCES", "0,1,2")
runs = env_int("CHOPOFF_PAMLESS_RUNS", 3)
warmups = env_int("CHOPOFF_PAMLESS_WARMUPS", 1)
mkpath(out_dir)

guides = LongDNA{4}.(filter(!isempty, strip.(readlines(guides_path))))
guide_bases = length(first(guides))
motif_kw(name) = name == "PAMless" ? (pamless = true,) : (motif = name,)
motif_obj(name, d) = name == "PAMless" ?
    CHOPOFF.prefix_hash_scan_pamless_motif(guide_bases, d) : Motif(name; distance = d)
case_path(name, d, mode) = joinpath(out_dir, "$(name)_d$(d)_$(mode).csv")
say(msg) = (println(Dates.format(now(), "HH:MM:SS"), " ", msg); flush(stdout))

function run_case(name, d, mode)
    return @elapsed search_prefixHashScan(
        guides, genome, case_path(name, d, mode); motif_kw(name)...,
        distance = d, output = mode, scan_threads = threads,
        early_stopping = fill(UNLIMITED, d + 1))
end

data_rows(path) = countlines(path) - 1

cases = [(d, mode) for d in count_distances for mode in (:counts, :detail)
    if mode == :counts || d in detail_distances]
timings = DataFrame()
stats_rows = DataFrame()
for (d, mode) in cases
    names = ["PAMless", "Cas9"]
    for _ in 1:warmups, name in names
        run_case(name, d, mode)
    end
    times = Dict(name => Float64[] for name in names)
    for r in 1:runs
        for name in (isodd(r) ? names : reverse(names))
            push!(times[name], run_case(name, d, mode))
        end
    end
    for name in names
        push!(timings, (motif = name, distance = d, output = String(mode),
            threads = threads, runs = runs, median_s = median(times[name]),
            min_s = minimum(times[name]), max_s = maximum(times[name]),
            rows = data_rows(case_path(name, d, mode))))
        say("$(name) d$(d) $(mode): median $(round(median(times[name]); digits = 3)) s")
    end
    CSV.write(joinpath(out_dir, "timings.csv"), timings)

    mode == :counts || continue
    for name in names
        stats = CHOPOFF.PrefixHashScanStats()
        search_prefixHashScan(guides, genome, motif_obj(name, d),
            joinpath(out_dir, "stats_$(name)_d$(d).csv");
            distance = d, output = :counts, scan_threads = threads,
            early_stopping = fill(UNLIMITED, d + 1), stats = stats)
        push!(stats_rows, (motif = name, distance = d,
            motif_candidates = stats.motif_candidates, prefix_hits = stats.prefix_hits,
            guide_pairs = stats.guide_pairs, path_rows = stats.path_rows,
            query_hashes = stats.query_hashes,
            query_build_s = stats.query_build_ns / 1e9))
    end
    CSV.write(joinpath(out_dir, "prefixhashscan_stats.csv"), stats_rows)
end

# Every Cas9 off-target is a PAMless off-target; only `start` may shift by a constant per strand.
function containment(d)
    cas9 = DataFrame(CSV.File(case_path("Cas9", d, :detail); types = Dict(:chromosome => String)))
    pamless = DataFrame(CSV.File(case_path("PAMless", d, :detail); types = Dict(:chromosome => String)))
    key(r, shift) = (r.guide, r.alignment_guide, r.alignment_reference, r.distance,
        r.chromosome, r.strand, r.start + shift)
    pamless_keys = Set(key(r, 0) for r in eachrow(pamless))
    rows = NamedTuple[]
    for strand in ("+", "-")
        sub = cas9[cas9.strand .== strand, :]
        shift = argmax(s -> count(r -> key(r, s) in pamless_keys, eachrow(sub)), -3:3)
        found = count(r -> key(r, shift) in pamless_keys, eachrow(sub))
        push!(rows, (distance = d, strand = strand, start_shift = shift,
            cas9_rows = nrow(sub), contained = found, missing = nrow(sub) - found))
    end
    return rows
end

contain = DataFrame(reduce(vcat, [containment(d) for d in detail_distances]))
CSV.write(joinpath(out_dir, "containment.csv"), contain)
println(timings)
println(stats_rows)
println(contain)
