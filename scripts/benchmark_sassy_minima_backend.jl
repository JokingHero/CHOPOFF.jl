#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Random
using Statistics

const TEST_DIR = mktempdir()

function build_genome(genome_seq::String; chrom::String = "chr1")
    genome_path = joinpath(TEST_DIR, "bench_genome.fa")
    len = length(genome_seq)
    open(genome_path, "w") do f
        write(f, ">$chrom\n$genome_seq\n")
    end
    header_len = length(">$chrom\n")
    open(genome_path * ".fai", "w") do f
        write(f, "$chrom\t$len\t$header_len\t$len\t$(len + 1)\n")
    end
    return genome_path
end

function run_once(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif;
    distance::Int,
    force_safe_minima::Bool,
)
    output_path = joinpath(TEST_DIR, "bench_$(randstring(10)).csv")
    elapsed = @elapsed CHOPOFF.search_sassy(
        guides,
        genome_path,
        motif,
        output_path;
        distance = distance,
        force_safe_minima = force_safe_minima,
    )
    if isfile(output_path)
        rm(output_path; force = true)
    end
    return elapsed
end

function main()
    runs = parse(Int, get(ENV, "SASSY_BENCH_RUNS", "7"))
    enforce = get(ENV, "SASSY_BENCH_ENFORCE_NO_REGRESSION", "0") == "1"

    pad = repeat("TACG", 64)
    guide1 = "TTTTTTTTTTTTTTTTTTTT"
    guide2 = "ACGTACGTACGTACGTACGT"
    guides = [LongDNA{4}(guide1), LongDNA{4}(guide2)]

    hit1 = guide1 * "AGG"
    hit2 = "ATTTTTTTTTTTTTTTTTTT" * "AGG"
    hit3 = guide2 * "TGG"
    hit4 = "ACGTACGTACGTACGTACGA" * "AGG"
    genome_seq = pad * hit1 * pad * hit2 * pad * hit3 * pad * hit4 * pad
    genome_path = build_genome(genome_seq)
    motif = Motif("Cas9"; distance = 3)

    println("Sassy minima backend benchmark")
    println("CPU architecture: $(Sys.ARCH)")
    println("BMI2 available: $(CHOPOFF.Sassy.can_use_bmi2_pext())")
    println("Runs per backend: $runs")

    # Warmup compile + I/O once per backend.
    run_once(guides, genome_path, motif; distance = 3, force_safe_minima = true)
    run_once(guides, genome_path, motif; distance = 3, force_safe_minima = false)

    safe_times = Float64[]
    auto_times = Float64[]
    for _ in 1:runs
        push!(safe_times, run_once(guides, genome_path, motif; distance = 3, force_safe_minima = true))
        push!(auto_times, run_once(guides, genome_path, motif; distance = 3, force_safe_minima = false))
    end

    safe_median = median(safe_times)
    auto_median = median(auto_times)
    delta_pct = safe_median == 0.0 ? 0.0 : (auto_median - safe_median) / safe_median * 100

    println("Safe minima median (s): ", round(safe_median; digits = 6))
    println("Auto minima median (s): ", round(auto_median; digits = 6))
    println("Delta (auto-safe) %: ", round(delta_pct; digits = 2))

    if enforce && auto_median > safe_median
        error("No-regression gate failed: auto minima is slower than safe minima.")
    end

    try
        rm(TEST_DIR; recursive = true, force = true)
    catch
    end
end

main()
