#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Random
using Statistics

const TEST_DIR = mktempdir()
const PARITY_COLS = [
    :guide, :alignment_guide, :alignment_reference, :distance,
    :chromosome, :start, :strand,
]

function build_fixture(n_bases::Int, n_guides::Int)
    rng = MersenneTwister(0x5a55)
    genome_seq = String(rand(rng, ['A', 'C', 'G', 'T'], n_bases))
    guides = LongDNA{4}.([
        String(rand(rng, ['A', 'C', 'G', 'T'], 20)) for _ in 1:n_guides
    ])
    genome_path = joinpath(TEST_DIR, "bench_genome.fa")
    open(genome_path, "w") do io
        write(io, ">chr1\n", genome_seq, "\n")
    end
    open(genome_path * ".fai", "w") do io
        write(io, "chr1\t$n_bases\t6\t$n_bases\t$(n_bases + 1)\n")
    end
    return genome_path, guides
end

function run_once(guides, genome_path, motif, backend::Symbol, run_id::Int)
    output_path = joinpath(TEST_DIR, "$(backend)_$(run_id).csv")
    elapsed = @elapsed search_sassy(
        guides, genome_path, motif, output_path;
        distance = motif.distance,
        backend = backend,
    )
    frame = CSV.read(output_path, DataFrame)
    return elapsed, sort(select(frame, PARITY_COLS), PARITY_COLS)
end

function supported_backends()
    backends = Symbol[:auto, :avx2_safe]
    CHOPOFF.Sassy.can_use_bmi2_pext() && push!(backends, :avx2_pext)
    CHOPOFF.Sassy.can_use_avx512() &&
        CHOPOFF.Sassy.can_use_bmi2_pext() && push!(backends, :avx512)
    return backends
end

function main()
    runs = parse(Int, get(ENV, "SASSY_BENCH_RUNS", "7"))
    n_bases = parse(Int, get(ENV, "SASSY_BENCH_BASES", "8000000"))
    n_guides = parse(Int, get(ENV, "SASSY_BENCH_GUIDES", "8"))
    tolerance = parse(Float64, get(ENV, "SASSY_BENCH_TOLERANCE", "0.03"))
    enforce = get(ENV, "SASSY_BENCH_ENFORCE_NO_REGRESSION", "0") == "1"
    output_path = get(ENV, "SASSY_BENCH_OUT", "")

    genome_path, guides = build_fixture(n_bases, n_guides)
    motif = Motif("Cas9"; distance = 3)
    backends = supported_backends()
    resolved_auto = CHOPOFF.Sassy.resolve_sassy_backend()

    println("SASSY backend benchmark")
    println("cpu: ", Sys.CPU_NAME)
    println("threads: ", Threads.nthreads())
    println("avx2: ", CHOPOFF.Sassy.can_use_avx2())
    println("bmi2: ", CHOPOFF.Sassy.can_use_bmi2_pext())
    println("avx512f_bw: ", CHOPOFF.Sassy.can_use_avx512())
    println("auto backend: ", resolved_auto)
    println("fixture: $(n_bases) bases, $(n_guides) guides")

    for backend in backends
        run_once(guides, genome_path, motif, backend, 0)
    end

    rows = NamedTuple[]
    expected = nothing
    timings = Dict(backend => Float64[] for backend in backends)
    for run_id in 1:runs
        order = isodd(run_id) ? backends : reverse(backends)
        for backend in order
            elapsed, frame = run_once(guides, genome_path, motif, backend, run_id)
            expected === nothing && (expected = frame)
            parity = frame == expected
            parity || error("Output parity failed for backend=$backend run=$run_id")
            push!(timings[backend], elapsed)
            push!(rows, (
                cpu = Sys.CPU_NAME,
                threads = Threads.nthreads(),
                requested_backend = backend,
                resolved_backend = CHOPOFF.Sassy.resolve_sassy_backend(backend),
                run = run_id,
                elapsed_s = elapsed,
                bases = n_bases,
                guides = n_guides,
                parity = parity,
            ))
        end
    end

    for backend in backends
        println(backend, " median_s=", round(median(timings[backend]); digits = 6))
    end

    reference = :avx2_pext
    if reference in backends && resolved_auto in backends
        ratio = median(timings[resolved_auto]) / median(timings[reference])
        println("auto/reference ratio: ", round(ratio; digits = 4))
        if enforce && ratio > 1 + tolerance
            error("Backend regression: ratio=$ratio exceeds $(1 + tolerance)")
        end
    end

    isempty(output_path) || CSV.write(output_path, DataFrame(rows))
    rm(TEST_DIR; recursive = true, force = true)
end

main()
