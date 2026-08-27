#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Printf
using Statistics

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_GENOME = "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
const DEFAULT_CAS9_GUIDES = joinpath(
    ROOT_DIR, "test", "local_human", "data", "guides_for_tests.txt")
const DEFAULT_CAS12A_GUIDES = joinpath(
    ROOT_DIR, "test", "local_human", "data", "cas12a_guides_for_tests.txt")

parse_int_env(name, default) = parse(Int, strip(get(ENV, name, string(default))))
load_guides(path) = LongDNA{4}.(filter(!isempty, strip.(readlines(path))))

function alternating_times(f_avx2, f_avx512, runs::Int, warmups::Int)
    for _ in 1:warmups
        f_avx2()
        f_avx512()
    end
    times = Dict(:avx2 => Float64[], :avx512 => Float64[])
    values = Dict{Symbol, Any}()
    for run in 1:runs
        order = isodd(run) ? (:avx2, :avx512) : (:avx512, :avx2)
        for backend in order
            f = backend == :avx2 ? f_avx2 : f_avx512
            elapsed = @elapsed values[backend] = f()
            push!(times[backend], elapsed)
        end
    end
    return times, values
end

function scanner_case(raw, motif::Motif, label::String, runs::Int, warmups::Int)
    geometry = CHOPOFF.resolve_prefix_scan_geometry(motif, 3, 16)
    geometry === nothing && error("No optimized geometry for $label")
    kind = CHOPOFF.prefix_scan_kind(geometry)
    auto_eligible = kind in (:cas9, :cas12a)
    query = CHOPOFF.PrefixHashScanBitmaskQuery(Dict{UInt32, UInt64}(), 1)
    candidate_last = length(raw) - CHOPOFF.prefix_scan_candidate_last_offset(geometry)
    run(backend) = CHOPOFF.scan_prefix_hits_raw_range(
        geometry, raw, query, 1, candidate_last, 1, candidate_last,
        1, candidate_last, Val(backend))
    times, values = alternating_times(
        () -> run(:avx2), () -> run(:avx512), runs, warmups)
    values[:avx2][3] == values[:avx512][3] ||
        error("Scanner candidate parity failed for $label")
    avx2 = median(times[:avx2])
    avx512 = median(times[:avx512])
    return (
        stage = "scanner", workload = label, scan_kind = string(kind),
        auto_eligible, threads = Threads.nthreads(),
        avx2_median_s = avx2, avx512_median_s = avx512,
        speedup = avx2 / avx512, parity = true,
        passed = !auto_eligible || avx2 / avx512 >= 1.05,
    )
end

function end_to_end_case(
    genome, guides, motif::Motif, label::String, runs::Int, warmups::Int)

    geometry = CHOPOFF.resolve_prefix_scan_geometry(motif, 3, 16)
    geometry === nothing && error("No optimized geometry for $label")
    kind = CHOPOFF.prefix_scan_kind(geometry)
    auto_eligible = kind in (:cas9, :cas12a)
    tdir = mktempdir(prefix = "chopoff_avx512_")
    outputs = Dict(
        :avx2 => joinpath(tdir, "avx2.csv"),
        :avx512 => joinpath(tdir, "avx512.csv"),
    )
    run(backend) = begin
        search_prefixHashScan(
            guides, genome, outputs[backend]; motif, distance = 3,
            early_stopping = fill(1_000_000, 4),
            simd_backend = backend, output = :counts)
        return read(outputs[backend])
    end
    times, values = alternating_times(
        () -> run(:avx2), () -> run(:avx512), runs, warmups)
    parity = values[:avx2] == values[:avx512]
    avx2 = median(times[:avx2])
    avx512 = median(times[:avx512])
    return (
        stage = "end_to_end", workload = label, scan_kind = string(kind),
        auto_eligible,
        threads = Threads.nthreads(), avx2_median_s = avx2,
        avx512_median_s = avx512, speedup = avx2 / avx512,
        parity, passed = parity && (!auto_eligible || avx512 <= avx2 * 1.03),
    )
end

function main()
    CHOPOFF.can_use_prefix_hash_scan_avx512() ||
        error("AVX-512F/BW and BMI2 are required")
    genome = abspath(get(ENV, "CHOPOFF_AVX512_GENOME", DEFAULT_GENOME))
    cas9_path = abspath(get(ENV, "CHOPOFF_AVX512_CAS9_GUIDES", DEFAULT_CAS9_GUIDES))
    cas12a_path = abspath(get(
        ENV, "CHOPOFF_AVX512_CAS12A_GUIDES", DEFAULT_CAS12A_GUIDES))
    runs = parse_int_env("CHOPOFF_AVX512_RUNS", 11)
    warmups = parse_int_env("CHOPOFF_AVX512_WARMUPS", 2)
    raw_mb = parse_int_env("CHOPOFF_AVX512_RAW_MB", 32)
    out = abspath(get(
        ENV, "CHOPOFF_AVX512_OUT",
        "/tmp/prefix_hash_scan_avx512_threads$(Threads.nthreads()).csv"))
    isfile(genome) || error("Genome not found: $genome")
    isfile(genome * ".fai") || error("FASTA index not found: $(genome).fai")

    motifs = [
        ("Cas9", Motif("Cas9"; distance = 3), load_guides(cas9_path)),
        ("Cas12a", Motif("Cas12a"; distance = 3), load_guides(cas12a_path)),
        ("Cas9_NGA", Motif("Cas9_NGA"; distance = 3), load_guides(cas9_path)),
    ]
    seed = Vector{UInt8}(codeunits("ACGTACGTGCGTACGTTTACGGAT"))
    raw = repeat(seed, cld(raw_mb * 1024^2, length(seed)))
    resize!(raw, raw_mb * 1024^2)
    rows = NamedTuple[]
    for (label, motif, _) in motifs
        push!(rows, scanner_case(raw, motif, label, runs, warmups))
    end
    for (label, motif, guides) in motifs
        push!(rows, end_to_end_case(
            genome, guides, motif, label, runs, warmups))
    end
    frame = DataFrame(rows)
    mkpath(dirname(out))
    CSV.write(out, frame)
    for row in eachrow(frame)
        @printf("stage=%s workload=%s kind=%s auto=%s threads=%d avx2=%.6fs avx512=%.6fs speedup=%.3fx parity=%s gate=%s\n",
            row.stage, row.workload, row.scan_kind, row.auto_eligible,
            row.threads, row.avx2_median_s, row.avx512_median_s,
            row.speedup, row.parity, row.passed)
    end
    println("wrote $out")
    all(frame.passed) || error("AVX-512 qualification gate failed")
end

main()
