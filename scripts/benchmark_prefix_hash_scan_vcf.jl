# julia --project=. -t 8 scripts/benchmark_prefix_hash_scan_vcf.jl GENOME VCF GUIDES CHROM|all OUTDIR [RUNS]
# Compare uncapped Cas9 distance-3 searches. VCF preparation is timed on every run.
using CHOPOFF, BioSequences, CSV, DataFrames, Statistics

function main(args)
    length(args) in (5, 6) || error("Expected GENOME VCF GUIDES CHROM OUTDIR [RUNS]")
    genome, vcf, guide_file, chrom, output = args[1:5]
    runs = length(args) == 6 ? parse(Int, args[6]) : 5
    runs > 0 || error("RUNS must be positive")
    guides = LongDNA{4}.(filter(!isempty, strip.(readlines(guide_file))))
    length(guides) <= 64 || error("This per-batch benchmark accepts at most 64 guides")
    mkpath(output)
    motif = Motif("Cas9"; distance = 3)
    expected = Dict{Bool, Vector{UInt8}}()
    measurements = NamedTuple[]
    function scan(with_vcf, repetition)
        stats = PrefixHashScanStats()
        path = joinpath(output, with_vcf ? "variants.csv" : "reference.csv")
        GC.gc()
        timing = @timed CHOPOFF.search_prefixHashScan(guides, genome, motif, path;
            distance = 3, vcf = with_vcf ? vcf : nothing, region = chrom == "all" ? nothing : chrom,
            early_stopping = fill(typemax(Int), 4), stats = stats)
        bytes = read(path)
        if repetition == 0
            expected[with_vcf] = bytes
        else
            bytes == expected[with_vcf] || error("Output differs between repetitions")
            push!(measurements, (; repetition, with_vcf, seconds = timing.time,
                allocated_bytes = timing.bytes, peak_process_rss = Sys.maxrss(),
                preparation_s = stats.enumeration_ns / 1e9,
                worker_io_s = stats.variant_io_ns / 1e9,
                stats.variants_read, stats.fragments, stats.fragment_bases,
                stats.descriptor_bytes, stats.variant_candidates, stats.anchor_rejected,
                stats.variant_rows))
        end
        println("repetition=$repetition vcf=$with_vcf seconds=$(round(timing.time; digits=3))")
        flush(stdout)
    end
    scan(false, 0)
    scan(true, 0)
    reference = DataFrame(CSV.File(joinpath(output, "reference.csv")))
    variants = DataFrame(CSV.File(joinpath(output, "variants.csv")))
    select(variants[ismissing.(variants.variants), :], names(reference)) == reference ||
        error("Reference rows differ with VCF enabled")
    for repetition in 1:runs
        for with_vcf in (isodd(repetition) ? (false, true) : (true, false))
            scan(with_vcf, repetition)
        end
    end
    CSV.write(joinpath(output, "timings.csv"), DataFrame(measurements))
    baseline = median(m.seconds for m in measurements if !m.with_vcf)
    variant = median(m.seconds for m in measurements if m.with_vcf)
    println("median reference=$(round(baseline; digits=3))s vcf=$(round(variant; digits=3))s ratio=$(round(variant / baseline; digits=3))")
    println("Reference rows match; repeated outputs are byte-identical. Peak RSS is process-wide, including warmup.")
end

main(ARGS)
