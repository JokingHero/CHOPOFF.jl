using Random

# BED rows are (chrom, start, end), 0-based half-open.
function vtest_bed(path, rows)
    open(path, "w") do io
        println(io, "# comment")
        println(io, "track name=targets")
        foreach(r -> println(io, join(r, '\t'), "\textra\t0\t+"), rows)
    end
    return path
end

vtest_bed_overlaps(rows, chrom, lo, hi) = any(((c, a, b),) -> c == chrom && a + 1 <= hi && lo <= b, rows)

# Reference rows use the independent span; variant rows use their lifted endpoints.
function vtest_row_span(row, motif)
    hasproperty(row, :variants) && !ismissing(row.variants) && return minmax(row.start, row.end)
    return vtest_span(row, motif)
end

function vtest_region_rows(path)
    df = DataFrame(CSV.File(path; types = String))
    return sort!([Tuple(coalesce(x, "") for x in r) for r in eachrow(df)])
end

function vtest_expected_region_rows(path, bed, motif)
    df = DataFrame(CSV.File(path))
    keep = [vtest_bed_overlaps(bed, string(r.chromosome), vtest_row_span(r, motif)...) for r in eachrow(df)]
    filtered = joinpath(dirname(path), "expected_" * basename(path))
    CSV.write(filtered, df[keep, :])
    return vtest_region_rows(filtered)
end

function vtest_random_region_case(rng, name, n_guides)
    n = name == "Cas9" ? 20 : 21
    guides = [vtest_random_dna(rng, n) for _ in 1:n_guides]
    chromosomes = Pair{String, Tuple{String, Vector}}[]
    for chrom in ("chr1", "chr2")
        parts = [vtest_random_dna(rng, 30)]
        for _ in 1:6
            g = collect(rand(rng, guides))
            for _ in 1:rand(rng, 0:2)
                p = rand(rng, eachindex(g))
                g[p] = only(vtest_alternative(g[p]))
            end
            site = vtest_site(name, String(g))
            push!(parts, rand(rng, Bool) ? site : vtest_rc(site), vtest_random_dna(rng, rand(rng, 5:40)))
        end
        reference = join(parts)
        edits = []
        positions = sort(rand(rng, 2:(length(reference) - 70), 8))
        for (i, p) in enumerate(positions)
            kind = rand(rng, 1:4)
            id = "$(chrom)_e$i"
            push!(edits, kind == 1 ? vtest_edit(reference, p, p + 1, vtest_alternative(reference[p]), id) :
                kind == 2 ? vtest_edit(reference, p, p, vtest_random_dna(rng, rand(rng, 1:4)), id) :
                kind == 3 ? vtest_edit(reference, p, p + rand(rng, 1:3), "", id) :
                vtest_edit(reference, p, p + rand(rng, 30:60), "", id))
        end
        push!(chromosomes, chrom => (reference, edits))
    end
    bed = Tuple{String, Int, Int}[]
    for (chrom, (reference, _)) in chromosomes, _ in 1:rand(rng, 1:3)
        a = rand(rng, 0:(length(reference) - 2))
        push!(bed, (chrom, a, min(length(reference), a + rand(rng, [1, 5, 40, 120]))))
    end
    return LongDNA{4}.(guides), chromosomes, bed
end

@testset "BED regions [D03]" begin
    @testset "Matches a whole-chromosome search filtered by site overlap, $name" for name in ("Cas9", "Cas12a")
        rng = MersenneTwister(name == "Cas9" ? 31 : 32)
        for case in 1:4, distance in 0:4
            guides, chromosomes, bed = vtest_random_region_case(rng, name, 3)
            motif = Motif(name; distance)
            mktempdir() do dir
                genome = vtest_fasta(joinpath(dir, "reference.fa"), [c => r for (c, (r, _)) in chromosomes])
                records = reduce(vcat, (vtest_records(e; chrom = c) for (c, (_, e)) in chromosomes); init = Tuple[])
                vcf = vtest_vcf(joinpath(dir, "input.vcf"), sort(records; by = r -> (r[1], r[2])))
                bed_path = vtest_bed(joinpath(dir, "targets.bed"), bed)
                limits = fill(typemax(Int), distance + 1)
                full_ref, full_vcf = joinpath(dir, "full_ref.csv"), joinpath(dir, "full_vcf.csv")
                search_prefixHashScan(guides, genome, motif, full_ref; distance, early_stopping = limits)
                search_prefixHashScan(guides, genome, motif, full_vcf; distance, vcf, early_stopping = limits)
                expected_ref = vtest_expected_region_rows(full_ref, bed, motif)
                expected_vcf = vtest_expected_region_rows(full_vcf, bed, motif)
                for (chunk, threads) in ((64, 1), (2 * 1024 * 1024, min(4, Threads.nthreads())))
                    out_ref, out_vcf = joinpath(dir, "bed_ref.csv"), joinpath(dir, "bed_vcf.csv")
                    search_prefixHashScan(guides, genome, motif, out_ref; distance, region = bed_path,
                        early_stopping = limits, stream_chunk_bases = chunk, scan_threads = threads)
                    search_prefixHashScan(guides, genome, motif, out_vcf; distance, vcf, region = bed_path,
                        early_stopping = limits, stream_chunk_bases = chunk, scan_threads = threads)
                    @test vtest_region_rows(out_ref) == expected_ref
                    @test vtest_region_rows(out_vcf) == expected_vcf
                end
            end
        end
    end

    @testset "A deletion brings a distant anchor into the interval" begin
        rng = MersenneTwister(33)
        g = VCF_TEST_GUIDE
        filler = vtest_random_dna(rng, 100)
        reference = "T"^40 * g[1:10] * filler * g[11:20] * "AGG" * "T"^40
        snp = vtest_edit(reference, 43, 44, vtest_alternative(reference[43]), "far_snp")
        deletion = vtest_edit(reference, 51, 151, "", "bridge")
        variant_guide = collect(g); variant_guide[3] = only(vtest_alternative(g[3]))
        mktempdir() do dir
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => reference])
            vcf = vtest_vcf(joinpath(dir, "input.vcf"), vtest_records([snp, deletion]))
            # The interval lies past the PAM-side bases, far right of SNP + W.
            bed = vtest_bed(joinpath(dir, "targets.bed"), [("chr1", 155, 160)])
            output = joinpath(dir, "out.csv")
            search_prefixHashScan([LongDNA{4}(String(variant_guide))], genome, Motif("Cas9"; distance = 0), output;
                distance = 0, vcf, region = bed)
            rows = DataFrame(CSV.File(output))
            @test nrow(rows) == 1
            @test occursin("far_snp", rows.variants[1]) && occursin("bridge", rows.variants[1])
        end
    end

    @testset "Sites outside the region spend no early-stopping budget" begin
        g = VCF_TEST_GUIDE
        site = g * "AGG"
        reference = "T"^40 * site * site * "T"^40
        mktempdir() do dir
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => reference])
            bed = vtest_bed(joinpath(dir, "targets.bed"), [("chr1", 63, 64)])   # first base of site 2
            full, inside = joinpath(dir, "full.csv"), joinpath(dir, "inside.csv")
            search_prefixHashScan([LongDNA{4}(g)], genome, Motif("Cas9"; distance = 2), full;
                distance = 2, early_stopping = [1, 1, 1])
            search_prefixHashScan([LongDNA{4}(g)], genome, Motif("Cas9"; distance = 2), inside;
                distance = 2, early_stopping = [1, 1, 1], region = bed)
            full_rows, inside_rows = DataFrame(CSV.File(full)), DataFrame(CSV.File(inside))
            @test full_rows.start == [63]
            @test inside_rows.start == [86]
        end
    end

    @testset "Enumeration blocks match across thread counts [D07]" begin
        rng = MersenneTwister(34)
        reference = vtest_random_dna(rng, 600)
        edits = [vtest_edit(reference, p, p + 1, vtest_alternative(reference[p]), "s$p") for p in 20:7:560]
        append!(edits, [vtest_edit(reference, p, p + 4, "", "d$p") for p in 25:31:540])
        sort!(edits; by = v -> v.pos)
        mktempdir() do dir
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => reference])
            vcf = vtest_vcf(joinpath(dir, "input.vcf"), vtest_records(edits))
            dbi, lengths = CHOPOFF.prefix_hash_scan_dbinfo(genome, Motif("Cas9"; distance = 2))
            data = CHOPOFF.prepare_prefix_hash_scan_variants(vcf, dbi, lengths, trues(1), String[])
            blocks = CHOPOFF.enumerate_variant_chromosome(data, dbi, lengths, 1, 4)
            @test length(blocks) == 16
            # Materialized fragments in block order equal one serial block's.
            sequences(bs) = [String(CHOPOFF.materialize_variant_fragment!(UInt8[], data, b, f, 1, f.bases))
                for b in bs for f in b.fragments]
            single = CHOPOFF.enumerate_variant_block(data, dbi, lengths, data.anchors[1], 1,
                length(dbi.motif) + 2, 0, nothing)
            @test length(single.fragments) > 500
            @test sequences(blocks) == sequences([single])
            @test reduce(vcat, (b.per_anchor for b in blocks)) == single.per_anchor
            outputs = map((1, 4)) do threads
                path = joinpath(dir, "out$threads.csv")
                search_prefixHashScan([LongDNA{4}(reference[100:119])], genome, Motif("Cas9"; distance = 2), path;
                    distance = 2, vcf, scan_threads = threads)
                read(path, String)
            end
            @test outputs[1] == outputs[2]
        end
    end
end
