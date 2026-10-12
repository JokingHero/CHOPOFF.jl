using Random, Logging

vtest_site(name, g) = name == "Cas9" ? g * "AGG" : "TTTA" * g
vtest_random_dna(rng, n) = String(rand(rng, collect("ACGT"), n))

@testset "VCF coverage gaps" begin
    g, r = VCF_TEST_GUIDE, VCF_TEST_REFERENCE
    guide = LongDNA{4}(g)

    @testset "Multiple chromosomes; empty chromosome between; edge fragment packed after another" begin
        edge = vtest_rc(g * "AGG") * "T"^60
        chromosomes = [
            "chr1" => (r, [vtest_edit(r, 45, 46, "A", "c1snp"), vtest_edit(r, 51, 51, "C", "c1ins")]),
            "chr2" => ("C"^30 * g * "AGG" * "C"^30, []),
            "chr3" => (edge, [vtest_edit(edge, 5, 6, edge[5] == 'A' ? "C" : "A", "c3edge"),
                vtest_edit(edge, 10, 12, "", "c3del")]),
        ]
        for chunk in (64, 2 * 1024 * 1024), threads in unique([1, min(4, Threads.nthreads())])
            mktempdir() do dir
                rows = vtest_compare_multi(dir, chromosomes, [guide], Motif("Cas9"; distance = 2), 2;
                    chunk, threads, case_id = "multi_chrom_$chunk")
                @test Set(rows.chromosome) == Set(["chr1", "chr2", "chr3"])
                @test any(v -> !ismissing(v) && occursin("c1snp", v), rows.variants)
                @test any(v -> !ismissing(v) && occursin("c3edge", v), rows.variants)
            end
        end
        mktempdir() do dir
            only_last = ["chr1" => (r, []), "chr2" => (r, [vtest_edit(r, 45, 46, "A", "late")])]
            rows = vtest_compare_multi(dir, only_last, [guide], Motif("Cas9"; distance = 1), 1; case_id = "variants_on_chr2")
            @test all(row -> ismissing(row.variants) || row.chromosome == "chr2", eachrow(rows))
            @test any(!ismissing, rows.variants)
        end
    end

    @testset "Several guides on a random background, $name" for name in ("Cas9", "Cas12a")
        rng = MersenneTwister(name == "Cas9" ? 11 : 12)
        n = name == "Cas9" ? 20 : 21
        g1, g2 = vtest_random_dna(rng, n), vtest_random_dna(rng, n)
        left = vtest_random_dna(rng, 40)
        reference = left * vtest_site(name, g1) * vtest_random_dna(rng, 30) *
            vtest_rc(vtest_site(name, g2)) * vtest_random_dna(rng, 40)
        offset = name == "Cas9" ? 40 : 44                     # first guide base of site 1
        p1 = offset + 5
        alt1 = vtest_alternative(reference[p1])
        site2 = 40 + length(vtest_site(name, g1)) + 30          # base before site 2
        middle2 = site2 + length(vtest_site(name, g2)) ÷ 2
        q = offset + n + 8
        edits = [vtest_edit(reference, p1, p1 + 1, alt1, "s1"),
            vtest_edit(reference, middle2, middle2, "GA", "ins2"),
            vtest_edit(reference, offset + 12, offset + 14,
                join(vtest_alternative(c) for c in reference[(offset + 12):(offset + 13)]), "mnp1"),
            vtest_edit(reference, q, q + 2, "", "flankdel")]
        variant_guide = collect(g1); variant_guide[5] = only(alt1)
        iupac_guide = collect(g2); iupac_guide[9] = 'N'
        guides = LongDNA{4}.([g1, String(variant_guide), g2, String(iupac_guide), vtest_random_dna(rng, n)])
        mktempdir() do dir
            rows = vtest_compare(dir, reference, edits, guides, Motif(name; distance = 2), 2;
                chunk = 97, threads = min(4, Threads.nthreads()), case_id = "multiguide_$name")
            @test length(unique(rows.guide[.!ismissing.(rows.variants)])) >= 3
            @test any(v -> !ismissing(v) && occursin("ins2", v), rows.variants)
        end
    end

    @testset "Reference ambiguity budget $budget with packed fragments" for budget in 1:3
        block(i) = begin
            b = collect(r)
            b[42] = 'N'                                        # inside the guide
            b[45 - 26] = 'N'                                   # first fragment base (W = 26)
            b[45 + 26] = 'N'                                   # last fragment base
            budget >= 2 && (b[50] = 'R')
            String(b)
        end
        reference = join(block(i) for i in 1:3)
        n = length(r)
        edits = [vtest_edit(reference, (i - 1) * n + 45, (i - 1) * n + 46, "A", "b$i") for i in 1:3]
        mktempdir() do dir
            rows = vtest_compare(dir, reference, edits, [guide],
                Motif("Cas9"; distance = 3, ambig_max = budget), 3; case_id = "ambig_packed_$budget")
            @test count(v -> !ismissing(v), rows.variants) > 0
        end
    end

    @testset "Trimming paths" begin
        # GCG>GG: shared prefix and suffix. Deleting C creates the PAM, so the
        # trimmed suffix base is a site endpoint and must keep its coordinate.
        t = "T"^40 * g * "AGCG" * "T"^60
        mktempdir() do dir
            edits = [vtest_padded(t, 62, 3, "GG", 63, 64, "", "suffix")]
            @test t[62:64] == "GCG"
            rows = vtest_compare(dir, t, edits, [guide], Motif("Cas9"; distance = 2), 2; case_id = "suffix_trim")
            @test any(row -> row.distance == 0 && !ismissing(row.variants) &&
                occursin("suffix|chr1:62:GCG>GG", row.variants) &&
                ismissing(row.start_origin) && ismissing(row.end_origin), eachrow(rows))
        end
        # Multiallelic STR record ATT>A,AT,ATTT; the reader trims each ALT on its own.
        s = "T"^40 * "ACGTATTTGCAGTACGATCG" * "AGG" * "T"^60
        mktempdir() do dir
            @test s[45:47] == "ATT"
            edits = [vtest_padded(s, 45, 3, "A", 46, 48, "", "str"; group = "str"),
                vtest_padded(s, 45, 3, "AT", 47, 48, "", "str"; group = "str"),
                vtest_padded(s, 45, 3, "ATTT", 48, 48, "T", "str"; group = "str")]
            rows = vtest_compare(dir, s, edits, [LongDNA{4}(s[41:60])], Motif("Cas9"; distance = 2), 2; case_id = "str_multiallelic")
            for alt in ("A", "AT", "ATTT")
                @test any(v -> !ismissing(v) && occursin("str|chr1:45:ATT>$alt", v), rows.variants)
            end
            @test all(v -> ismissing(v) || !occursin(';', v), rows.variants)
        end
        # MNP keeping an interior reference base.
        mktempdir() do dir
            alt = vtest_alternative(r[46]) * r[47] * vtest_alternative(r[48])
            rows = vtest_compare(dir, r, [vtest_edit(r, 46, 49, alt, "interior")], [guide],
                Motif("Cas9"; distance = 3), 3; case_id = "mnp_interior")
            @test any(v -> !ismissing(v) && occursin("interior", v), rows.variants)
        end
        # IUPAC ALT inside an insertion: realizations share one origin axis.
        mktempdir() do dir
            edits = [vtest_edit(r, 52, 52, "GA", "iupac_ins"; original = "RA", group = "iupac_ins"),
                vtest_edit(r, 52, 52, "AA", "iupac_ins"; original = "RA", group = "iupac_ins")]
            rows = vtest_compare(dir, r, edits, [guide], Motif("Cas9"; distance = 2), 2; case_id = "iupac_insertion")
            for concrete in ("GA", "AA")
                @test any(v -> !ismissing(v) && occursin("resolved=" * r[51] * concrete, v), rows.variants)
            end
        end
    end

    @testset "Insertion gaps at reference-site boundaries" begin
        # The d0 reference site covers 41:63 (guide plus PAM).
        edits = [vtest_edit(r, gap, gap, "C", "gap$gap") for gap in (41, 42, 63, 64)]
        mktempdir() do dir
            rows = vtest_compare(dir, r, edits, [guide], Motif("Cas9"; distance = 0), 0; case_id = "gap_boundaries_overlap")
            reference_row = only(eachrow(rows[ismissing.(rows.variants), :]))
            ids = Set(first(split(t, '|')) for t in vtest_tokens(reference_row.overlapping_variants))
            @test ids == Set(["gap42", "gap63"])
        end
    end

    @testset "Chromosome edges: $name $(strand) $(side)" for name in ("Cas9", "Cas12a"),
            strand in ("+", "-"), side in (:start, :end)
        n = name == "Cas9" ? 20 : 21
        rng = MersenneTwister(hash((name, strand, side)))
        gg = vtest_random_dna(rng, n)
        site = vtest_site(name, gg)
        site = strand == "+" ? site : vtest_rc(site)
        filler = "C"^50
        reference = side == :start ? site * filler : filler * site
        base = side == :start ? 0 : length(filler)
        p = base + length(site) ÷ 2
        edits = [vtest_edit(reference, p, p + 1, vtest_alternative(reference[p]), "edge"),
            vtest_edit(reference, base + 2, base + 3, "", "edge_del")]
        for chunk in (64, 2 * 1024 * 1024)
            mktempdir() do dir
                rows = vtest_compare(dir, reference, edits, [LongDNA{4}(gg)], Motif(name; distance = 3), 3;
                    chunk, case_id = "edge_$(name)_$(strand)_$(side)")
                @test any(!ismissing, rows.variants)
            end
        end
    end

    @testset "Soft-masked reference and lowercase VCF alleles" begin
        mktempdir() do dir
            edits = [vtest_edit(r, 45, 46, "A", "low_snp"), vtest_edit(r, 51, 51, "C", "low_ins")]
            masked = lowercase(r[1:50]) * r[51:end]
            genome = vtest_fasta(joinpath(dir, "masked.fa"), ["chr1" => masked])
            records = [(c, p, id, lowercase(ref), lowercase(alt), q, f, i)
                for (c, p, id, ref, alt, q, f, i) in vtest_records(edits)]
            vcf = vtest_vcf(joinpath(dir, "lower.vcf"), records)
            rows = vtest_compare(dir, r, edits, [guide], Motif("Cas9"; distance = 2), 2;
                genome, vcf, case_id = "soft_masked")
            @test any(v -> !ismissing(v) && occursin("low_snp", v), rows.variants)
        end
    end

    @testset "Variant early stopping is per guide" begin
        g1 = g
        g2 = "TGCA" * g[5:end]
        g3 = "GACT" * g[5:end]
        block(x) = "T"^40 * x * "AGG" * "T"^40
        reference = block(g1) * block(g1) * block(g2) * block(g3)
        n = length(block(g1))
        edits = [vtest_edit(reference, (i - 1) * n + 50, (i - 1) * n + 51, "A", "v$i") for i in 1:4]
        guides = LongDNA{4}.([begin x = collect(s); x[10] = 'A'; String(x) end for s in (g1, g2, g3)])
        mktempdir() do dir
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => reference])
            vcf = vtest_vcf(joinpath(dir, "input.vcf"), vtest_records(edits))
            full, capped = joinpath(dir, "full.csv"), joinpath(dir, "capped.csv")
            CHOPOFF.search_prefixHashScan(guides, genome, Motif("Cas9"; distance = 1), full;
                distance = 1, vcf, early_stopping = fill(typemax(Int), 2), stream_chunk_bases = 64)
            stats = PrefixHashScanStats()
            @test_logs (:warn, r"variant early-stopping") match_mode = :any CHOPOFF.search_prefixHashScan(
                guides, genome, Motif("Cas9"; distance = 1), capped;
                distance = 1, vcf, early_stopping = [1, 1000], stream_chunk_bases = 64, stats)
            @test stats.variant_retired_guides == 1
            f, c = DataFrame(CSV.File(full)), DataFrame(CSV.File(capped))
            @test issubset(Set(vtest_signature.(eachrow(c))), Set(vtest_signature.(eachrow(f))))
            variant_rows(df, s) = Set(vtest_signature(row) for row in eachrow(df)
                if row.guide == string(s) && !ismissing(row.variants))
            @test count(row -> !ismissing(row.variants) && row.distance == 0, eachrow(f[f.guide .== string(guides[1]), :])) == 2
            @test count(row -> !ismissing(row.variants) && row.distance == 0, eachrow(c[c.guide .== string(guides[1]), :])) == 1
            for s in guides[2:3]
                @test variant_rows(c, s) == variant_rows(f, s)
                @test !isempty(variant_rows(c, s))
            end
        end
    end

    @testset "Dense combination warning" begin
        mktempdir() do dir
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            out = joinpath(dir, "out.csv")
            for (count, warns, expected) in ((13, false, 4096), (14, true, 8192))
                edits = [vtest_edit(r, p, p + 1, vtest_alternative(r[p]), "d$p") for p in 41:(40 + count)]
                vcf = vtest_vcf(joinpath(dir, "dense$count.vcf"), vtest_records(edits))
                stats = PrefixHashScanStats()
                logs, _ = Test.collect_test_logs() do
                    CHOPOFF.search_prefixHashScan([guide], genome, Motif("Cas9"; distance = 0), out;
                        distance = 0, vcf, stats, vcf_max_alleles = count)
                end
                messages = [string(l.message) for l in logs]
                @test any(m -> occursin("exceeds the default", m), messages)
                @test any(m -> occursin("combinations are dense", m), messages) == warns
                @test stats.max_anchor_fragments == expected
                @test stats.vcf_truncated_anchors == 0
            end
        end
    end

    @testset "ALT allele cap" begin
        # Ten SNPs inside one site: the cap keeps haplotypes with at most k ALTs.
        dense = [vtest_edit(r, p, p + 1, vtest_alternative(r[p]), "c$p") for p in 42:2:60]
        g10 = collect(g)
        for v in dense
            g10[v.first - 40] = only(v.alt)
        end
        guides = LongDNA{4}.([g, String(g10)])
        @testset "SNP cluster k=$k" for k in 1:4
            mktempdir() do dir
                stats = PrefixHashScanStats()
                rows = vtest_compare(dir, r, dense, guides, Motif("Cas9"; distance = 3), 3;
                    max_alleles = k, case_id = "cap_snp_$k")
                @test maximum(row -> length(vtest_tokens(row.variants)), eachrow(rows)) <= k
                # At d=3 the reference guide reaches at most three ALT SNPs.
                @test any(row -> length(vtest_tokens(row.variants)) == min(k, 3), eachrow(rows))
            end
        end
        @testset "Mixed cluster k=$k" for k in 2:3
            mixed = [vtest_edit(r, 44, 45, "A", "snp"), vtest_edit(r, 47, 47, "C", "ins"),
                vtest_edit(r, 50, 52, "", "del"),
                vtest_edit(r, 54, 55, "A", "amb"; original = "R", group = "amb"),
                vtest_edit(r, 54, 55, "G", "amb"; original = "R", group = "amb"),
                vtest_edit(r, 57, 58, "T", "snp2")]
            mktempdir() do dir
                vtest_compare(dir, r, mixed, [guide], Motif("Cas9"; distance = 3), 3;
                    max_alleles = k, case_id = "cap_mixed_$k")
            end
        end
        @testset "Clusters farther apart than 2W stay complete at k=1" begin
            far = r * r * r
            n = length(r)
            edits = [vtest_edit(far, (i - 1) * n + 45, (i - 1) * n + 46, "A", "far$i") for i in 1:3]
            mktempdir() do dir
                rows = vtest_compare(dir, far, edits, [guide], Motif("Cas9"; distance = 2), 2;
                    max_alleles = 1, case_id = "cap_far")
                @test Set(first(split(only(vtest_tokens(v)), '|')) for v in skipmissing(rows.variants)) ==
                    Set(["far1", "far2", "far3"])
            end
        end
        @testset "Warnings, stats and validation" begin
            mktempdir() do dir
                genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
                vcf = vtest_vcf(joinpath(dir, "dense.vcf"), vtest_records(dense))
                out = joinpath(dir, "out.csv")
                search(; kw...) = CHOPOFF.search_prefixHashScan(guides, genome, Motif("Cas9"; distance = 3), out;
                    distance = 3, vcf, kw...)
                stats = PrefixHashScanStats()
                @test_logs (:warn, r"truncated at vcf_max_alleles=2") match_mode = :any search(; vcf_max_alleles = 2, stats)
                @test stats.vcf_truncated_anchors > 0
                stats = PrefixHashScanStats()
                logs, _ = Test.collect_test_logs() do
                    search(; vcf_max_alleles = 10, stats)
                end
                @test stats.vcf_truncated_anchors == 0
                @test !any(l -> occursin("truncated", string(l.message)), logs)
                @test any(l -> occursin("exceeds the default", string(l.message)), logs)
                @test_throws ErrorException search(; vcf_max_alleles = 0)
                @test_throws ErrorException search_prefixHashScan([guide], genome, out; distance = 1, vcf_max_alleles = 2)
                @test_throws ErrorException CHOPOFF.search_prefixHashScan([guide], genome, Motif("Cas9"; distance = 1), out;
                    distance = 1, vcf_max_alleles = 2)
            end
        end
    end
    @testset "IUPAC expansion limit" begin
        mktempdir() do dir
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            out = joinpath(dir, "out.csv")
            for (alt, skipped) in (("NNN", false), ("NNNN", true))
                vcf = vtest_vcf(joinpath(dir, "iupac$alt.vcf"), [("chr1", 45, "amb", "GTCA", alt, ".", ".", ".")])
                stats = PrefixHashScanStats()
                logs, _ = Test.collect_test_logs() do
                    CHOPOFF.search_prefixHashScan([guide], genome, Motif("Cas9"; distance = 2), out;
                        distance = 2, vcf, stats)
                end
                @test any(l -> occursin("IUPAC realizations ignored", string(l.message)), logs) == skipped
                rows = DataFrame(CSV.File(out))
                @test any(!ismissing, rows.variants) == !skipped
            end
        end
    end
end
