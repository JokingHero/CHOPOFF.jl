@testset "VCF regressions" begin
    mktempdir() do dir
        guide = dna"ACGTGTCAGTACGATCGTAC"
        reference = repeat("T", 40) * string(guide) * "AGG" * repeat("T", 80)

        @testset "reference cap retirement and variant deduplication" begin
            genome = vtest_fasta(joinpath(dir, "retire.fa"), ["chr1" => reference * reference])
            vcf = vtest_vcf(joinpath(dir, "retire.vcf"),
                [vtest_record(vtest_edit(reference, 45, 46, "A", "variant"))])
            output = joinpath(dir, "retire.csv")
            stats = PrefixHashScanStats()
            CHOPOFF.search_prefixHashScan([guide], genome, Motif("Cas9"; distance = 1), output;
                distance = 1, vcf = vcf, early_stopping = [1, 1], stream_chunk_bases = 64, stats = stats)
            @test stats.variant_rows == 1
            @test count(!ismissing, DataFrame(CSV.File(output)).variants) == 1
            # Right-hand choices outside a hit produce duplicate fragments,
            # but must not spend the variant cap or retire the guide.
            vcf = vtest_vcf(joinpath(dir, "duplicate_fragments.vcf"),
                [vtest_record(vtest_edit(reference, 45, 46, "A", "site")),
                    vtest_record(vtest_edit(reference, 68, 69, "C", "flank"))])
            CHOPOFF.search_prefixHashScan([guide], genome, Motif("Cas9"; distance = 1), output;
                distance = 1, vcf = vcf, early_stopping = [100, 1], stream_chunk_bases = 64, stats = stats)
            @test stats.variant_rows == 1
            @test stats.variant_retired_guides == 0
            vcf = vtest_vcf(joinpath(dir, "variant_retire.vcf"),
                [vtest_record(vtest_edit(reference * reference, p, p + 1, "A", "site$p"))
                    for p in (45, 45 + length(reference))])
            @test_logs (:warn, r"VCF search reached variant early-stopping limits") CHOPOFF.search_prefixHashScan(
                [guide], genome, Motif("Cas9"; distance = 1), output;
                distance = 1, vcf = vcf, early_stopping = [100, 1], stream_chunk_bases = 64, stats = stats)
            @test stats.variant_rows == 1
            @test stats.variant_retired_guides == 1
        end

        @testset "example-genome fixture and CLI" begin
            root = joinpath(@__DIR__, "..", "sample_data")
            vcf = joinpath(root, "vcf", "small_norm.vcf")
            genome = joinpath(root, "genome", "semirandom.fa")
            fixture_guide = LongDNA{4}(strip(read(joinpath(root, "vcf", "guide.txt"), String)))
            output = joinpath(dir, "fixture.csv")
            search_prefixHashScan([fixture_guide], genome, output; distance = 2, vcf = vcf, region = "semirandom1")
            rows = DataFrame(CSV.File(output))
            @test any(!ismissing, rows.variants)
            args = CHOPOFF.parse_commandline(["search", "--guides", "guides.txt", "--output", output,
                "prefixHashScan", "--genome", genome, "--vcf", vcf, "--region", "semirandom1", "--vcf_annotate", "INFO/AF"])
            @test args["search"]["prefixHashScan"]["vcf"] == vcf
            summary = joinpath(dir, "summary.csv")
            CHOPOFF.main(["summarize", "--detail_file", output, "--output", summary, "--variant_aware", "--distance", "2"])
            @test DataFrame(CSV.File(summary)) == summarize_offtargets(rows; distance = 2, variant_aware = true)
        end
    end
end
