@testset "VCF prefixHashScan" begin
    mktempdir() do dir
        guide = dna"ACGTGTCAGTACGATCGTAC"
        reference = repeat("T", 40) * string(guide) * "AGG" * repeat("T", 80)
        motif = Motif("Cas9"; distance = 2)

        @testset "independent haplotype oracle" begin
            edits = [vtest_edit(reference, 45, 46, "A", "snp"),
                vtest_edit(reference, 51, 51, "C", "ins"),
                vtest_edit(reference, 54, 55, "", "del")]
            vtest_compare(dir, reference, edits, [guide], motif, 2; chunk = 64)
            overlapping = [vtest_edit(reference, 44, 48, "", "del1"),
                vtest_edit(reference, 45, 46, "A", "inside"),
                vtest_edit(reference, 47, 50, "", "del2")]
            vtest_compare(dir, reference, overlapping, [guide], Motif("Cas9"; distance = 4), 4)
            boundary = [vtest_edit(reference, 46, 46, "C", "gap"),
                vtest_edit(reference, 46, 47, "G", "base"),
                vtest_edit(reference, 48, 49, "", "adj1"),
                vtest_edit(reference, 49, 50, "", "adj2")]
            vtest_compare(dir, reference, boundary, [guide], motif, 2; chunk = 64)
            replacement = [vtest_edit(reference, 46, 48, "CA", "mnp"),
                vtest_edit(reference, 53, 55, "CAA", "replace")]
            vtest_compare(dir, reference, replacement, [guide], Motif("Cas9"; distance = 4), 4)
        end

        @testset "PAM and IUPAC choices" begin
            edit = vtest_edit(reference, 62, 63, "A", "pam")
            rows = vtest_compare(dir, reference, [edit], [guide], motif, 2)
            @test any(!ismissing, rows.overlapping_variants)
            @test all(ismissing, rows.variants)
            ambiguous = [vtest_edit(reference, 45, 46, "A", "amb"; original = "R", group = "amb"),
                vtest_edit(reference, 45, 46, "G", "amb"; original = "R", group = "amb")]
            vtest_compare(dir, reference, ambiguous, [guide], motif, 2)
            ambiguous = [vtest_edit(reference, 45, 46, string(base), "any"; original = "N", group = "any") for base in "ACGT"]
            vtest_compare(dir, reference, ambiguous, [guide], motif, 2)
            ambiguous_reference = reference[1:42] * "N" * reference[44:end]
            for budget in 0:1
                vtest_compare(dir, ambiguous_reference,
                    [vtest_edit(ambiguous_reference, 45, 46, "A", "near_reference_N")],
                    [guide], Motif("Cas9"; distance = 2, ambig_max = budget), 2)
            end
            for distance in 0:4
                altguide = copy(guide)
                altguide[5] = DNA_A
                vtest_compare(dir, reference, [vtest_edit(reference, 45, 46, "A", "snp")],
                    [altguide], Motif("Cas9"; distance = distance), distance)
            end
        end

        @testset "long insertions, both strands and Cas12a" begin
            inserted = string(guide) * "AGG" * repeat("T", 35) * string(guide) * "AGG"
            rows = vtest_compare(dir, reference, [vtest_edit(reference, 90, 90, inserted, "long")],
                [guide], motif, 2; chunk = 64, threads = 2)
            inside = rows[.!ismissing.(rows.start_origin), :]
            @test length(unique(inside.start_origin)) >= 2
            filtered = filter_overlapping(copy(rows), 2)
            @test summarize_offtargets(rows; distance = 2, variant_aware = true) ==
                summarize_offtargets(filtered; distance = 2, variant_aware = true)
            reversed = string(reverse_complement(LongDNA{4}(reference)))
            vtest_compare(dir, reversed, [vtest_edit(reversed, 92, 93, "A", "reverse")], [guide], motif, 2)
            guide12 = dna"ACGTGTCAGTACGATCGTACA"
            ref12 = repeat("G", 40) * "TTTA" * string(guide12) * repeat("G", 60)
            vtest_compare(dir, ref12, [vtest_edit(ref12, 48, 49, "C", "snp12"),
                vtest_edit(ref12, 53, 53, "A", "ins12")], [guide12], Motif("Cas12a"; distance = 2), 2; chunk = 64)
        end

        @testset "validation, annotations, compression and chromosomes" begin
            genome = vtest_fasta(joinpath(dir, "two.fa"), ["chr1" => reference, "chr2" => reference])
            vcf = joinpath(dir, "contract.vcf")
            output = joinpath(dir, "contract.csv")
            run(; kwargs...) = search_prefixHashScan([guide], genome, output; distance = 1, vcf = vcf, kwargs...)
            vtest_vcf(vcf, [("wrong", 45, ".", "G", "A", ".", ".", ".")])
            @test_throws ErrorException run()
            vtest_vcf(vcf, [("chr1", 45, ".", "T", "A", ".", ".", ".")])
            @test_throws ErrorException run()
            vtest_vcf(vcf, [vtest_record(vtest_edit(reference, 50, 51, "A", "later")),
                vtest_record(vtest_edit(reference, 45, 46, "A", "earlier"))])
            @test_throws ErrorException run()
            vtest_vcf(vcf, [("chr1", 45, ".", reference[45:45], "<DEL>", ".", ".", ".")])
            @test_throws ErrorException run()
            vtest_vcf(vcf, [])
            @test_throws ErrorException run(output = :counts)
            @test_throws ErrorException run(region = "chr1:1-100")
            @test_throws ErrorException run(vcf_annotate = ["INFO/absent"])
            run(region = "chr1")
            @test all(==("chr1"), DataFrame(CSV.File(output)).chromosome)
            base = joinpath(dir, "base.csv")
            search_prefixHashScan([guide], genome, base; distance = 1, region = "chr1")
            @test select(DataFrame(CSV.File(output)), names(DataFrame(CSV.File(base)))) == DataFrame(CSV.File(base))
            vtest_vcf(vcf, [("chr1", 45, "multi", reference[45:45], "A,C", ".", ".", "AF=0.1,0.2;NOTE=a,b")];
                header = ["##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">",
                    "##INFO=<ID=NOTE,Number=.,Type=String,Description=\"Text\">"])
            run(vcf_annotate = ["INFO/AF", "INFO/NOTE"], region = "chr1")
            rows = DataFrame(CSV.File(output))
            @test any(v -> !ismissing(v) && occursin("INFO/AF=0.1", v) && occursin("INFO/NOTE=a,b", v), rows.variants)
            @test any(v -> !ismissing(v) && occursin("INFO/AF=0.2", v), rows.variants)
            zipped = vcf * ".gz"
            open(zipped, "w") do io
                stream = GzipCompressorStream(io)
                write(stream, read(vcf))
                close(stream)
            end
            zipout = joinpath(dir, "compressed.csv")
            search_prefixHashScan([guide], genome, zipout; distance = 1, vcf = zipped,
                vcf_annotate = ["INFO/AF", "INFO/NOTE"], region = "chr1")
            @test read(output) == read(zipout)
            # BGZF consists of concatenated gzip members. Split within a row
            # so the parser must see one continuous decompressed stream.
            payload = read(vcf)
            middle = length(payload) - 10
            open(zipped, "w") do io
                for block in (payload[1:middle], payload[(middle + 1):end])
                    write(io, transcode(GzipCompressor, block))
                end
            end
            search_prefixHashScan([guide], genome, zipout; distance = 1, vcf = zipped,
                vcf_annotate = ["INFO/AF", "INFO/NOTE"], region = "chr1")
            @test read(output) == read(zipout)
            CHOPOFF.search_prefixHashScan([guide], genome, Motif("Cas9"; distance = 1), base;
                distance = 1, region = "chr2", scan_backend = :legacy)
            @test all(==("chr2"), DataFrame(CSV.File(base)).chromosome)
            search_prefixHashScan([guide], genome, base; distance = 1, region = "chr1", output = :counts)
            @test DataFrame(CSV.File(base)).D0 == [1]
        end

        @testset "caps and filter semantics" begin
            genome = vtest_fasta(joinpath(dir, "cap.fa"), ["chr1" => reference])
            vcf = vtest_vcf(joinpath(dir, "cap.vcf"),
                [vtest_record(vtest_edit(reference, 45, 46, "A", "variant"))])
            output = joinpath(dir, "cap.csv")
            search_prefixHashScan([guide], genome, output; distance = 1, vcf = vcf, early_stopping = [1, 1])
            rows = DataFrame(CSV.File(output))
            @test any(ismissing, rows.variants)
            @test any(!ismissing, rows.variants)
            # Global locus assignment must precede per-allele filtering.
            df = DataFrame(guide = fill("g", 3), chromosome = fill("chr1", 3), strand = fill("+", 3),
                start = [10, 13, 7], distance = [1, 2, 0], variants = ["A", "A", "B"],
                start_origin = fill("", 3), alignment_guide = fill("g", 3), alignment_reference = fill("a", 3))
            filtered = filter_overlapping(copy(df), 3)
            @test nrow(filtered) == 3
            @test summarize_offtargets(df; distance = 3, variant_aware = true) ==
                summarize_offtargets(filtered; distance = 3, variant_aware = true)
            @test sum(Matrix(select(summarize_offtargets(df; distance = 3, variant_aware = true), Not(:guide)))) == 2
        end

        @testset "deletion reachability and chromosome edges" begin
            distant = repeat("T", 40) * string(guide)[1:10] * repeat("C", 70) * string(guide)[11:20] * "AGG" * repeat("T", 40)
            rows = vtest_compare(dir, distant, [vtest_edit(distant, 46, 47, "A", "distant_snp"),
                vtest_edit(distant, 51, 121, "", "long_del")], [guide], motif, 2; chunk = 64)
            @test any(v -> !ismissing(v) && occursin("distant_snp", v) && occursin("long_del", v), rows.variants)
            edge = string(guide) * "AGG" * repeat("T", 50)
            vtest_compare(dir, edge, [vtest_edit(edge, 1, 1, "C", "before_first"),
                vtest_edit(edge, 1, 2, "", "delete_first"),
                vtest_edit(edge, 3, 4, "A", "near_edge")], [guide], motif, 2)
            edge = repeat("T", 50) * string(guide) * "AGG"
            vtest_compare(dir, edge, [vtest_edit(edge, 71, 72, "G", "last_pam")], [guide], motif, 2)
        end

        @testset "sibling ALTs and neutral allele combinations" begin
            genome = vtest_fasta(joinpath(dir, "siblings.fa"), ["chr1" => reference])
            vcf = vtest_vcf(joinpath(dir, "siblings.vcf"),
                [("chr1", 45, "siblings", "G", "A,GC", ".", ".", ".")])
            output = joinpath(dir, "siblings.csv")
            search_prefixHashScan([guide], genome, output; distance = 2, vcf = vcf)
            rows = DataFrame(CSV.File(output))
            @test any(!ismissing, rows.variants)
            @test all(v -> ismissing(v) || !occursin(';', v), rows.variants)
            # Two independently represented changes can reconstruct reference.
            edit = [vtest_edit(reference, 47, 48, "", "remove"),
                vtest_edit(reference, 47, 47, reference[47:47], "restore")]
            rows = vtest_compare(dir, reference, edit, [guide], motif, 2)
            @test any(row -> !ismissing(row.variants) && occursin("remove", row.variants) &&
                occursin("restore", row.variants) && row.distance == 0, eachrow(rows))
        end

        @testset "2bit, portable SIMD, annotations and batching" begin
            genome = vtest_fasta(joinpath(dir, "formats.fa"), ["chr1" => reference])
            twobit = joinpath(dir, "formats.2bit")
            writer = TwoBit.Writer(open(twobit, "w"), ["chr1"])
            write(writer, TwoBit.Record("chr1", LongDNA{4}(reference)))
            close(writer)
            vcf = vtest_vcf(joinpath(dir, "formats.vcf"),
                [("chr1", 45, ".", "G", "A", ".", ".", ".", "AF", "0.25")];
                    header = ["##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Fraction\">"], samples = ["sample%1"])
            output = joinpath(dir, "formats.csv")
            search_prefixHashScan([guide], genome, output; distance = 1, vcf = vcf, vcf_annotate = ["FORMAT/AF"])
            expected = read(output)
            search_prefixHashScan([guide], twobit, output; distance = 1, vcf = vcf,
                vcf_annotate = ["FORMAT/AF"], simd_backend = :portable)
            @test read(output) == expected
            @test occursin("FORMAT/AF=sample%251=0.25", String(copy(expected)))
            guides = [LongDNA{4}(join("ACGT"[1 + ((i >> (2 * j)) & 3)] for j in 0:3) * string(guide)[5:end]) for i in 0:64]
            search_prefixHashScan(guides, genome, output; distance = 1, vcf = vcf)
            @test count(line -> startswith(line, "guide,"), readlines(output)) == 1
            @test nrow(DataFrame(CSV.File(output))) > 0
            stats = PrefixHashScanStats()
            CHOPOFF.search_prefixHashScan([guide], genome, Motif("Cas9"; distance=1), output;
                distance = 1, vcf = vcf, stats = stats)
            @test stats.variants_read == 1
            @test stats.fragments == 1
            @test stats.variant_rows > 0
            @test stats.descriptor_bytes > 0
            @test_throws ErrorException CHOPOFF.search_prefixHashScan([guide], genome, motif, output;
                distance = 2, vcf = vcf, scan_backend = :legacy)
            write(output, "preserve me")
            @test_throws ErrorException search_prefixHashScan([guide], genome, output; vcf = vcf, vcf_annotate = ["INFO/absent"])
            @test read(output, String) == "preserve me"
        end

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
