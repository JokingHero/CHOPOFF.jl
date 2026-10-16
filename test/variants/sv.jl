using Random
include(joinpath(@__DIR__, "..", "helpers", "vcf_sv.jl"))

@testset "Symbolic structural variants [D02a]" begin
    rng = MersenneTwister(20261015)
    g = VCF_TEST_GUIDE
    motif0 = Motif("Cas9"; distance = 0)
    @testset "Long DEL bridges a site, with a flank SNP" begin
        mktempdir() do dir
            left = g[1:2] * "T" * g[4:10]          # base 3 needs the SNP
            r = "T"^40 * left * vtest_sv_random(rng, 12_000) * g[11:20] * "AGG" * "T"^40
            del = vtest_sv(r, "DEL", 50, 50 + 12_000, "del")
            rows, stats = vtest_sv_compare(dir, r, [del], [LongDNA{4}(g)], Motif("Cas9"; distance = 1), 1)
            @test stats.vcf_sv_searched == 1
            @test any(row -> !ismissing(row.variants) && occursin("<DEL>/12050", row.variants), eachrow(rows))
        end
        # Exact sites need the SNP and the deletion together.
        mktempdir() do dir
            left = g[1:2] * "T" * g[4:10]
            r = "T"^40 * left * vtest_sv_random(rng, 12_000) * g[11:20] * "AGG" * "T"^40
            del = vtest_sv(r, "DEL", 50, 50 + 12_000, "del")
            snp = vtest_edit(r, 43, 44, string(g[3]), "snp")
            svs = [del, (snp, vtest_record(snp), vtest_key(snp, "chr1"))]
            rows = first(vtest_sv_compare(dir, r, svs, [LongDNA{4}(g)], motif0, 0))
            @test any(row -> row.distance == 0 && length(vtest_tokens(row.variants)) == 2, eachrow(rows))
        end
    end
    @testset "Tandem DUP longer than the chunk creates a junction site" for alt in ("<DUP>", "<DUP:TANDEM>")
        mktempdir() do dir
            region = g[11:20] * "AGG" * vtest_sv_random(rng, 3_000) * g[1:10]
            r = "T"^40 * region * "T"^40
            dup = vtest_sv(r, "DUP", 40, 40 + length(region), "dup"; alt)
            rows, stats = vtest_sv_compare(dir, r, [dup], [LongDNA{4}(g)], motif0, 0)
            @test stats.vcf_sv_searched == 1
            @test any(row -> row.distance == 0 && !ismissing(row.variants), eachrow(rows))
        end
    end
    @testset "INV longer than the chunk restores a site" begin
        mktempdir() do dir
            region = "A" * vtest_rc(g * "AGG") * vtest_sv_random(rng, 2_000) * "A"
            r = "T"^40 * region * "T"^40
            inv = vtest_sv(r, "INV", 40, 40 + length(region), "inv")
            rows, _ = vtest_sv_compare(dir, r, [inv], [LongDNA{4}(g)], motif0, 0)
            @test any(row -> row.distance == 0 && !ismissing(row.variants), eachrow(rows))
        end
    end
    @testset "Upstream walk taking a long DUP (segment truncation)" begin
        mktempdir() do dir
            region = g[11:20] * "AGG" * vtest_sv_random(rng, 2_000) * g[1:10]
            r = "T"^40 * region * "T"^40
            dup = vtest_sv(r, "DUP", 40, 40 + length(region), "dup")
            snp_edit = vtest_edit(r, 35, 36, "A", "snp")
            snp = (snp_edit, vtest_record(snp_edit), vtest_key(snp_edit, "chr1"))
            vtest_sv_compare(dir, r, [dup, snp], [LongDNA{4}(g)], Motif("Cas9"; distance = 2), 2)
        end
    end
    @testset "Deep DUP copy sites are their own locus" begin
        mktempdir() do dir
            region = vtest_sv_random(rng, 30) * g * "AGG" * vtest_sv_random(rng, 30)
            r = "T"^40 * region * "T"^40
            dup = vtest_sv(r, "DUP", 40, 40 + length(region), "dup")
            rows, _ = vtest_sv_compare(dir, r, [dup], [LongDNA{4}(g)], motif0, 0)
            @test any(row -> !ismissing(row.start_origin) && occursin("<DUP>/", row.start_origin), eachrow(rows))
            @test summarize_offtargets(rows; distance = 0, variant_aware = true).D0 == [2]
            @test nrow(filter_overlapping(copy(rows), 0)) == 2
        end
    end
    @testset "SVLEN-only and multiallelic records" begin
        mktempdir() do dir
            r = "T"^40 * g[1:10] * vtest_sv_random(rng, 500) * g[11:20] * "AGG" * "T"^40
            del = vtest_sv(r, "DEL", 50, 550, "del"; info = "SVLEN=-500")
            vtest_sv_compare(dir, r, [del], [LongDNA{4}(g)], motif0, 0)
        end
        mktempdir() do dir
            r = "T"^40 * g[1:10] * vtest_sv_random(rng, 500) * g[11:20] * "AGG" * "T"^40
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            vcf = vtest_vcf(joinpath(dir, "multi.vcf"), [("chr1", 50, "m", r[50:50], "<DEL>,<DUP>", ".", "PASS", "SVLEN=-500,500")])
            stats = PrefixHashScanStats()
            out = joinpath(dir, "multi.csv")
            CHOPOFF.search_prefixHashScan([LongDNA{4}(g)], genome, motif0, out; distance = 0, vcf, stats)
            rows = DataFrame(CSV.File(out))
            @test stats.vcf_sv_searched == 2
            @test any(row -> !ismissing(row.variants) && occursin("<DEL>/550", row.variants), eachrow(rows))
        end
    end
    @testset "Skipped kinds are counted, not errors" begin
        mktempdir() do dir
            r = VCF_TEST_REFERENCE
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            row(pos, alt, info) = ("chr1", pos, "s$pos", r[pos:pos], alt, ".", "PASS", info)
            records = [row(10, "<INS>", "SVLEN=100"), row(11, "<CNV>", "END=20"),
                row(12, "T]chr1:30]", "SVTYPE=BND"), row(13, "<DEL>", "END=20;IMPRECISE"),
                row(14, "<DUP:INT>", "END=20"), row(15, "<FOO>", "END=20"), row(16, "T.", "SVTYPE=BND"),
                row(17, "<DEL:ME:ALU>", "END=20"), row(18, "<DEL>", "END=20;CIPOS=-5,5")]
            vcf = vtest_vcf(joinpath(dir, "skip.vcf"), records)
            stats = PrefixHashScanStats()
            out = joinpath(dir, "skip.csv")
            CHOPOFF.search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)], genome, Motif("Cas9"; distance = 1), out;
                distance = 1, vcf, stats)
            @test stats.vcf_sv_skipped == 7
            @test stats.vcf_sv_searched == 2
        end
    end
    @testset "Malformed meta lines are skipped" begin
        mktempdir() do dir
            r = VCF_TEST_REFERENCE
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            header = ["###INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strands [++,+-]\">",
                "##not a key value line", "##ALT=<ID=DEL,Description=\"Deletion\"",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT frequency\">"]
            vcf = vtest_vcf(joinpath(dir, "meta.vcf"), [("chr1", 45, "snp", "G", "A", ".", "PASS", "AF=0.5")]; header)
            out = joinpath(dir, "meta.csv")
            CHOPOFF.search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)], genome, Motif("Cas9"; distance = 1), out;
                distance = 1, vcf, vcf_annotate = ["INFO/AF"])
            rows = DataFrame(CSV.File(out))
            @test any(v -> !ismissing(v) && occursin("INFO/AF=0.5", v), rows.variants)
        end
    end
    @testset "Malformed symbolic records fail" begin
        mktempdir() do dir
            r = VCF_TEST_REFERENCE
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            for (info, message) in (("SVTYPE=DEL", "needs integer INFO/END"), ("END=10", "ends outside"),
                    ("END=$(length(r) + 1)", "ends outside"))
                vcf = vtest_vcf(joinpath(dir, "bad.vcf"), [("chr1", 10, "bad", r[10:10], "<DEL>", ".", "PASS", info)])
                err = try
                    CHOPOFF.search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)], genome, Motif("Cas9"; distance = 1),
                        joinpath(dir, "bad.csv"); distance = 1, vcf)
                    nothing
                catch e
                    e
                end
                @test err isa ErrorException && occursin(message, err.msg)
            end
        end
    end
    @testset "Genotype mode: het DEL carriers and lost reference sites" begin
        mktempdir() do dir
            r = "T"^40 * g[1:10] * "C"^30 * g * "AGG" * "T"^20 * g[11:20] * "AGG" * "T"^40
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            stop = 50 + 30 + 23 + 20
            vcf = vtest_vcf(joinpath(dir, "gt.vcf"), [("chr1", 50, "del", r[50:50], "<DEL>", ".", "PASS", "END=$stop", "GT", "0|1")];
                samples = ["S1"], header = ["##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"])
            out = joinpath(dir, "gt.csv")
            CHOPOFF.search_prefixHashScan([LongDNA{4}(g)], genome, motif0, out; distance = 0, vcf, vcf_mode = :genotype)
            rows = DataFrame(CSV.File(out))
            variant = filter(row -> !ismissing(row.variants), rows)
            reference = filter(row -> ismissing(row.variants), rows)
            @test nrow(variant) >= 1 && all(==("S1:2"), variant.carriers)
            @test nrow(reference) == 1 && only(reference.lost) == "S1:2"
        end
    end
end
