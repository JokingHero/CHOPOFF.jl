using Random
isdefined(@__MODULE__, :vtest_random_case) || include(joinpath(@__DIR__, "..", "helpers", "vcf_sweeps.jl"))
isdefined(@__MODULE__, :VCF_TEST_REFERENCE) || include(joinpath(@__DIR__, "..", "sample_data", "vcf", "cases.jl"))
include(joinpath(@__DIR__, "..", "helpers", "vcf_genotypes.jl"))

@testset "VCF genotype mode [D01]" begin
    gt_seed = parse(Int, get(ENV, "CHOPOFF_VCF_GT_SEED", "20261014"))
    rng = MersenneTwister(gt_seed)
    r = VCF_TEST_REFERENCE
    snp(p, id) = vtest_edit(r, p, p + 1, vtest_alternative(r[p]), id)
    alt_guide(ps...) = (g = collect(VCF_TEST_GUIDE); foreach(p -> g[p - 40] = only(vtest_alternative(r[p])), ps); LongDNA{4}(String(g)))
    sample(name, calls...; ploidy = 2) = (; name, ploidy, calls = Dict{String, Any}(calls...))
    tokens(df, col) = [vtest_text(x) for x in df[!, col]]
    @testset "Cis, trans, unphased and split phase sets" begin
        edits = [snp(45, "a"), snp(50, "b")]
        samples = [sample("cis", "a" => ([1, 0], :chrom), "b" => ([1, 0], :chrom)),
            sample("trans", "a" => ([1, 0], :chrom), "b" => ([0, 1], :chrom)),
            sample("unph", "a" => ([0, 1], :unphased), "b" => ([0, 1], :unphased)),
            sample("ps", "a" => ([1, 0], 100), "b" => ([1, 0], 200))]
        mktempdir() do dir
            rows = vtest_compare_genotypes(dir, r, edits, samples, [alt_guide(45, 50)], Motif("Cas9"; distance = 0), 0)
            both = filter(row -> !ismissing(row.variants), rows)
            @test tokens(both, :carriers) == ["cis:1;unph:?;ps:?"]
            @test both.n_carriers == [1] && both.n_possible == [2]
        end
        mktempdir() do dir
            rows = vtest_compare_genotypes(dir, r, edits, samples, [alt_guide(45)], Motif("Cas9"; distance = 0), 0)
            @test tokens(filter(row -> !ismissing(row.variants), rows), :carriers) == ["trans:1;unph:?;ps:?"]
        end
        mktempdir() do dir
            rows = vtest_compare_genotypes(dir, r, edits, samples, [LongDNA{4}(VCF_TEST_GUIDE)], Motif("Cas9"; distance = 0), 0)
            reference = filter(row -> ismissing(row.variants), rows)
            @test tokens(reference, :lost) == ["cis:1;trans:1;trans:2;unph:*;unph:?;ps:*;ps:?"]
            @test reference.n_lost == [5] && reference.n_lost_possible == [2]
            vcf = joinpath(dir, "input.vcf")
            summary = summarize_offtargets(rows; distance = 0, per_sample = true, vcf)
            @test summary.sample == ["cis", "trans", "unph", "ps"]
            @test summary.D0 == [1, 0, 0, 0] && summary.P0 == [1, 0, 1, 1]
            detail, counts = joinpath(dir, "genotype.csv"), joinpath(dir, "per_sample.csv")
            CHOPOFF.main(["summarize", "--detail_file", detail, "--output", counts, "--per_sample",
                "--vcf", vcf, "--vcf_samples", "trans,cis", "--distance", "0"])
            @test DataFrame(CSV.File(counts)) == summary[1:2, :]
        end
        mktempdir() do dir
            rows = vtest_compare_genotypes(dir, r, edits, samples, [alt_guide(45, 50)], Motif("Cas9"; distance = 1), 1)
            summary = summarize_offtargets(rows; distance = 1, per_sample = true, vcf = joinpath(dir, "input.vcf"))
            # cis: both-ALT site at 0; trans: one-ALT site at 1; others possible only.
            @test summary.D0 == [1, 0, 0, 0] && summary.D1 == [0, 1, 0, 0]
            @test summary.P0 == [1, 0, 1, 1] && summary.P1 == [0, 1, 0, 0]
        end
    end
    @testset "A carried variant left of the anchor removes the reference-flank site" begin
        edits = [snp(45, "e"), snp(50, "a")]
        samples = [sample("hom", "e" => ([1, 1], :unphased), "a" => ([1, 1], :unphased))]
        mktempdir() do dir
            rows = vtest_compare_genotypes(dir, r, edits, samples, [alt_guide(50)], Motif("Cas9"; distance = 0), 0)
            @test all(ismissing, rows.variants)
            genome = joinpath(dir, "reference.fa")
            union = joinpath(dir, "union_only.csv")
            search_prefixHashScan([alt_guide(50)], genome, union; distance = 0, vcf = joinpath(dir, "input.vcf"))
            @test count(!ismissing, DataFrame(CSV.File(union)).variants) == 1
        end
    end
    @testset "Sample selection, CLI parity and input errors" begin
        edits = [snp(45, "a"), snp(50, "b")]
        samples = [sample("x", "a" => ([1, 0], :chrom), "b" => ([1, 0], :chrom)),
            sample("y", "a" => ([0, 1], :chrom), "b" => ([1, 0], :chrom))]
        mktempdir() do dir
            genome = vtest_fasta(joinpath(dir, "reference.fa"), ["chr1" => r])
            vcf = vtest_genotype_vcf(joinpath(dir, "input.vcf"), edits, samples)
            guides = joinpath(dir, "guides.txt")
            write(guides, VCF_TEST_GUIDE * "\n")
            api, cli = joinpath(dir, "api.csv"), joinpath(dir, "cli.csv")
            search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)], genome, api; distance = 1, vcf,
                vcf_mode = :genotype, vcf_samples = ["y"], scan_threads = 1)
            CHOPOFF.main(["search", "--guides", guides, "--output", cli, "--distance", "1",
                "prefixHashScan", "--genome", genome, "--vcf", vcf, "--vcf_mode", "genotype",
                "--vcf_samples", "y", "--scan_threads", "1"])
            @test read(cli) == read(api)
            df = DataFrame(CSV.File(api))
            @test names(df)[13:end] == ["carriers", "lost", "n_carriers", "n_possible", "n_lost", "n_lost_possible"]
            @test all(t -> !occursin("x:", t), vtest_text.(df.carriers))
            @test all(t -> !occursin("x:", t), vtest_text.(df.lost))
            message(; kw...) = try
                search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)], genome, api; distance = 1, kw...); ""
            catch err
                sprint(showerror, err)
            end
            @test occursin("absent from the VCF: z", message(; vcf, vcf_mode = :genotype, vcf_samples = ["z"]))
            @test occursin("requires vcf_mode=:genotype", message(; vcf, vcf_samples = ["x"]))
            @test occursin("vcf_mode must be", message(; vcf, vcf_mode = :haplotype))
            @test occursin("vcf_mode requires vcf", message(; vcf_mode = :genotype))
            sites = vtest_vcf(joinpath(dir, "sites.vcf"), vtest_records(edits))
            @test occursin("needs sample columns", message(; vcf = sites, vcf_mode = :genotype))
            nogt = vtest_vcf(joinpath(dir, "nogt.vcf"), [(row..., "DP", "3") for row in vtest_records(edits)];
                header = ["##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">"], samples = ["x"])
            @test occursin("FORMAT/GT", message(; vcf = nogt, vcf_mode = :genotype))
            deletion = vtest_edit(r, 44, 48, "", "d")
            # Overlapping alleles on one copy warn; the deletion hides the SNP.
            for (calls, warns) in ((((([1, 0], :chrom), ([1, 0], :chrom))), true),
                    ((([0, 1], :unphased), ([0, 1], :unphased)), false),
                    ((([1, 1], :unphased), ([0, 1], :unphased)), true))
                conflict = vtest_genotype_vcf(joinpath(dir, "conflict.vcf"), [deletion, snp(45, "a")],
                    [sample("c", "d" => calls[1], "a" => calls[2])])
                if warns
                    @test_logs (:warn, r"overlapping alleles on one copy") match_mode = :any search_prefixHashScan(
                        [LongDNA{4}(VCF_TEST_GUIDE)], genome, api; distance = 1, vcf = conflict, vcf_mode = :genotype)
                else
                    @test isempty(message(; vcf = conflict, vcf_mode = :genotype))
                end
            end
            # Phased on one copy: rows of the SNP alone never credit that copy.
            conflict = vtest_genotype_vcf(joinpath(dir, "conflict.vcf"), [deletion, snp(45, "a")],
                [sample("c", "d" => ([1, 0], :chrom), "a" => ([1, 0], :chrom))])
            snp_guide = (g = collect(VCF_TEST_GUIDE); g[5] = only(vtest_alternative(r[45])); LongDNA{4}(String(g)))
            search_prefixHashScan([snp_guide], genome, api; distance = 0, vcf = conflict, vcf_mode = :genotype)
            @test all(ismissing, DataFrame(CSV.File(api)).variants)
        end
    end
    @testset "Random genotype case $i" for i in 1:parse(Int, get(ENV, "CHOPOFF_VCF_GT_CASES", "20"))
        mktempdir() do dir
            r, edits, guides, motif, d = vtest_random_case(rng, i)
            samples = vtest_random_samples(rng, r, edits)
            write(joinpath(dir, "reproduce.txt"), "CHOPOFF_VCF_GT_SEED=$gt_seed CHOPOFF_VCF_GT_CASES=$i\n" *
                repr((r, edits, guides, motif, samples)))
            vtest_compare_genotypes(dir, r, edits, samples, guides, motif, d; case_id = "gt$i",
                chunk = isodd(i) ? 64 : 2 * 1024 * 1024)
        end
    end
end
