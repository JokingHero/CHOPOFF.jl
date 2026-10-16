isdefined(@__MODULE__, :vfixture_cases) || include(joinpath(@__DIR__, "..", "sample_data", "vcf", "cases.jl"))
include(joinpath(@__DIR__, "..", "sample_data", "vcf", "make_test_vcfs.jl"))

@testset "VCF fixture corpus [V04-V11,V16]" begin
    @testset "Reproducible files" begin
        @test isnothing(check_vcf_fixtures())
    end
    @testset "Hand-calculated oracle checks" begin
        s = "AACCGG"
        snp = vtest_edit(s,2,3,"T","snp")
        seq, coords, events = vtest_haplotype(s,[snp])
        @test seq == "ATCCGG"
        @test vtest_origin.(coords) == [(i,"") for i in 1:6]
        @test vtest_annotation(snp,"chr1") == "snp|chr1:2:A>T"
        insertion = vtest_edit(s,3,3,"T","ins")
        seq,coords,events = vtest_haplotype(s,[insertion])
        @test seq == "AATCCGG"
        @test vtest_origin.(coords) == [(1,""),(2,""),(2,"chr1:2:A>AT:1"),(3,""),(4,""),(5,""),(6,"")]
        deletion = vtest_edit(s,3,5,"","del")
        seq,coords,events = vtest_haplotype(s,[deletion])
        @test seq == "AAGG"
        @test vtest_origin.(coords) == [(1,""),(2,""),(5,""),(6,"")]
        @test (events[1][1],events[1][2]) == (3,2)
        @test vtest_haplotype(s,[deletion,vtest_edit(s,4,5,"T","inside")]) === nothing
        @test length(vtest_haplotypes(s,[deletion,insertion])) == 4
    end
    root = joinpath(@__DIR__, "..", "sample_data", "vcf")
    for case in vfixture_cases()
        @testset "$(case.id): $(case.purpose)" begin
            mktempdir() do dir
                guides = LongDNA{4}.(readlines(joinpath(root,"guides",case.id*".txt")))
                @test guides == case.guides
                @test length(vtest_haplotypes(case.reference,case.edits)) <= 64
                rows = vtest_compare(dir,case.reference,case.edits,guides,vfixture_motif(case),case.distance;
                    chunk=64, threads=min(4,Threads.nthreads()), case_id=case.id,
                    genome=joinpath(root,"references",case.id*".fa"),vcf=joinpath(root,"valid",case.id*".vcf"))
                if case.variant
                    @test any(eachrow(rows)) do row
                        ids = sort([first(split(token,'|')) for token in vtest_tokens(row.variants)])
                        ids == sort(case.witness)
                    end
                else
                    @test !isempty(rows.distance)
                    @test all(ismissing,rows.variants)
                    @test any(!ismissing,rows.overlapping_variants)
                end
                if startswith(case.id,"long_insertion")
                    @test any(row -> !ismissing(row.start_origin) && !ismissing(row.end_origin),eachrow(rows))
                    filtered = filter_overlapping(copy(rows),case.distance)
                    @test summarize_offtargets(rows;distance=case.distance,variant_aware=true) ==
                        summarize_offtargets(filtered;distance=case.distance,variant_aware=true)
                elseif case.id == "endpoint_insertion"
                    @test any(row -> !ismissing(row.start_origin) && ismissing(row.end_origin),eachrow(rows))
                elseif case.id == "pam_create"
                    @test all(!ismissing,rows.variants)
                end
            end
        end
    end
end
