@testset "VCF execution and output [V12-V15]" begin
    @testset "Formats, SIMD, chunks and workers" begin
        backends = [:portable, :auto]
        CHOPOFF.can_use_prefix_hash_scan_avx2() && push!(backends, :avx2)
        CHOPOFF.can_use_prefix_hash_scan_avx512() && push!(backends, :avx512)
        @info "VCF execution matrix" backends threads=Threads.nthreads()
        for id in ("mixed_choices", "long_insertion_both")
            case = only(filter(c -> c.id == id, vfixture_cases()))
            mktempdir() do dir
                motif = vfixture_motif(case)
                guides = LongDNA{4}.(case.guides)
                expected = vtest_oracle(dir, case.reference, case.edits, guides, motif, case.distance)
                vcf = vtest_vcf(joinpath(dir,"input.vcf"), vtest_records(case.edits))
                genomes = [vtest_fasta(joinpath(dir,"wrap$width.fa"), ["chr1" => case.reference]; width)
                    for width in (7,60)]
                twobit = joinpath(dir,"reference.2bit")
                writer = TwoBit.Writer(open(twobit,"w"), ["chr1"])
                write(writer, TwoBit.Record("chr1",LongDNA{4}(case.reference),nothing))
                close(writer)
                push!(genomes,twobit)
                baseline = nothing
                for genome in genomes, backend in backends, chunk in (64,97,2*1024*1024),
                        threads in unique([1,min(4,Threads.nthreads())])
                    @testset "$id $(basename(genome)) $backend chunk=$chunk workers=$threads" begin
                        out = joinpath(dir,"actual.csv")
                        CHOPOFF.search_prefixHashScan(guides,genome,motif,out;
                            distance=case.distance,vcf,simd_backend=backend,stream_chunk_bases=chunk,
                            scan_threads=threads,early_stopping=fill(typemax(Int),case.distance+1))
                        vtest_assert_rows(DataFrame(CSV.File(out)),expected,dir;case_id=id)
                        baseline === nothing && (baseline = read(out))
                        @test read(out) == baseline
                    end
                end
            end
        end
    end
    @testset "Guide batches at 64, 65 and 129" begin
        mktempdir() do dir
            g = VCF_TEST_GUIDE
            reference = VCF_TEST_REFERENCE
            anchors = Dict(1=>g,64=>"T"*g[2:end],65=>"C"*g[2:4]*"A"*g[6:end],
                129=>"G"*g[2:4]*"A"*g[6:end])
            pool = [join("ACGT"[1+((i >> (2*j)) & 3)] for j in 0:3)*g[5:end] for i in 0:255]
            filter!(x -> !(x in values(anchors)),pool)
            guides = [LongDNA{4}(haskey(anchors,i) ? anchors[i] : popfirst!(pool)) for i in 1:129]
            @test length(unique(guides)) == 129
            genome = vtest_fasta(joinpath(dir,"reference.fa"),["chr1"=>reference])
            vcf = vtest_vcf(joinpath(dir,"input.vcf"),vtest_records([vtest_edit(reference,45,46,"A","snp")]))
            for n in (64,65,129)
                out = joinpath(dir,"batch.csv")
                search_prefixHashScan(guides[1:n],genome,out;distance=2,vcf,
                    early_stopping=fill(typemax(Int),3))
                expected = Set{NamedTuple}()
                for first in 1:32:n
                    single = joinpath(dir,"small_batch.csv")
                    CHOPOFF.search_prefixHashScan(guides[first:min(first+31,n)],genome,
                        Motif("Cas9";distance=2),single;distance=2,vcf,early_stopping=fill(typemax(Int),3))
                    union!(expected, (vtest_signature(row) for row in eachrow(DataFrame(CSV.File(single)))))
                end
                rows = vtest_assert_rows(DataFrame(CSV.File(out)),expected,dir;case_id="batch$n")
                @test count(line -> startswith(line,"guide,"),readlines(out)) == 1
                @test any(row -> row.guide == string(guides[n]) && !ismissing(row.variants),eachrow(rows))
            end
            stats = PrefixHashScanStats()
            CHOPOFF.search_prefixHashScan(guides[1:1],genome,Motif("Cas9";distance=1),joinpath(dir,"stats.csv");
                distance=1,vcf,stats)
            @test stats.variants_read == 1
            @test stats.fragments == 1
            @test stats.variant_rows > 0
            @test stats.descriptor_bytes > 0
        end
    end
    @testset "Exact filtering and locus summaries" begin
        df = DataFrame(guide=fill("g",4),chromosome=fill("chr1",4),strand=fill("+",4),
            start=[10,13,7,8],distance=[1,2,0,2],variants=["A","A","B","A"],
            start_origin=fill("",4),alignment_guide=fill("g",4),alignment_reference=fill("a",4))
        expected = DataFrame(guide=["g"],D0=[1],D1=[0],D2=[1],D3=[0])
        filtered = filter_overlapping(copy(df),3)
        @test Set(zip(filtered.start,filtered.variants)) == Set([(7,"B"),(10,"A"),(13,"A")])
        @test summarize_offtargets(df;distance=3,variant_aware=true) == expected
        @test summarize_offtargets(filtered;distance=3,variant_aware=true) == expected
        @test isequal(filter_overlapping(copy(filtered),3),filtered)
        for order in ([4,3,2,1],[2,4,1,3])
            @test isequal(filter_overlapping(df[order,:],3),filtered)
        end
        df.start .= 100
        df.distance = [0,0,1,0]
        df.variants = ["","A","A","B"]
        # A 10-base insertion: offset 1 is within distance of its end and joins
        # the reference locus at POS; offsets 5 and 6 form one locus inside it.
        key = "chr1:100:C>CAAAAAAAAAA"
        df.start_origin = ["","$key:1","$key:5","$key:6"]
        @test summarize_offtargets(df;distance=1,variant_aware=true) == DataFrame(guide=["g"],D0=[2],D1=[0])
        @test nrow(filter_overlapping(copy(df),1)) == 4
        # The last inserted base also joins POS; offset 5 (d1) is then alone.
        df.start_origin[4] = "$key:10"
        @test summarize_offtargets(df;distance=1,variant_aware=true) == DataFrame(guide=["g"],D0=[1],D1=[1])
        @test nrow(filter_overlapping(df[1:0,:],2)) == 0
        empty = summarize_offtargets(df[1:0,:];distance=2,variant_aware=true)
        @test names(empty) == ["guide","D0","D1","D2"]
        @test nrow(empty) == 0
    end
    @testset "Search, filter and summarize CLI" begin
        root = joinpath(@__DIR__,"..","sample_data","vcf")
        genome = joinpath(root,"references","contracts.fa")
        vcf = joinpath(root,"valid","annotations.vcf")
        guidefile = joinpath(root,"guides","mixed_choices.txt")
        mktempdir() do dir
            api,cli,filtered,counts = [joinpath(dir,name*".csv") for name in ("api","cli","filtered","counts")]
            search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)],genome,api;distance=1,vcf,
                region="chr1",vcf_annotate=["INFO/AF","FORMAT/AD"],scan_threads=1)
            CHOPOFF.main(["search","--guides",guidefile,"--output",cli,"--distance","1",
                "prefixHashScan","--genome",genome,"--vcf",vcf,"--region","chr1",
                "--vcf_annotate","INFO/AF,FORMAT/AD","--scan_threads","1"])
            @test read(cli) == read(api)
            CHOPOFF.main(["filter","--detail_file",cli,"--output",filtered,"--distance","1"])
            CHOPOFF.main(["summarize","--detail_file",filtered,"--output",counts,"--variant_aware","--distance","1"])
            @test isequal(DataFrame(CSV.File(filtered)),filter_overlapping(DataFrame(CSV.File(api)),1))
            @test DataFrame(CSV.File(counts)) == DataFrame(guide=[VCF_TEST_GUIDE],D0=[1],D1=[0])
            nohits = vtest_fasta(joinpath(dir,"nohits.fa"),["chr1"=>"T"^100])
            emptyvcf = vtest_vcf(joinpath(dir,"empty.vcf"),[])
            search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)],nohits,api;distance=1,vcf=emptyvcf)
            rows = DataFrame(CSV.File(api))
            @test nrow(rows) == 0
            @test length(names(rows)) == 12
        end
    end
    @testset "Verbose progress covers VCF read, enumeration and scan" begin
        root = joinpath(@__DIR__,"..","sample_data","vcf")
        mktempdir() do dir
            logs, _ = Test.collect_test_logs() do
                search_prefixHashScan([LongDNA{4}(VCF_TEST_GUIDE)],joinpath(root,"references","contracts.fa"),
                    joinpath(dir,"out.csv");distance=1,vcf=joinpath(root,"valid","annotations.vcf"),
                    verbose=true,progress_interval=0)
            end
            labels = Set(get(Dict(l.kwargs),:batch,"") for l in logs if l.message == "prefixHashScan progress")
            @test issubset(Set(["VCF read","1/1 VCF enumerate","1/1 VCF"]),labels)
        end
    end
    @testset "A failure during the variant scan keeps the previous output [V14]" begin
        root = joinpath(@__DIR__,"..","sample_data","vcf")
        mktempdir() do dir
            g = VCF_TEST_GUIDE
            guides = [LongDNA{4}(join("ACGT"[1+((i >> (2*j)) & 3)] for j in 0:3)*g[5:end]) for i in 0:64]
            out = joinpath(dir,"out.csv")
            write(out,"preserve existing output")
            before = readdir(dir)
            # Batch 1 completes; batch 2 fails after its reference rows reach the staged file.
            calls = Ref(0)
            CHOPOFF.eval(:(stream_prefix_hash_scan_variants!(out::IOStream, args...; kwargs...) =
                ($calls[] += 1) == 2 ? (write(out,"partial\n"); error("injected variant scan failure")) :
                invoke(stream_prefix_hash_scan_variants!, Tuple{ntuple(_ -> Any, length(args) + 1)...},
                    out, args...; kwargs...)))
            override = only(filter(m -> m.sig.parameters[2] == IOStream, methods(CHOPOFF.stream_prefix_hash_scan_variants!)))
            try
                @test_throws "injected variant scan failure" Base.invokelatest(search_prefixHashScan,guides,
                    joinpath(root,"references","contracts.fa"),out;distance=1,vcf=joinpath(root,"valid","annotations.vcf"))
                @test calls[] == 2
                @test read(out,String) == "preserve existing output"
                @test readdir(dir) == before
            finally
                Base.delete_method(override)
            end
        end
    end
end
