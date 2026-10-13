@testset "VCF input and annotations [V01-V03,V06,V10,V14]" begin
    root = joinpath(@__DIR__,"..","sample_data","vcf")
    guide = LongDNA{4}(VCF_TEST_GUIDE)
    genome = joinpath(root,"references","contracts.fa")
    mktempdir() do dir
        out = joinpath(dir,"result.csv")
        run(vcf; kwargs...) = search_prefixHashScan([guide],genome,out;
            distance=1,vcf,early_stopping=fill(typemax(Int),2),kwargs...)
        @testset "$(case.id) preserves existing output" for case in vfixture_invalid()
            write(out,"preserve existing output")
            caught = try
                run(joinpath(root,"invalid",case.id*".vcf"))
                nothing
            catch err
                err
            end
            @test caught isa ErrorException
            @test occursin(case.message,sprint(showerror,caught))
            @test read(out,String) == "preserve existing output"
        end
        @testset "Input errors suggest a bcftools fix" begin
            message(vcf) = try run(vcf); "" catch err; sprint(showerror,err) end
            @test occursin("`bcftools view -t ^absent`",message(joinpath(root,"invalid","unknown_contig.vcf")))
            @test occursin("bcftools norm -f",message(joinpath(root,"invalid","ref_mismatch.vcf")))
            vcf=vtest_vcf(joinpath(dir,"renamed.vcf"),[("1",45,".","G","A",".",".",".")])
            @test occursin("names it chr1, rename with `bcftools annotate --rename-chrs`",message(vcf))
        end
        @testset "Header-only VCF and selected chromosomes" begin
            vcf=vtest_vcf(joinpath(dir,"empty.vcf"),[])
            run(vcf;region="chr2")
            rows=DataFrame(CSV.File(out))
            @test rows.chromosome == ["chr2"]
            @test rows.distance == [0]
            @test all(ismissing,rows.variants)
            @test length(names(rows)) == 12
            run(vcf;region=" chr1, chr2 ")
            @test Set(DataFrame(CSV.File(out)).chromosome) == Set(["chr1","chr2"])
            @test_throws ErrorException run(vcf;region="chr1:1-20")
            @test_throws ErrorException run(vcf;region="targets.bed")
            @test_throws ErrorException run(vcf;region="absent")
            @test_throws ErrorException run(vcf;output=:counts)
            @test_throws ErrorException run(vcf;vcf_annotate=["INFO/absent"])
            @test_throws ErrorException run(vcf;vcf_annotate=["AF"])
            vtest_vcf(vcf,[("chr2",45,"outside","T","Z",".",".",".")])
            run(vcf;region="chr1")
            @test DataFrame(CSV.File(out)).chromosome == ["chr1"]
            @test_throws ErrorException run(vcf;region="chr2")
            @test_throws ErrorException run(joinpath(root,"invalid","unknown_contig.vcf");region="chr1")
            # Reference rows match a search without VCF; region also works without VCF.
            base=joinpath(dir,"base.csv")
            search_prefixHashScan([guide],genome,base;distance=1,region="chr1")
            @test select(DataFrame(CSV.File(out)),names(DataFrame(CSV.File(base)))) == DataFrame(CSV.File(base))
            CHOPOFF.search_prefixHashScan([guide],genome,Motif("Cas9";distance=1),base;
                distance=1,region="chr2",scan_backend=:legacy)
            @test DataFrame(CSV.File(base)).chromosome == ["chr2"]
            search_prefixHashScan([guide],genome,base;distance=1,region="chr1",output=:counts)
            @test DataFrame(CSV.File(base)).D0 == [1]
            @test_throws ErrorException CHOPOFF.search_prefixHashScan([guide],genome,Motif("Cas9";distance=1),out;
                distance=1,vcf,scan_backend=:legacy)
        end
        @testset "Ignored markers preserve literal siblings and emit one warning" begin
            stats=PrefixHashScanStats()
            @test_logs (:warn,r"VCF non-sequence alleles ignored") CHOPOFF.search_prefixHashScan(
                [guide],genome,Motif("Cas9";distance=1),out;distance=1,
                vcf=joinpath(root,"valid","ignored_markers.vcf"),stats)
            @test stats.variants_read == 4
            @test stats.variants_skipped == 4
            rows=DataFrame(CSV.File(out))
            @test collect(skipmissing(rows.variants)) == ["mixed|chr1:47:C>A"]
            @test rows.distance[.!ismissing.(rows.variants)] == [1]
        end
        @testset "Allele cardinalities, sample labels and escaping" begin
            fields=["INFO/AF","INFO/AD","INFO/PAIR","INFO/FLAG","INFO/NOTE","INFO/MISSING",
                "FORMAT/GT","FORMAT/AF","FORMAT/AD","FORMAT/PL"]
            vcf=joinpath(root,"valid","annotations.vcf")
            run(vcf;vcf_annotate=fields)
            rows=DataFrame(CSV.File(out))
            expected=String[]
            for (alt,af,depth,fraction,reads) in (("A","0.1","10","0.25","2"),("C","0.2","20","0.75","6"))
                push!(expected,"id%25q%3Bid2|chr1:45:G>$alt|INFO/AF=$af|INFO/AD=100,$depth|INFO/PAIR=7,8|INFO/FLAG=|INFO/NOTE=a,b\"c%25%7C|INFO/MISSING=.|FORMAT/GT=sample1=0%7C1,sample%252=./.|FORMAT/AF=sample1=$fraction,sample%252=.|FORMAT/AD=sample1=8,$reads,sample%252=.|FORMAT/PL=sample1=0,1,2,3,4,5,sample%252=.")
            end
            @test Set(skipmissing(rows.variants)) == Set(expected)
            @test rows.overlapping_variants[1] == "id%25q%3Bid2|chr1:45:G>A;id%25q%3Bid2|chr1:45:G>C"
            original=read(out)
            @test occursin("\"\"c",String(copy(original)))
            for kind in (:gzip,:bgzf)
                compressed=joinpath(dir,"$kind.vcf.gz")
                stream=kind == :gzip ? GzipCompressorStream(open(compressed,"w")) :
                    CHOPOFF.VariantCallFormat.BCF.BGZFStreams.BGZFStream(compressed,"w")
                try
                    payload=read(vcf)
                    write(stream,payload[1:(end-10)])
                    flush(stream)
                    write(stream,payload[(end-9):end])
                finally
                    close(stream)
                end
                run(compressed;vcf_annotate=fields)
                @test read(out) == original
            end
            # GT, phase, FILTER and frequency metadata do not select alleles.
            run(vcf)
            baseline=read(out)
            text=read(vcf,String)
            for gt in ("1/1","./1","1","0|1|2|0|1|2")
                changed=joinpath(dir,"metadata.vcf")
                write(changed,replace(text,"0|1"=>gt,"PASS"=>"LowQual","AF=0.1,0.2"=>"AF=0,0"))
                run(changed)
                @test read(out) == baseline
            end
        end
        @testset "REF ambiguity accepts supersets only" begin
            vcf=vtest_vcf(joinpath(dir,"ref_iupac.vcf"),[("chr1",45,"refR","R","A",".",".",".")])
            run(vcf)
            @test "refR|chr1:45:R>A" in skipmissing(DataFrame(CSV.File(out)).variants)
            vtest_vcf(vcf,[("chr1",45,"refY","Y","A",".",".",".")])
            @test_throws ErrorException run(vcf)
        end
    end
end
