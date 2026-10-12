include(joinpath(@__DIR__,"..","helpers","vcf_sweeps.jl"))

@testset "VCF reach and ambiguity [V06-V09]" begin
    @testset "SNP spacing $name d=$distance reverse=$reverse_strand gap=$spacing" for
        name in ("Cas9","Cas12a"), distance in 0:4, reverse_strand in (false,true),
        spacing in 1:(length(Motif(name;distance))+distance+1)
        mktempdir() do dir
            reference,edits,guides,motif=vtest_spacing_case(name,distance,spacing,reverse_strand)
            rows=vtest_compare(dir,reference,edits,guides,motif,distance;
                case_id="spacing_$(name)_$(distance)_$(reverse_strand)_$spacing",chunk=64)
            @test any(row -> row.distance == 0 && length(vtest_tokens(row.variants)) == 1 &&
                startswith(only(vtest_tokens(row.variants)),"a|"),eachrow(rows))
        end
    end
    @testset "Maximum consumed span $name d=$distance reverse=$reverse_strand" for
        name in ("Cas9","Cas12a"), distance in 0:4, reverse_strand in (false,true)
        mktempdir() do dir
            g=name == "Cas9" ? VCF_TEST_GUIDE : VCF_TEST_GUIDE*"A"
            expanded=g[1:10]*"C"^distance*g[11:end]
            site=name == "Cas9" ? expanded*"AGG" : "TTTA"*expanded
            r="T"^40*site*"T"^40
            p=name == "Cas9" ? 41 : 40+length(site)
            if reverse_strand
                p=length(r)-p+1
                r=vtest_rc(r)
            end
            edits=[vtest_edit(r,p,p+1,r[p:p],"edge")]
            rows=vtest_compare(dir,r,edits,[LongDNA{4}(g)],Motif(name;distance),distance;chunk=64,case_id="maxspan")
            @test any(row -> !ismissing(row.variants) && abs(row.start-row.end)+1 == length(site),eachrow(rows))
        end
    end
    @testset "Reference ambiguity budget $budget with $ambiguous bases" for budget in 0:3, ambiguous in 0:3
        mktempdir() do dir
            r=VCF_TEST_REFERENCE
            chars=collect(r)
            chars[42:41+ambiguous].='N'
            r=String(chars)
            edits=[vtest_edit(r,45,46,"A","alt")]
            rows=vtest_compare(dir,r,edits,[LongDNA{4}(VCF_TEST_GUIDE)],Motif("Cas9";distance=1,ambig_max=budget),1;chunk=64)
            @test any(!ismissing,rows.variants) == (ambiguous <= budget)
        end
    end
    @testset "Ambiguous guide and literal ALT" begin
        mktempdir() do dir
            g=collect(VCF_TEST_GUIDE);g[5]='R'
            rows=vtest_compare(dir,VCF_TEST_REFERENCE,[vtest_edit(VCF_TEST_REFERENCE,45,46,"A","alt")],
                [LongDNA{4}(String(g))],Motif("Cas9";distance=0),0;chunk=64)
            @test any(row -> row.distance == 0 && !ismissing(row.variants),eachrow(rows))
        end
    end
end
