# julia --project=. -t 4 test/extended/variants.jl
# Optional: CHOPOFF_VCF_SEED, CHOPOFF_VCF_CASES, CHOPOFF_VCF_GT_SEED, CHOPOFF_VCF_GT_CASES.
include(joinpath(@__DIR__,"..","helpers","vcf.jl"))
include(joinpath(@__DIR__,"..","sample_data","vcf","cases.jl"))
include(joinpath(@__DIR__,"..","helpers","vcf_sweeps.jl"))
include(joinpath(@__DIR__,"..","helpers","vcf_sv.jl"))

seed=parse(Int,get(ENV,"CHOPOFF_VCF_SEED","20261009"))
cases=parse(Int,get(ENV,"CHOPOFF_VCF_CASES","200"))
cases>0 || error("CHOPOFF_VCF_CASES must be positive")
rng=MersenneTwister(seed)
println("Extended VCF suite: seed=$seed cases=$cases threads=$(Threads.nthreads())")
flush(stdout)
timing=@timed @testset "Extended VCF [V17] seed=$seed" begin
    @testset "Random case $i" for i in 1:cases
        mktempdir() do dir
            r,edits,guides,motif,d=vtest_random_case(rng,i)
            @test length(vtest_haplotypes(r,edits)) <= 256
            write(joinpath(dir,"reproduce.txt"),"CHOPOFF_VCF_SEED=$seed CHOPOFF_VCF_CASES=$i julia --project=. -t 4 test/extended/variants.jl\n" * repr((r,edits,guides,motif)))
            rows=vtest_compare(dir,r,edits,guides,motif,d;case_id="seed$(seed)_case$i",
                chunk=isodd(i) ? 64 : 97,threads=min(4,Threads.nthreads()),backend=i%3 == 0 ? :portable : :auto)
            @test any(row -> !ismissing(row.variants) && occursin("anchor|",row.variants),eachrow(rows))
        end
    end
    # Indels, MNPs, IUPAC ALTs and multiallelic records packed into one site.
    @testset "Dense mixed case $i" for i in 1:60
        mktempdir() do dir
            name=isodd(i) ? "Cas9" : "Cas12a"
            d=(i-1)%5
            n=name == "Cas9" ? 20 : 21
            g=String(rand(rng,collect("ACGT"),n))
            r=String(rand(rng,collect("ACGT"),40))*(name == "Cas9" ? g*"AGG" : "TTTA"*g)*String(rand(rng,collect("ACGT"),40))
            edits=NamedTuple[]
            while true
                empty!(edits)
                for j in 1:rand(rng,4:7)
                    q=rand(rng,38:(length(r)-38))
                    kind=rand(rng,1:5)
                    id="m$j"
                    if kind == 1
                        push!(edits,vtest_edit(r,q,q+1,vtest_alternative(r[q]),id))
                    elseif kind == 2
                        push!(edits,vtest_edit(r,q,q,String(rand(rng,collect("ACGT"),rand(rng,1:3)-1))*vtest_alternative(r[q-1]),id))
                    elseif kind == 3
                        len=rand(rng,1:3)
                        r[q-1] == r[q+len-1] && continue
                        push!(edits,vtest_edit(r,q,q+len,"",id))
                    elseif kind == 4
                        # IUPAC SNP: every realization differs from the reference.
                        code=r[q] == 'A' ? 'Y' : r[q] == 'C' ? 'R' : r[q] == 'G' ? 'Y' : 'R'
                        append!(edits,[vtest_edit(r,q,q+1,string(b),id;original=string(code),group=id) for b in VCF_TEST_IUPAC[code]])
                    else
                        # Multiallelic SNP plus insertion from one record (shared POS/REF).
                        push!(edits,vtest_edit(r,q,q+1,vtest_alternative(r[q]),id;group=id))
                        push!(edits,vtest_edit(r,q+1,q+1,vtest_alternative(r[q+1]),id;group=id))
                    end
                end
                # One record per group must share POS/REF for vtest_records.
                ok=all(unique(e.group for e in edits)) do grp
                    c=filter(e -> e.group == grp,edits)
                    all(e -> (e.pos,e.ref) == (c[1].pos,c[1].ref),c)
                end
                ok && length(unique(e.pos for e in edits)) == length(unique(e.group for e in edits)) &&
                    length(vtest_haplotypes(r,edits)) <= 1024 && break
            end
            write(joinpath(dir,"reproduce.txt"),"CHOPOFF_VCF_SEED=$seed dense mixed case $i\n"*repr((r,edits,g,name,d)))
            # Every fourth case also checks the ALT cap against the capped oracle.
            k=i%4 == 0 ? rand(rng,1:3) : vtest_exact_cap(edits)
            vtest_compare(dir,r,edits,[LongDNA{4}(g)],Motif(name;distance=d),d;case_id="seed$(seed)_dense$(i)_k$k",
                chunk=isodd(i) ? 2*1024*1024 : 64,threads=min(4,Threads.nthreads()),max_alleles=k)
        end
    end
    @testset "Dense $n-variant combinations" for n in (6,8,10)
        mktempdir() do dir
            r=VCF_TEST_REFERENCE
            edits=[vtest_edit(r,p,p+1,vtest_alternative(r[p]),"dense$p") for p in 41:40+n]
            g=collect(VCF_TEST_GUIDE)
            for i in 1:n
                g[i]=only(edits[i].alt)
            end
            @test length(vtest_haplotypes(r,edits)) == 2^n
            rows=vtest_compare(dir,r,edits,[LongDNA{4}(String(g))],Motif("Cas9";distance=0),0;
                chunk=64,threads=min(4,Threads.nthreads()),case_id="dense$n")
            @test any(row -> row.distance == 0 && length(vtest_tokens(row.variants)) == n,eachrow(rows))
        end
    end
    # Symbolic DEL, DUP and INV with nearby SNPs: oracle-checked literal rows equal symbolic rows.
    @testset "Random symbolic SV case $i" for i in 1:90
        mktempdir() do dir
            r,svs,guides,d=vtest_sv_random_case(rng,i)
            vtest_sv_compare(dir,r,svs,guides,Motif("Cas9";distance=d),d;chunk=isodd(i) ? 64 : 97)
        end
    end
end
# Genotype mode: CHOPOFF_VCF_GT_SEED, CHOPOFF_VCF_GT_CASES (default 200).
haskey(ENV,"CHOPOFF_VCF_GT_CASES") || (ENV["CHOPOFF_VCF_GT_CASES"]="200")
include(joinpath(@__DIR__,"..","variants","genotype_corpus.jl"))
println("Elapsed $(round(timing.time;digits=2)) s; allocated $(timing.bytes) bytes; process peak RSS $(Sys.maxrss()) bytes")
