using BioSequences
using Test
using ARTEMIS

@testset "motif_path_templates.jl" begin

    @testset "templates_to_sequences" begin
        dist = 1
        motif = Motif("test"; distance = dist)
        template = ARTEMIS.build_PathTemplates(length_noPAM(motif), motif.distance)
        guide = dna"AAA"
        
        # might have redundant sequences
        pat = ARTEMIS.templates_to_sequences(guide, template; dist = dist)
        @test length(pat) == 21
        # 
        pat2 = ARTEMIS.templates_to_sequences(guide, template, motif)
        @test length(pat) == 21

        pat


        # expands to the maximal alignment length
        pat3 = ARTEMIS.templates_to_sequences_extended(guide, template)
        pat3_union = Set(string.(collect(union(pat3...))))
        
        # test templates _extended
        ref = ARTEMIS.comb_of_d1_extended_ref("AAA")
        @test setdiff(pat3_union, ref) == Set(["AAAA"])


        ref = map(x -> ARTEMIS.align(guide, LongDNA{4}(x), dist), collect(ref))
        function remove_gap(x::String)
            gaps = map(x -> x == '-', collect(x))
            x = collect(x)[.!gaps]
            return LongDNA{4}(join(x))
        end

        ref = map(x -> remove_gap(x.ref), ref)
        ref = Set(ref)
        pat = Set(map(x -> x.seq, pat))
        @test pat == ref
        =#
    end
end