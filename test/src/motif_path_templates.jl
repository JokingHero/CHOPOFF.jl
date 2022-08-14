using BioSequences
using Test
using ARTEMIS

@testset "motif_path_templates.jl" begin

    @testset "templates_to_sequences" begin
        motif = Motif("test")
        template = ARTEMIS.build_motifTemplates(motif)

        guide = dna"AAA"
        dist = 1

        pat = ARTEMIS.templates_to_sequences(guide, template; dist = dist)
        @test length(pat) == 8

        ref = ARTEMIS.comb_of_d1_extended_ref("AAAA")
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
    end
end