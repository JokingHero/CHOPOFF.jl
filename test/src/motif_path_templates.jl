using BioSequences
using Test
using ARTEMIS
using Combinatorics

@testset "motif_path_templates.jl" begin

    @testset "as UInt64" begin
        for i in 1:100
            guide = ARTEMIS.getseq(20)
            guide_uint8 = ARTEMIS.guide_to_template_format(guide; alphabet = ARTEMIS.ALPHABET_TWOBIT)[1:20] # first 20 are guide
            @test ARTEMIS.asUInt64(guide_uint8) == convert(UInt64, guide)
        end
    end

    @testset "templates_to_sequences" begin
        dist = 2
        motif = Motif("test"; distance = dist)
        g_len = length_noPAM(motif)
        template = ARTEMIS.build_PathTemplates(motif)
        template_minus1 = ARTEMIS.restrictDistance(template, dist - 1)
        @test_throws String ARTEMIS.restrictDistance(template, -1)
        @test_throws String ARTEMIS.restrictDistance(template, 4)
        @test maximum(template_minus1.distances) == (dist - 1)
        @test size(template_minus1.paths)[2] == dist + length_noPAM(motif)

        # all possible combinations for guide + extension (dist = 1) with 3 letters
        all_comb = [join(x) for x in multiset_permutations(repeat(['A', 'C', 'T', 'G'], g_len + dist), g_len + dist)]
        all_comb = LongDNA{4}.(all_comb)
        iter = Int(floor(4^6 / 100)) # how many iterations of guides to test

        # y is the one with longer seqeunces
        function setdiff_by_size(x::Vector{LongDNA{4}}, y::Vector{LongDNA{4}})
            x = copy(x)
            y = copy(y)
            # all of x have to be inside y and all of y have to be inside x
            x_in = falses(length(x))
            y_in = falses(length(y))
            for (i, x_i) in enumerate(x)
                for (j, y_j) in enumerate(y)
                    if y_j[1:length(x_i)] == x_i
                        x_in[i] = true
                        y_in[j] = true
                    end
                end
            end
            return all(x_in) & all(y_in)
        end

        for i in 1:iter
            guide = ARTEMIS.getseq(g_len)
            guide_uint8 = ARTEMIS.guide_to_template_format(guide; alphabet = ARTEMIS.ALPHABET_UINT8)
            all_potential_ot_no_PAM = guide_uint8[template.paths]
            all_potential_ot_no_PAM = reinterpret(DNA, all_potential_ot_no_PAM) # back to DNA
            all_potential_ot_no_PAM = map(LongDNA{4}, eachrow(all_potential_ot_no_PAM))
            
            # do the actual alignment between all possible references and a guide
            all_comb_dist = [levenshtein(guide, x, dist) for x in all_comb]
            for d in 0:dist
                all_comb_d = all_comb[all_comb_dist .== d]
                # the simplest case uses expanded reference, is therefore the same as our all possible cominations
                @test setdiff(Set(all_comb_d), all_potential_ot_no_PAM[template.distances .== d]) == Set([])
            end
        end
    end
end