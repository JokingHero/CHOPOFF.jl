using Test

using ARTEMIS
using BioSequences

@testset "find_offtargets.jl" begin

    genome = joinpath(vcat(
        splitpath(dirname(pathof(ARTEMIS)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))

    # construct example DBInfo
    dbi = DBInfo(genome, "Cas9_semirandom_noVCF", Motif("Cas9"))
    # finally gather all off-targets
    guides = Vector{String}()
    ambig = gatherofftargets!(guides, dbi)

    guides2 = Vector{UInt64}()
    ambig2 = gatherofftargets!(guides2, dbi)
    @test ambig == ambig2
    guide_with_extension_len = length_noPAM(dbi.motif) + dbi.motif.distance

    guides2 = String.(LongDNA{4}.(guides2, guide_with_extension_len))
    @test guides == guides2
    
    guides3 = Vector{UInt128}()
    ambig3 = gatherofftargets!(guides3, dbi)
    @test ambig3 == ambig2

    guides3 = String.(LongDNA{4}.(guides3, guide_with_extension_len))
    @test guides == guides3
end