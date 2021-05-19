using Test

using CRISPRofftargetHunter: load, save
using BioSequences

@testset "persistence.jl" begin

    @testset "load and save" begin
        struct TestSeq
            field::Bool
            vec::Vector{LongSequence{DNAAlphabet{4}}}
        end
        
        tdir, io = mktemp()
        tseq = TestSeq(true, [dna"ACTG", dna"AAAA"])
        save(tseq, tdir)
        tseq2 = load(tdir)
        @test tseq.field == tseq2.field
        @test tseq.vec == tseq2.vec
    end
end