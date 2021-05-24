using Test

using BioSequences
using CRISPRofftargetHunter: safeadd, smallestutype, base_to_idx, 
    getseq, extension, levenshtein, comb_of_d1, comb_of_d, minkmersize, balance
using Combinatorics

@testset "utils.jl" begin

    @testset "safeadd" begin
        @test safeadd(typemax(UInt16), UInt16(1)) == typemax(UInt16)
        @test safeadd(typemax(UInt16) - UInt16(1), UInt16(1)) == typemax(UInt16)
    end

    @testset "smallestutype" begin
        @test smallestutype(UInt(6)) == UInt8
        @test smallestutype(Unsigned(14141311414141414114142)) == UInt128
    end

    @testset "base_to_idx" begin
        @test base_to_idx('A') == 1
        @test_throws String base_to_idx('N')
    end

    @testset "getseq" begin
        @test length(getseq()) == 20
        @test getseq(5, ['A']) == dna"AAAAA"
    end

    @testset "extension" begin
        @test extension("test.bin") == ".bin"
        @test extension("./dad/dada/test.bin") == ".bin"
        @test extension("./dad/dada/test.bin.bin") == ".bin"
    end

    @testset "minkmersize" begin
        @test minkmersize() == 4
        @test minkmersize(10, 10) == 0
        @test minkmersize(5, 1) == 2
    end

    @testset "balance" begin
        @test balance(collect(1:10)) == 5
        @test balance(collect(1:11)) == 6
        @test balance([1, 1, 1, 1, 2, 5, 5, 5, 5]) == 2
    end

    @testset "comb_of_d" begin
        d = 6
        iter = Int(floor(4^6 / 100))
        all_comb = [join(x) for x in multiset_permutations(repeat(['A', 'C', 'T', 'G'], d), d)]
        for i in 1:iter
            seq = rand(all_comb)
            all_comb_dist = [levenshtein(LongDNASeq(seq), LongDNASeq(x), 4) for x in all_comb]
            for dist in [0, 1, 2, 3]
                all_comb_d = all_comb[all_comb_dist .== dist]
                combd, combd_b = comb_of_d(seq, dist)
                @test Set(combd) == Set(all_comb_d)
            end
        end
    end
end