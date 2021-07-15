using Base: UInt128
using Test

using BioSequences
using CRISPRofftargetHunter: safeadd, smallestutype, base_to_idx, 
    getseq, extension, levenshtein, comb_of_d1, comb_of_d, minkmersize, balance,
    locate_telomeres, findall, expand_ambiguous
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


    @testset "locate_telomeres" begin
        @test locate_telomeres(dna"NACTGN") == (2, 5)
        @test locate_telomeres(dna"ACTGN") == (1, 4)
        @test locate_telomeres(dna"NACTG") == (2, 5)
        @test locate_telomeres(dna"ACTG") == (1, 4)
        @test locate_telomeres(dna"ANGN") == (1, 3)
    end

    @testset "findall" begin
        @test isempty(findall(dna"ACTG", dna"AANN"))
        @test findall(dna"ACTG", dna"NNNG") == [UnitRange(1:4)]
        @test findall(dna"AAANN", dna"ACTGAAAGACTG") == [UnitRange(5:9)]
        @test findall(dna"AAANN", dna"ACTGAAAGA") == [UnitRange(5:9)]

        @test findall(dna"AAANN", dna"ACTGAAAGA", 5) == [UnitRange(5:9)]
        @test findall(dna"AAANN", dna"ACTGAAAGACTG", 5, 9) == [UnitRange(5:9)]
        @test isempty(findall(dna"AAANN", dna"ACTGAAAGACTG", 6, 9))
        @test isempty(findall(dna"AAANN", dna"ACTGAAAGACTG", 2, 8))
    end

    @testset "UInt128 conversion" begin
        x = dna"AAANRAAATGCTACTG"
        y = convert(UInt128, x)
        @test x == convert(LongSequence{DNAAlphabet{4}}, y)
        @test_throws String convert(UInt128, dna"A-A")
        @test_throws String convert(UInt128, dna"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    end

    @testset "expand_ambiguous" begin
        @test expand_ambiguous(dna"ACTG") == [dna"ACTG"]
        @test expand_ambiguous(dna"AR") == [dna"AA", dna"AG"]
        @test expand_ambiguous(dna"WR") == [dna"AA", dna"TA", dna"AG", dna"TG"]
    end
end