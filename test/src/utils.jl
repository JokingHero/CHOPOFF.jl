using Base: UInt128
using Test

using BioSequences
using CRISPRofftargetHunter: safeadd, smallestutype, base_to_idx, 
    getseq, extension, levenshtein, comb_of_d1, comb_of_d, minkmersize, balance,
    locate_telomeres, findall, expand_ambiguous, convert, 
    all_kmers, as_bitvector_of_kmers, as_kmers
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

        # ambig limit
        @test length(findall(dna"AAA", dna"ACTGAAANACTG"; ambig_max = 1)) == 3
        @test isempty(findall(dna"AAANN", dna"ACTGAAANACTG"; ambig_max = 0))
    end

    @testset "UInt128 conversion" begin
        x = dna"AAANRAAATGCTACTG"
        y = convert(UInt128, x)
        @test x == convert(LongDNASeq, y)
        x = LongDNASeq("GGAAATGCCCCGCGAACAGG")
        @test convert(LongDNASeq, convert(UInt128, x)) == x
        @test_throws String convert(UInt128, dna"A-A")
        @test_throws String convert(UInt128, dna"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

        for i in 1:10000
            x = getseq(ceil(Int, rand()*24), ['N', 'A', 'C', 'G', 'T'])
            @test convert(LongDNASeq, convert(UInt128, x)) == x
            @test LongDNASeq(convert(UInt128, x), length(x)) == x
        end
    end

    @testset "expand_ambiguous" begin
        @test expand_ambiguous(dna"ACTG") == [dna"ACTG"]
        @test expand_ambiguous(dna"AR") == [dna"AA", dna"AG"]
        @test expand_ambiguous(dna"WR") == [dna"AA", dna"TA", dna"AG", dna"TG"]
    end

    @testset "all_kmers" begin
        for i in 1:6
            @test length(all_kmers(i)) == 4^i
        end
    end

    @testset "as_kmers" begin
        @test first(as_kmers(LongDNASeq(repeat('A', 20)), 3)) == dna"AAA"
        @test isempty(setdiff(
            Set(as_kmers(dna"ACTGR", 4)), 
            Set([dna"ACTG", dna"CTGA", dna"CTGG"])))
        @test isempty(setdiff(
            Set(as_kmers(dna"ACTGR", 3)), 
            Set([dna"ACT", dna"CTG", dna"TGA", dna"TGG"])))
    end

    @testset "as_bitvector_of_kmers" begin
        kmers = all_kmers(2)
        kmers = IdDict(zip(kmers, 1:length(kmers)))
        b = as_bitvector_of_kmers(dna"AAAAAA", kmers)
        @test sum(b) == 1
        @test b[kmers[DNAMer(dna"AA")]]

        kmers = all_kmers(3)
        kmers = IdDict(zip(kmers, 1:length(kmers)))
        b = as_bitvector_of_kmers(dna"ACTG", kmers)
        @test sum(b) == 2
        @test b[kmers[DNAMer(dna"ACT")]]
        @test b[kmers[DNAMer(dna"CTG")]]
    end
end