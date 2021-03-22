using CRISPRofftargetHunter
using BioSequences
using Test

@testset "distance_metrics.jl" begin

    @testset "isinclusive" begin
        @test isinclusive(DNA_A, DNA_A)
        @test isinclusive(DNA_N, DNA_A)
        @test isinclusive(DNA_N, DNA_W)
        @test isinclusive(DNA_A, DNA_N)
        @test isinclusive(DNA_W, DNA_N)
        @test isinclusive(DNA_V, DNA_S)

        @test !isinclusive(DNA_W, DNA_R)
        @test !isinclusive(DNA_R, DNA_W)
        @test !isinclusive(DNA_A, DNA_C)
    end

    @testset "commonprefix" begin
        @test commonprefix(dna"ACTGACTG", dna"ACTGACTG") == 8
        @test commonprefix(dna"GCTGACTG", dna"ACTGACTG") == 0
        @test commonprefix(dna"NCTGACTG", dna"ACTGACTG") == 8
        @test commonprefix(dna"RCTGACTG", dna"WCTGACTG") == 0

        @test commonprefix(dna"NCTGACTG", dna"ACTGACTG", iscompatible) == 8
        @test commonprefix(dna"RCTGACTG", dna"WCTGACTG", iscompatible) == 8
    end

    @testset "hamming" begin
        @test hamming(dna"ACTGACTG", dna"ACTGACTG") == 0
        @test hamming(dna"AGTGACGG", dna"ACTGACTG") == 2
        @test hamming(dna"AGTGACGG", dna"ACTGACTGGGG") == 2
        @test hamming(dna"AGTGACGGGGG", dna"ACTGACTG") == 2

        @test hamming(dna"RGTGACGG", dna"WCTGACTG") == 3
        @test hamming(dna"RGTGACGG", dna"ACTGACTG") == 2

        @test hamming(dna"RGTGACGG", dna"WCTGACTG", iscompatible) == 2
        @test hamming(dna"RGTGACGG", dna"ACTGACTG", iscompatible) == 2
    end

    @testset "levenshtein" begin
        # default k = 4
        # dist <k
        @test levenshtein(dna"ACTG",
                          dna"ACTG", 4) == 0
        @test levenshtein(dna"ACTG",
                          dna"ACTGAAA", 4) == 0
        @test levenshtein(dna"GCTG",
                          dna"ACTGAAA", 4) == 1
        @test levenshtein(dna"GCTGAAA",
                          dna"ACTG", 4) == 4
        @test levenshtein(dna"RCTG",
                          dna"WCTGAAA", 4) == 1
        @test levenshtein(dna"CCTG",
                          dna"NCTRAAA", 4) == 0

        # dist == k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 4) == 4
        @test levenshtein(dna"TGAGAAAAAAC",
                          dna"GGAGAAAAAAG", 2) == 2

        # dist > k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 2) == 3

        # random values tested with results from
        # pairalign(LevenshteinDistance(), s, t)
        # without endings
        @test levenshtein(dna"CATGGTCGTTTGCCAAATGG",
                          dna"GTTTTTTAGGACGTCCAGGTAGTG", 20) == 13
        @test levenshtein(dna"TAAGTGGGTTGATCTTGGAG",
                          dna"AACGACGTATCTGATCTATTCTAT", 20) == 12
        @test levenshtein(dna"TGAACTTGCATCTTTCCCGC",
                          dna"GGCGTGAAGATAAAGGCCCCGATA", 20) == 14
        @test levenshtein(dna"GAGACCAGGAGAGTTATCCC",
                          dna"TTCTATATCCATTCAGACCTGTCT", 20) == 14
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 20) == 9
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 4) == 5
    end

    @testset "prefix_ and suffix_ levenshtein" begin
        # default k = 4
        # dist <k
        function pa_sa(guide::LongDNASeq, ref::LongDNASeq, d::Int, prefixlen::Int)
            prefix = ref[1:prefixlen]
            suffix = ref[prefixlen+1:end]
            pa = prefix_levenshtein(guide, prefix, d)
            return pa.isfinal ? pa.dist : suffix_levenshtein(guide, suffix, pa, d)
        end

        @test levenshtein(dna"ACTG",
                          dna"ACTG", 4) == pa_sa(dna"ACTG", dna"ACTG", 4, 2)
        @test levenshtein(dna"ACTG",
                          dna"ACTGAAA", 4) == pa_sa(dna"ACTG", dna"ACTGAAAA", 4, 2)
        @test levenshtein(dna"GCTG",
                          dna"ACTGAAA", 4) == pa_sa(dna"GCTG", dna"ACTGAAAA", 4, 4)
        @test levenshtein(dna"GCTGAAA",
                          dna"ACTG", 4) == pa_sa(dna"GCTGAAA", dna"ACTG", 4, 2)
        @test levenshtein(dna"RCTG",
                          dna"WCTGAAA", 4) == pa_sa(dna"RCTG", dna"WCTGAAA", 4, 2)
        @test levenshtein(dna"CCTG",
                          dna"NCTRAAA", 4) == pa_sa(dna"CCTG", dna"NCTRAAA", 4, 2)

        # dist == k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 4) == pa_sa(dna"TGAGAA", dna"CATCAAAAA", 4, 3)
        @test levenshtein(dna"TGAGAAAAAAC",
                          dna"GGAGAAAAAAG", 2) == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 2, 5)
        @test pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 5) == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 6)
        @test pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 2) == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 3)
        @test pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 1) == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 7)

        # dist > k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 2) == pa_sa(dna"TGAGAA", dna"CATCAAAAA", 2, 5)

        # random values tested with results from
        # pairalign(LevenshteinDistance(), s, t)
        # without endings
        @test levenshtein(dna"CATGGTCGTTTGCCAAATGG",
                          dna"GTTTTTTAGGACGTCCAGGTAGTG", 20) == pa_sa(
                          dna"CATGGTCGTTTGCCAAATGG",
                          dna"GTTTTTTAGGACGTCCAGGTAGTG", 20, 7)
        @test levenshtein(dna"TAAGTGGGTTGATCTTGGAG",
                          dna"AACGACGTATCTGATCTATTCTAT", 20) == pa_sa(
                          dna"TAAGTGGGTTGATCTTGGAG",
                          dna"AACGACGTATCTGATCTATTCTAT", 20, 7)
        @test levenshtein(dna"TGAACTTGCATCTTTCCCGC",
                          dna"GGCGTGAAGATAAAGGCCCCGATA", 20) == pa_sa(
                          dna"TGAACTTGCATCTTTCCCGC",
                          dna"GGCGTGAAGATAAAGGCCCCGATA", 20, 7)
        @test levenshtein(dna"GAGACCAGGAGAGTTATCCC",
                          dna"TTCTATATCCATTCAGACCTGTCT", 20) == pa_sa(
                          dna"GAGACCAGGAGAGTTATCCC",
                          dna"TTCTATATCCATTCAGACCTGTCT", 20, 7)
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 20) == pa_sa(
                          dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 20, 7)
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 4) == pa_sa(
                          dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 4, 7)
    end
end


@testset "utils.jl" begin

    @testset "deleterange" begin
        @test deleterange(parse(UInt64, "01111001"; base=2), 2, 5) == parse(UInt64, "0111"; base=2)
        @test deleterange(parse(UInt64, "01111001"; base=2), 2, 2) == parse(UInt64, "0111101"; base=2)
        @test deleterange(parse(UInt64, "1001"; base=2), 2, 3) == parse(UInt64, "11"; base=2)
        @test deleterange(parse(UInt64, "01111000"; base=2), 1, 1) == parse(UInt64, "0111100"; base=2)
        @test deleterange(parse(UInt64, "1111111111111111111111111111111111111111111111111111111111111111"; base=2), 3, 64) == parse(UInt64, "011"; base=2)
        @test deleterange(parse(UInt64, "1111111111111111111111111111111111111111111111111111111111111111"; base=2), 1, 64) == parse(UInt64, "0"; base=2)
        @test deleterange(parse(UInt64, "01111000"; base=2), 64, 64) == parse(UInt64, "01111000"; base=2)

        @test deleterange(parse(UInt128, "01111001"; base=2), 2, 5) == parse(UInt128, "0111"; base=2)
        @test deleterange(parse(UInt128, "01111001"; base=2), 2, 2) == parse(UInt128, "0111101"; base=2)
        @test deleterange(parse(UInt128, "1001"; base=2), 2, 3) == parse(UInt128, "11"; base=2)
        @test deleterange(parse(UInt128, "01111000"; base=2), 1, 1) == parse(UInt128, "0111100"; base=2)
        @test deleterange(parse(UInt128, "1111111111111111111111111111111111111111111111111111111111111111"; base=2), 3, 128) == parse(UInt128, "011"; base=2)
        @test deleterange(parse(UInt128, "1111111111111111111111111111111111111111111111111111111111111111"; base=2), 1, 128) == parse(UInt128, "0"; base=2)
        @test deleterange(parse(UInt128, "01111000"; base=2), 128, 128) == parse(UInt128, "01111000"; base=2)
    end
end
