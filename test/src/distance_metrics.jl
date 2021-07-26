using Test

using CRISPRofftargetHunter: isinclusive, commonprefix, hamming, getseq,
    levenshtein, pa_sa, align, prefix_align, suffix_align!, suffix_align
using BioSequences

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

        @test !isinclusive(DNA_Gap, DNA_N)
        @test !isinclusive(DNA_Gap, DNA_A)
    end

    @testset "commonprefix" begin
        @test commonprefix(dna"ACTGACTG", dna"ACTGACTG") == 8
        @test commonprefix(dna"GCTGACTG", dna"ACTGACTG") == 0
        @test commonprefix(dna"NCTGACTG", dna"ACTGACTG") == 8
        @test commonprefix(dna"RCTGACTG", dna"WCTGACTG", isinclusive) == 0

        @test commonprefix(dna"ACTGACTG", dna"ACTGACTGAAAAA") == 8
        @test commonprefix(dna"GCTGACTG", dna"ACTGACTGAAAAA") == 0

        @test commonprefix(dna"NCTGACTG", dna"ACTGACTG", iscompatible) == 8
        @test commonprefix(dna"RCTGACTG", dna"WCTGACTG", iscompatible) == 8
    end

    @testset "hamming" begin
        @test hamming(dna"ACTGACTG", dna"ACTGACTG") == 0
        @test hamming(dna"AGTGACGG", dna"ACTGACTG") == 2
        @test hamming(dna"AGTGACGG", dna"ACTGACTGGGG") == 2
        @test hamming(dna"AGTGACGGGGG", dna"ACTGACTG") == 2

        @test hamming(dna"RGTGACGG", dna"WCTGACTG", isinclusive) == 3
        @test hamming(dna"RGTGACGG", dna"ACTGACTG", isinclusive) == 2

        @test hamming(dna"RGTGACGG", dna"WCTGACTG") == 2
        @test hamming(dna"RGTGACGG", dna"ACTGACTG") == 2
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
                          dna"WCTGAAA", 4, isinclusive) == 1
        @test levenshtein(dna"CCTG",
                          dna"NCTRAAA", 4, isinclusive) == 0

        # dist == k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 4) == 4
        @test levenshtein(dna"TGAGAAAAAAC",
                          dna"GGAGAAAAAAG", 2) == 2

        # dist > k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 2) == 3

        # test skip at the begining
        @test levenshtein(dna"TA", dna"CCTA", 1) == 2
        @test levenshtein(dna"TA", dna"CTA", 1) == 1
        @test levenshtein(dna"TA", dna"CCTA", 2) == 2
        @test levenshtein(dna"TAT", dna"CCTAT", 1) == 2

        # better to take 2 mismatches than 3 gaps and 2 matches!
        # TA
        # CCCTA
        @test levenshtein(dna"TA", dna"CCCTA", 3) == 2
        @test levenshtein(dna"TAT", dna"CCTAT", 1) == 2

        # capped at 2 due to the k < len(guide)
        @test levenshtein(dna"TA", dna"CCTACCC", 10) == 2
        @test levenshtein(dna"TAAAAA", dna"CTAAAAACCCC", 10) == 1

        # A  AGCA
        # AGGAGCA
        @test levenshtein(dna"AAGCA", dna"AGGAGCA", 5) == 2

        # AAGCA
        # AGG AGTT
        @test levenshtein(dna"AAGCA", dna"AGGAGTT", 5) == 2

        #   TGAGAA 
        # CATCAAAAA
        @test levenshtein(dna"TGAGAA", dna"CATCAAAAA", 6) == 4

        # test gaps in the ref
        # AGGACC
        # TGG-CCAA
        levenshtein(dna"AGGACC", dna"TGGCCAA", 4) == 2

        # AGGACCT
        # TGG-CCAA
        levenshtein(dna"AGGACCT", dna"TGGCCAA", 5) == 3

        # ACCC
        # -CCCGGG
        levenshtein(dna"ACCC", dna"CCCGGG", 5) == 1

        # random values tested with results from
        # pairalign(LevenshteinDistance(), s, t)
        # without endings
        @test levenshtein(dna"CATGGTCGTTTGCCAAATGG",
                          dna"GTTTTTTAGGACGTCCAGGTAGTG", 20) == 12
        @test levenshtein(dna"TAAGTGGGTTGATCTTGGAG",
                          dna"AACGACGTATCTGATCTATTCTAT", 20) == 10
        @test levenshtein(dna"TGAACTTGCATCTTTCCCGC",
                          dna"GGCGTGAAGATAAAGGCCCCGATA", 20) == 12
        @test levenshtein(dna"GAGACCAGGAGAGTTATCCC",
                          dna"TTCTATATCCATTCAGACCTGTCT", 20) == 13
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 20) == 9
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 4) == 5
    end

    @testset "prefix_ and suffix_ alignment" begin
        @test levenshtein(dna"ACTG",
                          dna"ACTG", 4) == pa_sa(dna"ACTG", dna"ACTG", 4, 2).dist
        @test levenshtein(dna"ACTG",
                          dna"ACTGAAA", 4) == pa_sa(dna"ACTG", dna"ACTGAAAA", 4, 2).dist
        @test levenshtein(dna"GCTG",
                          dna"ACTGAAA", 4) == pa_sa(dna"GCTG", dna"ACTGAAAA", 4, 4).dist
        @test levenshtein(dna"GCTGAAA",
                          dna"ACTG", 4) == pa_sa(dna"GCTGAAA", dna"ACTG", 4, 4).dist
        @test levenshtein(dna"RCTG",
                          dna"WCTGAAA", 4) == pa_sa(dna"RCTG", dna"WCTGAAA", 4, 2).dist
        @test levenshtein(dna"CCTG",
                          dna"NCTRAAA", 4) == pa_sa(dna"CCTG", dna"NCTRAAA", 4, 2).dist

        # dist == k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 4) == pa_sa(dna"TGAGAA", dna"CATCAAAAA", 4, 3).dist
        @test levenshtein(dna"TGAGAAAAAAC",
                          dna"GGAGAAAAAAG", 2) == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 2, 5).dist
        @test pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 5).dist == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 6).dist
        @test pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 2).dist == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 3).dist
        @test pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 1).dist == pa_sa(dna"TGAGAAAAAAC", dna"GGAGAAAAAAG", 4, 7).dist

        # dist > k
        @test levenshtein(dna"TGAGAA",
                          dna"CATCAAAAA", 2) == pa_sa(dna"TGAGAA", dna"CATCAAAAA", 2, 5).dist

        # random values tested with results from
        # pairalign(LevenshteinDistance(), s, t)
        # without endings
        @test levenshtein(dna"CATGGTCGTTTGCCAAATGG",
                          dna"GTTTTTTAGGACGTCCAGGTAGTG", 20) == pa_sa(
                          dna"CATGGTCGTTTGCCAAATGG",
                          dna"GTTTTTTAGGACGTCCAGGTAGTG", 20, 15).dist
        @test levenshtein(dna"TAAGTGGGTTGATCTTGGAG",
                          dna"AACGACGTATCTGATCTATTCTAT", 20) == pa_sa(
                          dna"TAAGTGGGTTGATCTTGGAG",
                          dna"AACGACGTATCTGATCTATTCTAT", 20, 13).dist
        @test levenshtein(dna"TGAACTTGCATCTTTCCCGC",
                          dna"GGCGTGAAGATAAAGGCCCCGATA", 20) == pa_sa(
                          dna"TGAACTTGCATCTTTCCCGC",
                          dna"GGCGTGAAGATAAAGGCCCCGATA", 20, 12).dist
        @test levenshtein(dna"GAGACCAGGAGAGTTATCCC",
                          dna"TTCTATATCCATTCAGACCTGTCT", 20) == pa_sa(
                          dna"GAGACCAGGAGAGTTATCCC",
                          dna"TTCTATATCCATTCAGACCTGTCT", 20, 15).dist
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 20) == pa_sa(
                          dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 20, 10).dist
        @test levenshtein(dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 4) == pa_sa(
                          dna"CCGTAGCCTGTCCTCCTATA",
                          dna"TCACATGCGCACGTCCTCATATCT", 4, 7).dist
    end


    @testset "all alignments friction test" begin
        iter = 1000000
        k = rand(collect(1:10), iter)
        guide_sizes = rand(collect(1:20), iter)
        prefix_len = rand(collect(1:10), iter)
        for i in 1:iter
            if k[i] > guide_sizes[i]
                k[i] = guide_sizes[i]
            end
            g = getseq(guide_sizes[i])
            ref = getseq(guide_sizes[i] + k[i])
            #@show "$i $g $ref " * string(k[i]) * " " * string(prefix_len[i])
            aln = align(g, ref, k[i])
            if aln.dist <= k[i]
                @test aln.dist == hamming(LongDNASeq(aln.guide), LongDNASeq(aln.ref), isequal)    
            end
            @test levenshtein(g, ref, k[i]) == aln.dist
            @test aln.dist == pa_sa(g, ref, k[i] + prefix_len[i], k[i]).dist
        end
    end


    @testset "prefix_ and suffix_! alignment" begin

        function mutate_suffix!(suffix::LongDNASeq, changes::Int = 3)
            min = length(suffix)
            for i in changes
                idx = Int(ceil(rand() * length(suffix)))
                suffix[idx] = getseq(1)[1]
                if min > idx
                    min = idx
                end
            end
            return min
        end

        guide = getseq() 
        prefix_len = 7
        d = 4
        ref = copy(guide)
        prefix = ref[1:prefix_len]
        suffix = ref[prefix_len+1:end]

        # prefix - ref = 0D 
        # suffix - ref = 0D 
        pa = prefix_align(guide, prefix, length(suffix), d)
        pa_old = prefix_align(guide, prefix, length(suffix), d)
        
        aln = suffix_align(suffix, pa_old)
        aln_ = suffix_align!(suffix, pa) # this will progresively mutate pa
        @test aln == aln_
        @test pa_old != pa

        # mutate suffix to have 2 changes
        change_at = mutate_suffix!(suffix, 2)
        aln = suffix_align(suffix, pa_old)
        aln_ = suffix_align!(suffix, pa, change_at)
        @test aln == aln_
        
    end
end