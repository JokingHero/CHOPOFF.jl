using Test

using ARTEMIS: DBInfo, Loc, decode, 
    Motif, length_noPAM, removepam, combinestrings, notX,
    AmbigIdx, findbits, setdist
using BioSequences

@testset "structures" begin

    cas9 = Motif("Cas9")
    cpf1 = Motif("Cpf1")
    cas9_d1 = setdist(cas9, 1)
    cpf1_d1 = setdist(cpf1, 1)

    cas9_2 = Motif("Cas9", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 1, true, 0)

    @testset "Motif" begin
        @test length_noPAM(cas9) == 20
        @test length_noPAM(cpf1) == 20
        @test removepam(dna"ACTNN", 1:3) == dna"NN"
        @test combinestrings("XXXACT", "ACTXXX") == "ACTACT"
        @test length(cas9) == 23
        @test length(cpf1) == 24
        @test length(cas9_d1) == 23
        @test length(cpf1_d1) == 24
        @test isequal(cas9_d1, cas9_2)
    end


    @testset "DBInfo & Loc" begin
        dbi = DBInfo("./sample_data/genome/semirandom.fa", "test", cas9)
        @test dbi.gi.is_fa == true
        @test length(dbi.gi.chrom) == 8
        loc = Loc{dbi.gi.chrom_type, dbi.gi.pos_type}(1, 10, true)
        @test "semirandom1,10,+" == decode(loc, dbi)
    end


    @testset "AmbigIdx" begin
        guides = [dna"ACTG", dna"NNAC", dna"GGAA", dna"GGAA"]
        annot = ["rs131;rs1", "1", "2", "3"]
        idx = AmbigIdx(guides, annot)
        @test sum(findbits(dna"AAAC", idx)) == 1
        @test sum(findbits(dna"GGAA", idx)) == 2
        @test sum(findbits(dna"GCAA", idx)) == 0
        @test sum(findbits(dna"GGA", idx)) == 3
        @test idx.annot[findbits(dna"ACTG", idx)][1] == "rs131;rs1"
    end
end