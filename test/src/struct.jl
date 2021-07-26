using Test

using CRISPRofftargetHunter: DBInfo, Loc, decode, 
    Motif, length_noPAM, removepam, combinestrings, notX
using BioSequences

@testset "structures" begin

    cas9 = Motif("Cas9")
    cpf1 = Motif("Cpf1")
    @testset "Motif" begin
        @test length_noPAM(cas9) == 20
        @test length_noPAM(cpf1) == 24
        @test removepam(dna"ACTNN", 1:3) == dna"NN"
        @test combinestrings("XXXACT", "ACTXXX") == "ACTACT"
    end


    @testset "DBInfo & Loc" begin
        dbi = DBInfo("./sample_data/genome/semirandom.fa", "test", cas9)
        @test dbi.is_fa == true
        @test length(dbi.chrom) == 8
        loc = Loc{dbi.chrom_type, dbi.pos_type}(1, 10, true)
        @test "semirandom1,10,+" == decode(loc, dbi)
    end
end