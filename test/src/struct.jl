using Test

using ARTEMIS: DBInfo, Loc, decode, 
    Motif, length_noPAM, removepam, combinestrings, notX,
    AmbigIdx, findbits, setdist, Offtarget, insert_offtarget!, display_motif
using BioSequences

@testset "structures" begin

    cas9 = Motif("Cas9")
    cpf1 = Motif("Cas12a")
    cas9_d1 = setdist(cas9, 1)
    cpf1_d1 = setdist(cpf1, 1)

    cas9_2 = Motif("Cas9", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 1, true, 0)

    @testset "Motif" begin
        @test length_noPAM(cas9) == 20
        @test length_noPAM(cpf1) == 21
        @test removepam(dna"ACTNN", 1:3) == dna"NN"
        @test combinestrings("XXXACT", "ACTXXX") == "ACTACT"
        @test length(cas9) == 23
        @test length(cpf1) == 25
        @test length(cas9_d1) == 23
        @test length(cpf1_d1) == 25
        @test isequal(cas9_d1, cas9_2)

        @test isnothing(Base.show(cas9))
        @test isnothing(Base.print(cas9))
        @test display_motif(cas9) == "Alias: Cas9\nMaximum search distance: 3\nNumber of allowed ambigous bp: 0\n20N-NGG"
    end


    @testset "DBInfo & Loc" begin
        cas9 = Motif("Cas9")
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

    function rand_offtarget()
        Offtarget(Loc(rand(UInt8.([1, 2, 3])), rand(UInt8.(1:10)), rand([true, false])), rand(1:10), "smth", "smth")
    end

    @testset "Offtarget" begin
        x = Vector{Offtarget}()
        r = 3 # distance of 1 gives range of 3

        # testing case with 0 elements
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(1), true), 2, "", ""), r)
        @test length(x) == 1 && isnothing(o) && i == 2

        # testing case with idx 1
        # in range, but less distance
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(1), true), 0, "", ""), r)
        @test length(x) == 1 && o == 2 && i == 0
        # in range, but now back larger distance
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(0), true), 2, "", ""), r)
        @test length(x) == 1 && isnothing(o) && isnothing(i)

        # testing cases with last idx
        # insert first
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(2), UInt8(1), true), 2, "", ""), r)
        @test length(x) == 2 && isnothing(o) && i == 2
        # in range, but less distance
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(2), UInt8(2), true), 0, "", ""), r)
        @test length(x) == 2 && o == 2 && i == 0
        # in range, but now back larger distance
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(2), UInt8(3), true), 2, "", ""), r)
        @test length(x) == 2 && isnothing(o) && isnothing(i)

        # testing cases with middle idx
        # test the case of insertion in the middle, outside of range
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(10), true), 2, "", ""), r)
        @test length(x) == 3 && isnothing(o) && i == 2
        # in range to middle element, but less distance
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(7), true), 1, "", ""), r)
        @test length(x) == 3 && o == 2 && i == 1
        # in range to middle element, but less distance
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(6), true), 0, "", ""), r)
        @test length(x) == 3 && o == 1 && i == 0
        # in range, but now back larger distance
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(5), true), 2, "", ""), r)
        @test length(x) == 3 && isnothing(o) && isnothing(i)
        (o, i) = insert_offtarget!(x, Offtarget(Loc(UInt8(1), UInt8(7), true), 2, "", ""), r)
        @test length(x) == 3 && isnothing(o) && isnothing(i)
    end
end