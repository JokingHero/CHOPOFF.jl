using Test

using CHOPOFF: load, save
using BioSequences

@testset "persistence.jl" begin

    @testset "load and save" begin
        struct TestSeq
            field::Bool
            vec::Vector{LongDNA{4}}
        end
        
        tdir, io = mktemp()
        tseq = TestSeq(true, [dna"ACTG", dna"AAAA"])
        save(tseq, tdir)
        tseq2 = load(tdir)
        @test tseq.field == tseq2.field
        @test tseq.vec == tseq2.vec
    end

    @testset "detail parts are private and published atomically" begin
        with_detail_parts = CHOPOFF.with_detail_parts
        dir = mktempdir()
        output = joinpath(dir, "result.csv")
        # Files that only share the old `detail_` naming must not be touched.
        decoys = Dict(
            joinpath(dir, "detail_user.csv") => "user data\n",
            joinpath(dir, "my_detail_notes.csv") => "notes\n")
        foreach(((path, text),) -> write(path, text), decoys)
        # Parts left behind by a crashed run must not be merged.
        stale = mkpath(joinpath(dir, ".chopoff_parts_stale"))
        write(joinpath(stale, "detail_A.csv"), "stale\n")

        with_detail_parts(output; first_line = "h\n") do parts_dir
            @test dirname(parts_dir) == dir
            write(joinpath(parts_dir, "detail_C.csv"), "c1\nc2")
            write(joinpath(parts_dir, "detail_A.csv"), "a1\n")
        end
        @test read(output, String) == "h\na1\nc1\nc2\n"
        @test all(((path, text),) -> read(path, String) == text, decoys)
        @test read(joinpath(stale, "detail_A.csv"), String) == "stale\n"
        @test count(startswith(".chopoff_parts_"), readdir(dir)) == 1

        # A second run into the same directory sees only its own parts.
        with_detail_parts(output; first_line = "h\n") do parts_dir
            write(joinpath(parts_dir, "detail_A.csv"), "second\n")
        end
        @test read(output, String) == "h\nsecond\n"

        # A failing run keeps the previous output and removes its parts.
        @test_throws ErrorException with_detail_parts(output) do parts_dir
            write(joinpath(parts_dir, "detail_A.csv"), "partial\n")
            error("worker failed")
        end
        @test read(output, String) == "h\nsecond\n"
        @test count(startswith(".chopoff_parts_"), readdir(dir)) == 1
        @test count(startswith("tmp"), readdir(dir)) == 0
    end
end