using Test

using CHOPOFF
using BioSequences
using CSV
using DataFrames

# Fast debug test for sassy Cas9 support
# This is a subset of the db.jl tests that can run quickly for debugging

genome = joinpath(dirname(pathof(CHOPOFF)), "..",
    "test", "sample_data", "genome", "semirandom.fa")
guides_s = Set(readlines(joinpath(dirname(pathof(CHOPOFF)), "..",
    "test", "sample_data", "guides.txt")))
guides = LongDNA{4}.(guides_s)
tdir = tempname()
mkpath(tdir)

@testset "sassy Cas9 debug" begin
    # Set JULIA_NUM_THREADS to 1 for easier debugging
    # Threads.resize_nthreads(1)
    # First, let's just check if sassy is finding anything for the first guide
    first_guide = first(guides)
    println("Testing guide: ", first_guide)
    println("Reversed guide: ", reverse(first_guide))
    println()

    motif = Motif("Cas9")

    # Build linearDB as gold standard
    sassy_ldb_path = joinpath(tdir, "sassy_linearDB")
    mkpath(sassy_ldb_path)
    build_linearDB("samirandom_sassy", genome, motif, sassy_ldb_path, 7)
    sassy_detail_path = joinpath(sassy_ldb_path, "detail.csv")

    for d in 1:3
        sassy_path = joinpath(tdir, "sassy_d$d.csv")
        search_sassy(guides, genome, motif, sassy_path; distance = d)
        sassy_df = DataFrame(CSV.File(sassy_path))

        search_linearDB(sassy_ldb_path, guides, sassy_detail_path; distance = d)
        ldb_df = DataFrame(CSV.File(sassy_detail_path))

        # Compare results - all linearDB results must be in sassy
        failed = antijoin(ldb_df, sassy_df, on = [:guide, :distance, :chromosome, :start, :strand])
        if nrow(failed) > 0
            println("Distance $d: linearDB has $(nrow(failed)) results not in sassy")
            println(first(failed, 5))
            # Check if sassy finds these with a different position
            for row in eachrow(first(failed, 3))
                matching = filter(s -> s.guide == row.guide && s.chromosome == row.chromosome && s.strand == row.strand, sassy_df)
                if nrow(matching) > 0
                    println("  sassy found $(nrow(matching)) matches for this guide/chrom/strand:")
                    println("  ", first(matching, 3))
                end
            end
        end
        @test nrow(failed) == 0

        # And all sassy results must be in linearDB
        failed2 = antijoin(sassy_df, ldb_df, on = [:guide, :distance, :chromosome, :start, :strand])
        if nrow(failed2) > 0
            println("Distance $d: sassy has $(nrow(failed2)) results not in linearDB")
            println(first(failed2, 5))
        end
        @test nrow(failed2) == 0
    end
end
