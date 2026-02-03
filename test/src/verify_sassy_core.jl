using Test

using CHOPOFF
using BioSequences
using CSV
using DataFrames

## SET WD when debugging
# cd("test")

@testset "databases" begin
    genome = joinpath(dirname(pathof(CHOPOFF)), "..", 
        "test", "sample_data", "genome", "semirandom.fa")
    guides_s = Set(readlines("./sample_data/guides.txt"))
    guides = LongDNA{4}.(guides_s)
    tdir = tempname()
    mkpath(tdir)
    # guide ACTCAATCATGTTTCCCGTC is on the border - depending on the motif definition
    # it can/can't be found by different methods

    len_noPAM = CHOPOFF.length_noPAM(Motif("Cas9"))

    @testset "linearDB vs sassy" begin
        # Test sassy algorithm matches linearDB for Cas9 (extends5=true)
        motif = Motif("Cas9")
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

            # Compare detail results directly by all columns
            failed = antijoin(ldb_df, sassy_df, on = [:guide, :distance, :chromosome, :start, :strand])
            @test nrow(failed) == 0

            failed2 = antijoin(sassy_df, ldb_df, on = [:guide, :distance, :chromosome, :start, :strand])
            @test nrow(failed2) == 0
        end
    end
end