using Test

using CRISPRofftargetHunter
using BioSequences
using CSV
using DataFrames

## SET WD when debugging
# cd("test")

## CRISPRitz compare functions - we test with up to 4 distance
function asguide(x::String)
    x = x[1:(length(x) - 3)]
    x = replace.(x, "-" => "")
    @assert length(x) == 20
    return x
end

function countspaces(x)
    return length(findall("-", x))
end

function ldb_start(pos, rrna_len, czdna_spac, strand)
    start = Vector{Int}()
    for i in 1:length(pos)
        if strand[i] == "+"
            i_start = pos[i] + rrna_len[i] - czdna_spac[i]
        else
            i_start = pos[i] + 1
        end
        push!(start, i_start)
    end
    return start
end

function rows_not_in(len::Int, rows_in::Vector{Int})
    all_rows = collect(1:len)
    rows_in_ = sort(unique(rows_in))
    if len == length(rows_in_)
        return []
    end
    deleteat!(all_rows, rows_in)
    return all_rows
end

@testset "databases" begin

    genome = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..", 
        "test", "sample_data", "genome", "semirandom.fa")
    guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
    guides = LongDNASeq.(guides_s)
    tdir = tempname()
    mkpath(tdir)

    # make and run default linearDB
    ldb_path = joinpath(tdir, "linearDB")
    mkpath(ldb_path)
    build_linearDB("samirandom", genome, Motif("Cas9"), ldb_path, 7)
    detail_path = joinpath(ldb_path, "detail.csv")
    ldb_res = search_linearDB(ldb_path, guides, 3; detail = detail_path)
    ldb = DataFrame(CSV.File(detail_path))

    #= @testset "linearDB against CRISPRitz" begin
        ## Files
        cz_file = "./sample_data/crispritz_results/guides.output.targets.txt"
        cz = DataFrame(CSV.File(cz_file))
        cz = cz[cz.Total .<= 3, :]
        cz.guide = asguide.(cz.crRNA)
        cz.start = ldb_start(cz.Position, length.(cz.crRNA), countspaces.(cz.DNA), cz.Direction)

        # list all alignments that are in linearDB and not in crizpritz output
        for g in guides_s
            czg = cz[cz.guide .== g, :]
            ldbg = ldb[ldb.guide .== g, :]

            isfound_in_czg = Vector{Int}() # indexes of rows of ldbg that are found in czg
            isfound_in_ldbg = Vector{Int}() # indexes of rows of czg that are found in ldbg
            for j in 1:nrow(ldbg)
                issame = (czg.Total .== ldbg.distance[j]) .&
                    (czg.Direction .== ldbg.strand[j]) .&
                    (czg.start .== ldbg.start[j]) .&
                    (czg.Chromosome .== ldbg.chromosome[j])
                if any(issame)
                    j_idx = findall(issame)
                    j_idx = j_idx[1]
                    push!(isfound_in_czg, j)
                    isfound_in_ldbg = vcat(isfound_in_ldbg, findall(
                        (czg.start .== czg.start[j_idx]) .&
                        (czg.Chromosome .== czg.Chromosome[j_idx]) .&
                        (czg.Direction .== czg.Direction[j_idx])))
                end
            end

            ldbg = ldbg[rows_not_in(nrow(ldbg), isfound_in_czg), :]
            czg = czg[rows_not_in(nrow(czg), isfound_in_ldbg), :]

            # crispritz can't handle anything too close to the telomeres NNN
            ldbg = ldbg[ldbg.start .> 125, :]

            # test that all guides are found in crispritz
            @test nrow(ldbg) <= 0
            if nrow(ldbg) > 0
                @info "Failed ldbg finding $g"
                @info "$ldbg"
            end

            # test that all guides are found in linear database
            @test nrow(czg) <= 0
            if nrow(czg) > 0
                @info "Failed cz finding $g"
                @info "$czg"
            end
        end
    end =#

    @testset "linearDB vs sketchDB" begin
        # make and run default sketchDB
        sdb_path = joinpath(tdir, "sketchDB")
        mkpath(sdb_path)
        build_sketchDB("samirandom", genome, Motif("Cas9"), sdb_path)

        sdb_res = search_sketchDB(sdb_path, guides, 2)
        @test nrow(sdb_res) == length(guides)
        @test all(sdb_res.guide .== guides)
        @test all(ldb_res.guide .== guides)
        ldb_res = convert(Array, ldb_res[:, 1:3])
        sdb_res = convert(Array, sdb_res[:, 1:3])
        for i in 1:length(guides)
            compare = ldb_res[i, :] .<= sdb_res[i, :]
            @test all(compare)
            if !all(compare)
                @info "Failed at guide $i " * string(guides[i])
                @info "linearDB result: " * string(ldb_res[i, :])
                @info "sketchDB result: " * string(sdb_res[i, :])
            end
        end
    end
end