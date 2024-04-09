using Test

using CHOPOFF
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
    for i in 1:lastindex(pos)
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
    deleteat!(all_rows, rows_in_)
    return all_rows
end


function compare_result(res::DataFrame, res2::DataFrame; less_or_equal::Bool = false)
    if (nrow(res) != nrow(res2)) 
        throw("Unequal row count to compare.") 
    end
    res.guide = LongDNA{4}.(res.guide)
    res2.guide = LongDNA{4}.(res2.guide)
    nres = propertynames(res)
    nres2 = propertynames(res2)
    if length(nres) > length(nres2)
        res = select(res, nres2)
        nres = nres2
    else
        res2 = select(res2, nres)
    end
    comb = outerjoin(res, res2, on = [:guide], makeunique = true)
    nres_noguide = filter(x -> x != :guide, nres)
    for col in nres_noguide
        col1 = Symbol(string(col) * ("_1"))
        if less_or_equal
            if !all(comb[:, col] .<= comb[:, col1])
                return false
            end
        else
            if !all(comb[:, col] .== comb[:, col1])
                return false
            end
        end
    end
    return true
end


@testset "databases" begin
    genome = joinpath(dirname(pathof(CHOPOFF)), "..", 
        "test", "sample_data", "genome", "semirandom.fa")
    guides_s = Set(readlines("./sample_data/guides.txt"))
    guides = LongDNA{4}.(guides_s)
    tdir = tempname()
    mkpath(tdir)
    # guide ACTCAATCATGTTTCCCGTC is on the border - depending on the motif definition
    # it can/can't be found by different methods

    @testset "vcfDB result is same as in saved file" begin
        vcf = joinpath(dirname(pathof(CHOPOFF)), "..", 
        "test", "sample_data", "artificial.vcf")
        vcf_path = joinpath(tdir, "vcfDB")
        mkpath(vcf_path)
        vcf_storage_path = joinpath(vcf_path, "vcfDB.bin")
        build_vcfDB(
            "samirandom", genome, vcf,
            Motif("Cas9"; distance = 2, ambig_max = 3), vcf_storage_path)
        
        detail_path_vcf = joinpath(vcf_path, "output.csv")
        search_vcfDB(vcf_storage_path, guides, detail_path_vcf; distance = 2, 
            early_stopping = [300, 300, 300])
        vcf_detail = DataFrame(CSV.File(detail_path_vcf))
        vcf_res = summarize_offtargets(vcf_detail; distance = 2)

        ar_file = joinpath(dirname(pathof(CHOPOFF)), "..", 
            "test", "sample_data", "artificial_results.csv")
        ar_detail = DataFrame(CSV.File(ar_file))
        ar = summarize_offtargets(ar_detail; distance = 2)
        @test compare_result(ar, vcf_res)
    end

    # make and run default linearDB
    ldb_path = joinpath(tdir, "linearDB")
    mkpath(ldb_path)
    build_linearDB("samirandom", genome, Motif("Cas9"), ldb_path, 7)
    detail_path = joinpath(ldb_path, "detail.csv")
    search_linearDB(ldb_path, guides, detail_path; distance = 3)
    ldb = DataFrame(CSV.File(detail_path))
    ldb_res = summarize_offtargets(ldb; distance = 3)

    # make and run default dictDB
    dictDB = build_dictDB(
        "samirandom", genome, 
        Motif("Cas9"; distance = 2))
    ddb_res = search_dictDB(dictDB, guides)

    # make and run default hashDB
    hashDB = build_hashDB(
        "samirandom", genome, 
        Motif("Cas9"; distance = 1, ambig_max = 0))
    hdb_res = search_hashDB(hashDB, guides, false)

    len_noPAM = CHOPOFF.length_noPAM(Motif("Cas9"))

    @testset "linearDB against CRISPRitz" begin
        ## Files
        cz_file = "./sample_data/crispritz_results/guides.output.targets.txt"
        cz = DataFrame(CSV.File(cz_file))
        cz = cz[cz.Total .<= 3, :]
        cz.guide = asguide.(String.(cz.crRNA))
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
    end


    @testset "linearDB vs dictDB" begin
        @test compare_result(ldb_res, ddb_res)
    end


    @testset "linearDB vs hashDB" begin
        @test compare_result(ldb_res, hdb_res)
    end

    
    @testset "hashDB vs dictDB" begin
        @test compare_result(ddb_res, hdb_res)
    end

    
    @testset "linearDB vs treeDB" begin
        tdb_path = joinpath(tdir, "treeDB")
        mkpath(tdb_path)
        build_treeDB("samirandom", genome, Motif("Cas9"), tdb_path, 7)
        detail_path = joinpath(tdb_path, "detail.csv")

        # this should work without errors
        inspect_treeDB(tdb_path; inspect_prefix = "CCGTCGC")

        for d in 1:3
            search_treeDB(tdb_path, guides, detail_path; distance = d)
            tdb = DataFrame(CSV.File(detail_path))
            tdb_res = summarize_offtargets(tdb; distance = d)
            @test compare_result(ldb_res, tdb_res)
        end

        # for final distance check also detail output
        tdb = DataFrame(CSV.File(detail_path))
        failed = antijoin(ldb, tdb, on = [:guide, :distance, :chromosome, :start, :strand])
        @test nrow(failed) == 0
    end


    @testset "linearDB vs motifDB" begin
        mdb_path = joinpath(tdir, "motifDB")
        mkpath(mdb_path)
        build_motifDB("samirandom", genome, Motif("Cas9"), mdb_path, 7)
        detail_path = joinpath(mdb_path, "detail.csv")
        
        for d in 1:3
            search_motifDB(mdb_path, guides, detail_path; distance = d)
            mdb = DataFrame(CSV.File(detail_path))

            search_linearDB(ldb_path, guides, detail_path; distance = d)
            ldb = DataFrame(CSV.File(detail_path))
            failed = antijoin(ldb, mdb, on = [:guide, :distance, :chromosome, :start, :strand])
            @test nrow(failed) == 0
        end
    end


    @testset "linearDB vs fmiDB" begin
        fmi_dir = joinpath(tdir, "fmifDB")
        mkpath(fmi_dir)
        build_fmiDB(genome, fmi_dir)

        # build a pamDB
        motif = Motif("Cas9"; distance = 2)
        pamDB = build_pamDB(fmi_dir, motif)

        # prepare PathTemplates
        mpt = build_PathTemplates(motif; withPAM = true)

        # prepare output folder
        res_dir = joinpath(tdir, "results")
        mkpath(res_dir)

        # finally, make results!
        res_path = joinpath(res_dir, "results.csv")
        search_fmiDB(guides, mpt, fmi_dir, res_path; distance = 2)
        res_fmiDB = DataFrame(CSV.File(res_path))
        res_fmiDB = filter_overlapping(res_fmiDB, 23)
        select!(res_fmiDB, Not([:alignment_guide, :alignment_reference]))

        search_linearDB(ldb_path, guides, detail_path; distance = 2)
        ldb = DataFrame(CSV.File(detail_path))
        ldb = filter_overlapping(ldb, 23)
        select!(ldb, Not([:alignment_guide, :alignment_reference]))

        # test outputs for brute force method!
        failed = antijoin(ldb, res_fmiDB, on = [:guide, :distance, :chromosome, :start, :strand])
        @test nrow(failed) == 0

        # finally, make results!
        res_path = joinpath(res_dir, "results_seed.csv")
        search_fmiDB_seed(guides, fmi_dir, genome, pamDB, res_path; distance = 2)
        res_fmiDB_seed = DataFrame(CSV.File(res_path))
        res_fmiDB_seed = filter_overlapping(res_fmiDB_seed, 23)
        select!(res_fmiDB_seed, Not([:alignment_guide, :alignment_reference]))

        # test outputs for seed method!
        failed = antijoin(ldb, res_fmiDB_seed, on = [:guide, :distance, :chromosome, :start, :strand])
        @test nrow(failed) == 0
    end

    @testset "linearDB vs bffDB" begin
        fmi_dir = joinpath(tdir, "fmifDB")
        mkpath(fmi_dir)
        build_fmiDB(genome, fmi_dir)

        # build a pamDB
        motif = Motif("Cas9"; distance = 2)

        bff_dir = joinpath(tdir, "bffDB")
        mkpath(bff_dir)
        build_binaryFuseFilterDB("testBFF", genome, motif, bff_dir)

        # prepare output folder
        res_dir = joinpath(tdir, "resultsBFF")
        mkpath(res_dir)

        # finally, make results!
        res_path = joinpath(res_dir, "results.csv")
        search_binaryFuseFilterDB(bff_dir, fmi_dir, genome, guides, res_path; distance = 2)
        
        res_bffDB = DataFrame(CSV.File(res_path))
        res_bffDB = filter_overlapping(res_bffDB, 23)
        select!(res_bffDB, Not([:alignment_guide, :alignment_reference]))

        search_linearDB(ldb_path, guides, detail_path; distance = 2)
        ldb = DataFrame(CSV.File(detail_path))
        ldb = filter_overlapping(ldb, 23)
        select!(ldb, Not([:alignment_guide, :alignment_reference]))

        # test outputs for brute force method!
        failed = antijoin(ldb, res_bffDB, on = [:guide, :distance, :chromosome, :start, :strand])
        @test nrow(failed) == 0
    end


    @testset "linearDB vs linearDB early stopped" begin
        ldb_filt = DataFrame(CSV.File(detail_path))
        ldb_filt = filter_overlapping(ldb_filt, 3*2 + 1)
        ldb_res_filt = summarize_offtargets(ldb_filt; distance = 3)

        detail_path_es = joinpath(ldb_path, "detail_es.csv")
        # find all offtargets with overlap filtering on the go
        search_linearDB_with_es(ldb_path, guides, detail_path_es; distance = 3, early_stopping = [50, 50, 50, 50])
        ldbes = DataFrame(CSV.File(detail_path_es))
        ldbes_res = summarize_offtargets(ldbes; distance = 3)
        @test compare_result(ldb_res_filt, ldbes_res; less_or_equal = true)

        # find all offtargets with es with overlap filtering
        search_linearDB_with_es(ldb_path, [dna"NNNNNNNNNNNNNNNNNNNN"], 
            detail_path_es; distance = 3, early_stopping = repeat([2], 4))
        ldbes = DataFrame(CSV.File(detail_path_es))
        # because of the Thread break there it is highly non-deterministic how many offtargets we will get
        @test nrow(ldbes) >= 2
    end


    @testset "linearDB vs prefixHashDB early stopped" begin
        # remember that this early stopping can find overlaps
        search_linearDB(ldb_path, guides, detail_path; distance = 3)
        ldb = DataFrame(CSV.File(detail_path))

        phdb_path = joinpath(tdir, "prefixHashDBes")
        mkpath(phdb_path)
        build_prefixHashDB("samirandom", genome, setdist(Motif("Cas9"), 3), phdb_path)
        detail_path_es = joinpath(phdb_path, "detail_es.csv")

        search_prefixHashDB(phdb_path, guides, detail_path_es; distance = 3, 
            early_stopping = [300, 300, 300, 300])
        pdbes = DataFrame(CSV.File(detail_path_es))
        failed = antijoin(ldb, pdbes, on = [:guide, :distance, :chromosome, :start, :strand])
        @test nrow(failed) == 0
        
        # find all offtargets with es with overlap filtering - we dont support ambigous guides anymore
        search_prefixHashDB(phdb_path, guides, 
            detail_path_es; distance = 0, early_stopping = [0])
        pdbes = DataFrame(CSV.File(detail_path_es))
        pdbes_res = summarize_offtargets(pdbes)
        @test nrow(pdbes) == 20
    end
end