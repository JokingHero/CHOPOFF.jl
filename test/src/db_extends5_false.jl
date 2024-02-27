using Test

using ARTEMIS
using BioSequences
using CSV
using DataFrames

## SET WD when debugging
# cd("test")

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
    genome = joinpath(dirname(pathof(ARTEMIS)), "..", 
        "test", "sample_data", "genome", "semirandom.2bit")

    guides = LongDNA{4}.(["TCGATTGTTTGGCTCTCTAAA", "GCAGGGGGACGCAAGTACGAA", "GGGCCGAAACGCGACACCGCC"])
    tdir = tempname()
    mkpath(tdir)

    # make and run default vcfDB
    vcf = joinpath(dirname(pathof(ARTEMIS)), "..", 
        "test", "sample_data", "artificial.vcf")
    vcf_db = build_vcfDB(
        "samirandom", genome, vcf,
        Motif("Cas12a"; distance = 1, ambig_max = 0))
    vcf_res = search_vcfDB(vcf_db, guides)

    @testset "vcfDB result is empty" begin
        @test all(vcf_res.D0 .== 0 .& vcf_res.D1 .== 0)
    end

    # make and run default linearDB
    ldb_path = joinpath(tdir, "linearDB")
    mkpath(ldb_path)
    build_linearDB("samirandom", genome, Motif("Cas12a"), ldb_path, 7)
    detail_path = joinpath(ldb_path, "detail.csv")
    search_linearDB(ldb_path, guides, detail_path; distance = 3)
    ldb = DataFrame(CSV.File(detail_path))
    ldb_res = summarize_offtargets(ldb; distance = 3)

    # make and run default dictDB
    dictDB = build_dictDB(
        "samirandom", genome, 
        Motif("Cas12a"; distance = 2))
    ddb_res = search_dictDB(dictDB, guides)

    # make and run default hashDB
    hashDB = build_hashDB(
        "samirandom", genome, 
        Motif("Cas12a"; distance = 1, ambig_max = 0))
    hdb_res = search_hashDB(hashDB, guides, false)

    # hashDB but with ambig
    hashDBambig = build_hashDB(
        "samirandom", genome, 
        Motif("Cas12a"; distance = 1, ambig_max = 1))
    hdb_res2 = search_hashDB(hashDBambig, guides, false)

    len_noPAM = ARTEMIS.length_noPAM(Motif("Cas12a"))

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
        build_treeDB("samirandom", genome, Motif("Cas12a"), tdb_path, 7)
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
        build_motifDB("samirandom", genome, Motif("Cas12a"), mdb_path, 7)
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
        motif = Motif("Cas12a"; distance = 1)

        # prepare PathTemplates
        mpt = build_PathTemplates(motif; withPAM = true)

        # prepare output folder
        res_dir = joinpath(tdir, "results")
        mkpath(res_dir)

        # finally, make results!
        res_path = joinpath(res_dir, "results.csv")
        search_fmiDB(guides, mpt, fmi_dir, res_path; distance = 1)
        res_fmiDB = DataFrame(CSV.File(res_path))
        res_fmiDB = filter_overlapping(res_fmiDB, 23)
        select!(res_fmiDB, Not([:alignment_guide, :alignment_reference]))

        search_linearDB(ldb_path, guides, detail_path; distance = 1)
        ldb = DataFrame(CSV.File(detail_path))
        ldb = filter_overlapping(ldb, 23)
        select!(ldb, Not([:alignment_guide, :alignment_reference]))

        # test outputs for brute force method!
        failed = antijoin(ldb, res_fmiDB, on = [:guide, :distance, :chromosome, :start, :strand])
        @test nrow(failed) == 0
    end


    @testset "linearDB vs bffDB" begin
        fmi_dir = joinpath(tdir, "bffDB")
        mkpath(fmi_dir)
        build_fmiDB(genome, fmi_dir)

        # build a pamDB
        motif = Motif("Cas12a"; distance = 2)

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


    @testset "linearDB vs linearHashDB" begin
        lhdb_path = joinpath(tdir, "linearHashDB")
        mkpath(lhdb_path)
        build_linearHashDB("samirandom", genome, Motif("Cas12a"), lhdb_path, 7)
        detail_path = joinpath(lhdb_path, "detail.csv")
        
        for d in 1:3
            search_linearHashDB(lhdb_path, guides, detail_path; distance = d)
            lhdb = DataFrame(CSV.File(detail_path))

            search_linearDB(ldb_path, guides, detail_path; distance = d)
            ldb = DataFrame(CSV.File(detail_path))
            failed = antijoin(ldb, lhdb, on = [:guide, :distance, :chromosome, :start, :strand])
            @test nrow(failed) == 0
        end
    end


    @testset "linearDB vs linearHahsDB early stopped" begin
        # remember that this early stopping can find overlaps
        search_linearDB(ldb_path, guides, detail_path; distance = 3)
        ldb = DataFrame(CSV.File(detail_path))

        lhdb_path = joinpath(tdir, "linearHashDBes")
        mkpath(lhdb_path)
        build_linearHashDB("samirandom", genome, Motif("Cas12a"), lhdb_path, 7)
        detail_path_es = joinpath(lhdb_path, "detail_es.csv")

        # find all offtargets with overlap filtering on the go
        # search_linearHashDB(lhdb_path, guides, detail_path_es; distance = 3)
        # lhdb = DataFrame(CSV.File(detail_path))

        search_linearHashDB_with_es(lhdb_path, guides, detail_path_es; distance = 3, early_stopping = [300, 300, 300, 300])
        ldbes = DataFrame(CSV.File(detail_path_es))
        failed = antijoin(ldb, ldbes, on = [:guide, :distance, :chromosome, :start, :strand])
        @test nrow(failed) == 0
        
        # find all offtargets with es with overlap filtering - we dont support ambigous guides anymore
        search_linearHashDB_with_es(lhdb_path, guides, 
            detail_path_es; distance = 3, early_stopping = repeat([0], 4))
        ldbes = DataFrame(CSV.File(detail_path_es))
        ldbes_res = summarize_offtargets(ldbes)
        @test nrow(ldbes) == 6 # I checked these results
    end
end