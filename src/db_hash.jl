struct HashDB{T<:Unsigned}
    dbi::DBInfo
    bins_d0::Vector{Union{BinaryFuseFilter{UInt8}, Nothing}}
    bins_d1::Vector{Union{BinaryFuseFilter{UInt8}, Nothing}}
    bins_d2::Vector{Union{BinaryFuseFilter{UInt8}, Nothing}}
    ambig::AmbigIdx # if length = 0 then no ambig found
    #vcf::Vector{AmbigIdx}
    #vcf_names::Vector{String}
    counts::Vector{T}
end


function estimate(db::HashDB, guide::LongDNASeq) # TODO add direction
    for i in length(db.counts):-1:1
        if convert(UInt64, guide) in db.bins[i]
            return db.counts[i]
        end
    end
    return 0
end


function estimate(db::HashDB, guide::String)
    return estimate(db, LongDNASeq(guide))
end


function build_hashDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String,
    dictDB::String = "";
    seed::UInt64 = UInt64(0x726b2b9d438b9d4d),
    max_iterations::Int = 10,
    max_count::Int = 10)

    # assume we will allow search within distance of 2
    # extend motif to include 2 additional bases in the tail
    dbi = DBInfo(genomepath, name, motif)
    if motif.distance > 2
        throw("Only distance below 3 is supported!")
    
    # gather all unique off-targets
    max_count_type = smallestutype(unsigned(max_count))
    guides = Vector{guide_type}()
    gatherofftargets!(guides, dbi)
    guides = sort(guides)
    guides, counts = ranges(guides)
    counts = convert.(max_count_type, min.(length.(counts), max_count))



    max_count = convert(max_count_type, max_count)
    dict_counts = collect(values(dict))
    dict_guides = collect(keys(dict))

    bins = Vector{BinaryFuseFilter{UInt8}}()
    counts = Vector{max_count_type}()
    for count in 1:max_count
        idx = findall(x -> x == count, dict_counts)
        if length(idx) != 0
            push!(counts, convert(max_count_type, count))
            push!(bins, 
                BinaryFuseFilter{UInt8}(dict_guides[idx]; seed = seed, max_iterations = max_iterations))
        end
    end

    # TODO ambiguous
    db = nothing #HashDB(dbi, bins, counts)
    save(db, joinpath(storagedir, "hashDB.bin"))

    # use full db to estimate error rate!
    conflict = 0
    error = Vector{Int}()
    len_noPAM = length_noPAM(dbi.motif)
    for (guide, real_count) in dict
        guide = LongDNASeq(guide, len_noPAM)
        if n_ambiguous(guide) == 0
            est_count = estimate(db, guide)
            if real_count >= max_count
                real_count = max_count
            end
            if  real_count != est_count
                conflict += 1
                push!(error, Int(est_count) - Int(real_count))
            end
        end
    end
    error_rate = conflict / length(dict)

    @info "Finished constructing hashDB in " * storagedir * " consuming "  * 
        string(round((filesize(joinpath(storagedir, "hashDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    @info "Estimated probability of miscounting an element in the bins is " * string(round(error_rate; digits = 6))
    if length(error) > 1
        @info "Mean error was " * string(mean(error))
    end
    @info "Database is consuming: " * 
        string(round((filesize(joinpath(storagedir, "hashDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


"""
`search_hashDB(
    storagedir::String,
    guides::Vector{LongDNASeq},
    dist::Int = 1)`

Results are estimations of offtarget counts in the genome.

If CMS column is 0, it is guaranteed this guide has 
no 0-distance off-targets in the genome!

Also, maximum count for each off-target in the database is capped,
to `max_count` specified during building of binsDB. This means 
that counts larger than `max_count` (default 50) are no longer
estimating correctly, but we know at least 50 off-targets exist in 
the genome. Its likely you would not care for those guides anyhow.

In principle, you can use `distance` of any size, but for each increasing distance value
we have recurent solution that slows down significantly from distance 3 and upwards. These estimations
should be used normally just as prefiltering step for the guideRNA's, therefore checking just distance 0,
and distance 1 off-targets should be enough to rank the guides. For increasing distances partial errors in
estimations are stacking and overwhelm real estimations quickly.

## Result

Returns .csv table with the following columns:

D0, D1 ... D`distance` - these are estimated counts of off-targets for given guideRNA, 
each column for increasing `distance`, these is also a sum of the corresponding DN and DB columns

DN1, DN2, ... DN`distance`- these values are estimated counts of off-targets for given guideRNA,
each column for increasing `distance`, but these values are for guaranteed "within" distance off-targets

DB1, DB2, ... DB`distance` - same as above, but these values are representing "borderline cases" where 
genomic extension is not known and only asummed the worse case. These values are assuming worst case, 
and therefore are overestimating extensively with increasing `distance`.
"""
function search_hashDB(
    storagedir::String,
    guides::Vector{LongDNASeq},
    dist::Int = 1)

    if any(n_ambiguous.(guides) .> 0)
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    bdb = load(joinpath(storagedir, "hashDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if bdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), (dist + 1) * 3 - 2)
    for (i, s) in enumerate(guides_)
        res[i, 1] += estimate(bdb, s) # 0 distance
        for d in 1:dist
            norm_d, border_d = comb_of_d(string(s), d)
            norm_d_res = ThreadsX.sum(estimate(bdb, sd) for sd in norm_d)
            border_d_res = ThreadsX.sum(estimate(bdb, sd) for sd in border_d)
            res[i, d + 1] = norm_d_res + border_d_res
            res[i, dist + d + 1] = norm_d_res
            res[i, dist * 2 + d + 1] = border_d_res
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:dist]
    all_col_d = vcat(col_d, [Symbol("DN$i") for i in 1:dist], [Symbol("DB$i") for i in 1:dist])
    rename!(res, all_col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end