struct BinDB{T<:Unsigned}
    bins::Vector{BloomFilter}
    counts::Vector{T}
    dbi::DBInfo
end


function estimate(db::BinDB, guide::UInt64)
    for i in length(db.counts):-1:1
        if guide in db.bins[i]
            return db.counts[i]
        end
    end
    return 0
end


function estimate(db::BinDB, guide::Union{LongDNASeq, String})
    guide = unsigned(DNAMer(guide))
    return estimate(db, guide)
end


function build_binDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String,
    dictDB::String = "";
    probability_of_error::Float64 = 0.001,
    max_count::Int = 10)

    dbi = DBInfo(genomepath, name, motif)

    # first we measure how many unique guides there are
    max_count_type = smallestutype(unsigned(max_count))
    if dictDB == ""
        @info "Building Dictionary..."
        dict = IdDict{UInt64, max_count_type}()
        gatherofftargets!(dict, dbi)
    else 
        dict = load(joinpath(dictDB, "dictDB.bin"))
        dict = dict.dict
        # TODO make sure dict type is max_count_type
    end
    max_count = convert(max_count_type, max_count)

    counts = IdDict{UInt8, Int}()
    for v in values(dict)
        if v >= max_count
            counts[max_count] = get(counts, max_count, 0) + 1
        else
            counts[v] = get(counts, v, 0) + 1
        end
    end

    # we sort size of the key from smallest to biggest
    order = sortperm(collect(keys(counts)))
    counts = collect(counts)
    counts = counts[order]

    # now for every single count we build BloomFilters
    bins = Vector{BloomFilter}()
    for c in counts
        params = constrain(
            BloomFilter, 
            fpr = probability_of_error, 
            capacity = c.second)
        push!(bins, BloomFilter(params.m, params.k))
    end
    counts = [x.first for x in counts]

    for (guide, real_count) in dict
        if real_count > max_count
            real_count = max_count
        end
        idx = findfirst(x -> x == real_count, counts)
        push!(bins[idx], guide)
    end
    
    db = BinDB(bins, counts, dbi)
    save(db, joinpath(storagedir, "binDB.bin"))

    # use full db to estimate error rate!
    conflict = 0
    error = Vector()
    for (guide, real_count) in dict
        est_count = CRISPRofftargetHunter.estimate(db, guide)
        if real_count >= max_count
            real_count = max_count
        end
        if  real_count != est_count
            conflict += 1
            push!(error, Int(est_count) - Int(real_count))
        end
    end
    error_rate = conflict / length(dict)
    mean_err = mean(error)

    @info "Finished constructing binDB in " * storagedir
    @info "Estimated probability of miscounting an element in the bins is " * string(round(error_rate; digits = 6))
    @info "Mean error was " * string(mean_err)
    @info "Database is consuming: " * 
        string(round((filesize(joinpath(storagedir, "binDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


"""
`search_binDB(
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
function search_binDB(
    storagedir::String,
    guides::Vector{LongDNASeq},
    dist::Int = 1)

    bdb = load(joinpath(storagedir, "binDB.bin"))
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