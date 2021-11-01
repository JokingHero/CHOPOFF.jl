struct NoHashDB
    dbi::DBInfo
    guides::Vector{UInt64}
    counts::Vector{UInt8}
    ambig::Union{AmbigIdx, Nothing}
end


function get_guides_by_counts(dbi::DBInfo)
    guides = Vector{UInt64}()
    ambig = gatherofftargets!(guides, dbi)
    guides = sort(guides)
    guides, counts = ranges(guides)
    counts = convert.(UInt8, min.(length.(counts), 255))
    return guides, counts, ambig
end


function build_noHashDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String)

    if motif.ambig_max != 0
        throw("We don't support unambigous searches yet.")
    end
    if motif.distance != 1
        throw("Motif for distances not 1 is not supported.")
    end
    dbi = DBInfo(genomepath, name, motif)

    # first we measure how many unique guides there are
    @info "Building Dictionary..."
    guides, counts, ambig = get_guides_by_counts(dbi)
    ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing

    db = NoHashDB(dbi, guides, counts, ambig)
    save(db, joinpath(storagedir, "noHashDB.bin"))
    @info "Finished constructing noHashDB in " * storagedir
    @info "Database size is:" *
        "\n length -> " * string(length(db.guides)) *
        "\n consuming: " * string(round((filesize(joinpath(storagedir, "noHashDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


@inline function is_equal(guide::UInt64, query::UInt64, shift::Int64)
    return count_ones(~((guide >> shift) ⊻ query)) == 64
end


function search_noHashDB(
    storagedir::String,
    guides::Vector{LongDNASeq})

    if any(n_ambiguous.(guides) .> 0)
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    sdb = load(joinpath(storagedir, "noHashDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), 2)
    len_noPAM = length_noPAM(sdb.dbi.motif)
    for (i, s) in enumerate(guides_)
        idx0 = ThreadsX.findfirst(x -> is_equal(x, convert(UInt64, s), 2), sdb.guides)
        res[i, 1] += isnothing(sdb.ambig) ? 0 : sum(findbits(s, sdb.ambig))
        if !isnothing(idx0)
            res[i, 1] += sdb.counts[idx0] # 0 distance
        end
        
        d1_combs = LongDNASeq.(comb_of_d1_extended(string(s))) # 1 distance
        idx1 = Set{Int64}()
        for d1_comb in d1_combs
            idx1 = union(idx1, ThreadsX.findall(x -> is_equal(x, convert(UInt64, d1_comb), (len_noPAM - length(d1_comb)) * 2), sdb.guides))
        end
        
        if !isnothing(idx0)
            idx1 = setdiff(idx1, Set(idx0)) # remove 0D index inside 1D indexes
        end
        
        idx1 = collect(idx1)
        res[i, 2] = sum(sdb.counts[idx1])
        if !isnothing(sdb.ambig)
            res[i, 2] += sum(reduce(|, map(x -> findbits(x, sdb.ambig), d1_combs)))
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:1]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end