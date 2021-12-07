struct NoHashDB
    dbi::DBInfo
    guides::Vector{BitVector}
    counts::Vector{UInt8}
    ambig::Union{AmbigIdx, Nothing}
end


"
Prepare column wise storage of the guides.
"
function column_transform(guides::Vector{UInt64}, len::Int)
    guides = sort(guides)
    cols = Vector{BitVector}()
    for i in (len*2):-1:1
        push!(cols, ((guides .<< (64 - i)) .>> 63) .!= 0)
    end
    return cols
end


"Find which guides are matched with the `seq`."
function findbits(seq::LongDNASeq, cols::Vector{BitVector})
    seq_len = length(seq) * 2
    seq = BitVector(digits(convert(UInt64, seq), base=2, pad=64) |> reverse)
    idx = BitVector(ones(Bool, length(cols[1])))
    j = 1
    for i in seq_len:-1:1
        if seq[65 - i] # is true
            idx = idx .& cols[j]
        else
            idx = idx .& (.!cols[j])
        end
        j += 1
    end
    return idx
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

    if (motif.distance != 1)
        @info "Distance enforced to 1."
        motif = setdist(motif, 1)
    end

    dbi = DBInfo(genomepath, name, motif)

    # first we measure how many unique guides there are
    @info "Building Dictionary..."
    guides, counts, ambig = get_guides_by_counts(dbi)
    ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing
    guides = column_transform(guides, motif.distance + length_noPAM(motif))

    db = NoHashDB(dbi, guides, counts, ambig)
    save(db, joinpath(storagedir, "noHashDB.bin"))
    @info "Finished constructing noHashDB in " * storagedir
    @info "Database size is:" *
        "\n length -> " * string(length(db.guides)) *
        "\n consuming: " * string(round((filesize(joinpath(storagedir, "noHashDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


function search_noHashDB(
    storagedir::String,
    guides::Vector{LongDNASeq})

    if any(isambig.(guides))
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
        idx0 = findbits(s, sdb.guides)
        res[i, 1] += isnothing(sdb.ambig) ? 0 : sum(findbits(s, sdb.ambig))
        res[i, 1] += sum(sdb.counts[idx0]) # 0 distance
        
        d1_combs = LongDNASeq.(comb_of_d1_extended(string(s))) # 1 distance
        idx1 = reduce(.|, map(x -> findbits(x, sdb.guides), d1_combs))
        idx1 = idx1 .& .!idx0 # remove 0D index inside 1D indexes
        res[i, 2] = sum(sdb.counts[idx1])
        if !isnothing(sdb.ambig)
            res[i, 2] += sum(reduce(.|, map(x -> findbits(x, sdb.ambig), d1_combs)))
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:1]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end