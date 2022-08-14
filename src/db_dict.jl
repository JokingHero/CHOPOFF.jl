struct DictDB
    dict::IdDict
    mtp::MotifPathTemplates
    dbi::DBInfo
    ambig::Union{AmbigIdx, Nothing}
end


function build_guide_dict(dbi::DBInfo, max_count::Int, guide_type::Type{T}) where T <: Union{UInt64, UInt128}
    max_count_type = smallestutype(unsigned(max_count))
    guides = Vector{guide_type}()
    ambig = gatherofftargets!(guides, dbi)
    ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing
    guides = sort(guides)
    guides, counts = ranges(guides)
    counts = convert.(max_count_type, min.(length.(counts), max_count))
    return IdDict{guide_type, max_count_type}(zip(guides, counts)), ambig
end


function build_dictDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String)

    if motif.dist > 2 
        @info "Searches on distances of more than 2 will be very slow! Be warned!"
    end

    dbi = DBInfo(genomepath, name, motif)
    @info "Building Dictionary..."
    dict, ambig = build_guide_dict(dbi, Int(typemax(UInt32)), UInt64)
    @info "Building Motif templates..."
    mtp = build_motifTemplates(motif)

    db = DictDB(dict, mtp, dbi, ambig)
    save(db, joinpath(storagedir, "dictDB.bin"))
    @info "Finished constructing dictDB in " * storagedir
    @info "Database size is:" *
        "\n length -> " * string(length(db.dict)) *
        "\n consuming: " * string(round((filesize(joinpath(storagedir, "dictDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


function expand_path(path::Path, len::Int)
    path_seq = path.seq * repeat(dna"N", len - path.seq_len)
    return expand_ambiguous(path_seq)
end


function search_dictDB(
    storagedir::String,
    guides::Vector{LongDNA{4}})

    sdb = load(joinpath(storagedir, "dictDB.bin"))
    guides_ = copy(guides)
    dist = sdb.mtp.motif.distance # use maximal distance as the performance is always bottlenecked by that

    if any(length_noPAM(sdb.dbi.motif) .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(sdb.dbi.motif))
    end
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = zeros(Int, length(guides_), dist + 1)
    for (i, s) in enumerate(guides_)

        pat = ARTEMIS.templates_to_sequences_by_dist(s, sdb.mtp)
        for di in 1:(dist + 1)
            res[i, di] = Base.mapreduce(x -> get(sdb.dict, convert(UInt64, x), 0), +, pat[di])

            if !isnothing(sdb.ambig)
                res[i, di] += sum(Base.mapreduce(x -> findbits(x, sdb.ambig), .|, pat[di]))
            end
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end