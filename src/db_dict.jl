struct DictDB
    dict::IdDict
    mtp::PathTemplates
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
    motif::Motif;
    storage_path::String = "")

    if motif.distance > 2 
        @info "Searches on distances of more than 2 will be very slow! Be warned!"
    end

    dbi = DBInfo(genomepath, name, motif)
    @info "Building Dictionary..."
    dict, ambig = build_guide_dict(dbi, Int(typemax(UInt32)), UInt64)
    @info "Building Motif templates..."
    mtp = build_PathTemplates(length_noPAM(motif), motif.distance)

    db = DictDB(dict, mtp, dbi, ambig)
    if storage_path != ""
        save(db, storage_path)
        @info "Finished constructing dictDB in " * storage_dir
        @info "Database size is:" *
            "\n length -> " * string(length(db.dict)) *
            "\n consuming: " * string(round((filesize(storage_path) * 1e-6); digits = 3)) * 
            " mb of disk space."
    end
    return db
end


function search_dictDB(
    db::DictDB,
    guides::Vector{LongDNA{4}})

    guides_ = copy(guides)
    dist = db.mtp.distance # use maximal distance as the performance is always bottlenecked by that

    if any(length_noPAM(db.dbi.motif) .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(db.dbi.motif))
    end
    # reverse guides so that PAM is always on the left
    if db.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = zeros(Int, length(guides_), dist + 1)
    for (i, s) in enumerate(guides_)

        pat = ARTEMIS.templates_to_sequences_extended(s, db.mtp)
        for di in 1:(dist + 1)
            res[i, di] = Base.mapreduce(x -> get(db.dict, convert(UInt64, x), 0), +, pat[di])

            if !isnothing(db.ambig)
                res[i, di] += sum(Base.mapreduce(x -> findbits(x, db.ambig), .|, pat[di]))
            end
        end
    end

    res = format_DF(res, dist, guides)
    return res
end