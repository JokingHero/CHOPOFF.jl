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

    if motif.distance != 1
        @info "Distance set to 1."
        motif = setdist(motif, 1)
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

    if any(isambig.(guides))
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    sdb = load(joinpath(storagedir, "dictDB.bin"))
    guides_ = copy(guides)
    len_noPAM_noEXT = length_noPAM(sdb.dbi.motif)
    len = len_noPAM_noEXT + sdb.dbi.motif.distance

    if any(len_noPAM_noEXT .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(sdb.dbi.motif))
    end
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = zeros(Int, length(guides_), 2)
    for (i, s) in enumerate(guides_)

        pat = CRISPRofftargetHunter.templates_to_sequences(s, sdb.mtp; dist = 1, reducible = false)
        d0 = Set(expand_path(pat[1], len))
        d1 = Set(mapreduce(x -> expand_path(x, len), vcat, pat[2:end]))
        d1 = setdiff(d1, d0) # remove d0 from d1
        d0 = collect(d0)
        d1 = collect(d1)

        res[i, 1] = Base.mapreduce(x -> get(sdb.dict, convert(UInt64, x), 0), +, d0)
        res[i, 2] = Base.mapreduce(x -> get(sdb.dict, convert(UInt64, x), 0), +, d1)

        if !isnothing(sdb.ambig)
            bits_mapped = Base.map(x -> findbits(x, sdb.ambig), d0)
            res[i, 1] += sum(reduce(.|, bits_mapped))
            bits_mapped = Base.map(x -> findbits(x, sdb.ambig), d1)
            res[i, 2] +=  sum(reduce(.|, bits_mapped))
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D0"), Symbol("D1")]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end