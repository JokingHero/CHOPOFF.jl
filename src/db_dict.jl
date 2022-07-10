struct DictDB
    dict::IdDict
    mtp::MotifPathTemplates
    dbi::DBInfo
end


function build_guide_dict(dbi::DBInfo, max_count::Int, guide_type::Type{T}) where T <: Union{UInt64, UInt128}
    max_count_type = smallestutype(unsigned(max_count))
    guides = Vector{guide_type}()
    gatherofftargets!(guides, dbi) # we ignore ambig
    guides = sort(guides)
    guides, counts = ranges(guides)
    counts = convert.(max_count_type, min.(length.(counts), max_count))
    return IdDict{guide_type, max_count_type}(zip(guides, counts))
end


function build_dictDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String)

    if motif.distance != 1 || motif.ambig_max != 0
        @info "Distance set to 1 and ambig_max set to 0."
        motif = setdist(motif, 1)
        motif = setambig(motif, 0)
    end
    dbi = DBInfo(genomepath, name, motif)
    @info "Building Dictionary..."
    dict = build_guide_dict(dbi, Int(typemax(UInt32)), UInt128)
    @info "Building Motif templates..."
    mtp = build_motifTemplates(motif)

    db = DictDB(dict, mtp, dbi)
    save(db, joinpath(storagedir, "dictDB.bin"))
    @info "Finished constructing dictDB in " * storagedir
    @info "Database size is:" *
        "\n length -> " * string(length(db.dict)) *
        "\n consuming: " * string(round((filesize(joinpath(storagedir, "dictDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


function get_path_count(dict::IdDict, path::Path, len::Int)
    path_seq = path.seq * repeat(dna"N", len - path.seq_len)
    ext_path = expand_ambiguous(path_seq)
    return mapreduce(x -> get(dict, convert(UInt128, x), 0), +, ext_path)
end


function search_dictDB(
    storagedir::String,
    guides::Vector{LongDNA{4}})

    if any(isambig.(guides))
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    sdb = load(joinpath(storagedir, "dictDB.bin"))
    guides_ = copy(guides)
    len = length(sdb.dbi.motif)
    len_noPAM_noEXT = length_noPAM(sdb.dbi.motif)

    if any(len_noPAM_noEXT .!= length(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(motif))
    end
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = zeros(Int, length(guides_), 2)
    for (i, s) in enumerate(guides_)

        pat = CRISPRofftargetHunter.templates_to_sequences(s, sdb.mtp; dist = 1)
        counts = Base.map(x -> get_path_count(sdb.dict, x, len), pat)
        # we need to keep track of reducibles
        # reducible should not include counts from that path
        # this is simpler as we only have distances of 0 and 1 to worry about
        d0_loc = 0
        for (j, p) in enumerate(pat)
            if p.reducible != 0 && (p.dist < pat[p.reducible].dist)
                counts[p.reducible] -= counts[j]
            end

            if p.dist == 0
                d0_loc = j
            end
        end
        
        res[i, 1] = counts[d0_loc]
        res[i, 2] = sum(counts) - counts[d0_loc]
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D0"), Symbol("D1")]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end