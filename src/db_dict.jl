struct DictDB
    dict::IdDict
    dbi::DBInfo
end


function add_guides!(dict::IdDict, guides::Vector{LongDNASeq})
    ktype = valtype(dict)
    for value in guides
        if isambiguous(value)
            value = expand_ambiguous(value)
            for v in value
                dict[v] = safeadd(get(dict, v, convert(ktype, 0)), convert(ktype, 1))
            end
        else
            dict[value] = safeadd(get(dict, value, convert(ktype, 0)), convert(ktype, 1))
        end
    end
    return dict
end


function build_guide_dict(dbi::DBInfo, max_count::Int)
    max_count_type = smallestutype(unsigned(max_count))
    guides = Vector{LongDNASeq}()
    gatherofftargets!(guides, dbi)
    guides = sort(guides)
    guides, counts = ranges(guides)
    counts = convert.(max_count_type, min.(length.(counts), max_count))
    return IdDict{LongDNASeq, max_count_type}(zip(guides, counts))
end


function build_dictDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String;
    max_count::Int = 255)

    dbi = DBInfo(genomepath, name, motif)

    # first we measure how many unique guides there are
    @info "Building Dictionary..."
    dict = build_guide_dict(dbi, max_count)

    db = DictDB(dict, dbi)
    save(db, joinpath(storagedir, "dictDB.bin"))
    @info "Finished constructing dictDB in " * storagedir
    @info "Database size is:" *
        "\n length -> " * string(length(db.dict)) *
        "\n consuming: " * string(round((filesize(joinpath(storagedir, "dictDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


function search_dictDB(
    storagedir::String,
    guides::Vector{LongDNASeq},
    dist::Int = 1)

    if any(isambiguous.(guides))
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    sdb = load(joinpath(storagedir, "dictDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), (dist + 1) * 3 - 2)
    for (i, s) in enumerate(guides_)
        res[i, 1] += get(sdb.dict, s, 0) # 0 distance
        for d in 1:dist
            norm_d, border_d = comb_of_d(string(s), d)
            norm_d_res = ThreadsX.sum(get(sdb.dict, LongDNASeq(sd), 0) for sd in norm_d)
            border_d_res = ThreadsX.sum(get(sdb.dict, LongDNASeq(sd), 0) for sd in border_d)
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