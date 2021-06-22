struct DictDB
    dict::IdDict
    dbi::DBInfo
end

function add_guides!(dict::IdDict, guides::Vector{LongDNASeq})
    for value in guides
        value = unsigned(DNAMer(value))
        dict[value] = get(dict, value, 0) + 1
    end
    return dict
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
    max_count_type = smallestutype(unsigned(max_count))
    dict = IdDict{UInt64, max_count_type}()
    gatherofftargets!(dict, dbi)
    
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

    sdb = load(joinpath(storagedir, "dictDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), (dist + 1) * 3 - 2)
    for (i, s) in enumerate(guides_)
        res[i, 1] += get(sdb.dict, unsigned(DNAMer(s)), 0) # 0 distance
        for d in 1:dist
            norm_d, border_d = comb_of_d(string(s), d)
            norm_d_res = ThreadsX.sum(get(sdb.dict, unsigned(DNAMer(LongDNASeq(sd))), 0) for sd in norm_d)
            border_d_res = ThreadsX.sum(get(sdb.dict, unsigned(DNAMer(LongDNASeq(sd))), 0) for sd in border_d)
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