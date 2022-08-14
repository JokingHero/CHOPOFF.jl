struct BinDB{T<:Unsigned}
    bins::Vector{BloomFilter}
    counts::Vector{T}
    mtp::MotifPathTemplates
    dbi::DBInfo
    ambig::Union{AmbigIdx, Nothing}
end


function estimate(db::BinDB, guide::LongDNA{4}, right::Bool)
    direction = right ? (1:length(db.bins)) : (length(db.bins):-1:1)
    for i in direction
        if guide in db.bins[i]
            return db.counts[i]
        end
    end
    return 0
end


function estimate(db::BinDB, guide::String, right::Bool)
    return estimate(db, LongDNA{4}(guide), right)
end


"""
```
build_binDB(
    name::String,
    genomepath::String, 
    motif::Motif,
    storagedir::String,
    dictDB::String = "";
    probability_of_error::Float64 = 0.001,
    max_count::Int = 10)
```

Prepare binDB index for future searches using `search_binDB`.


# Arguments
`name` - Your prefered name for this index to ease future identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif defines what kind of gRNA to search for.

`storagedir`  - Folder path to the where index will be saved with name `binDB.bin`.

`dictDB`  - Optional. Supply folder directory where dictDB.bin is located, bypasses normal DB build and reuses dictDB dataset (quicker).

`probability_of_error` - What is the false positive rate that we would like to achieve? 
             
`max_count`  - Above this count we put all unique off-targets into one bin. 
               Put number here that is the minimum number of off-targets that you think is fine
               for the distance of 1.


# Examples
```julia-repl
# make a temporary directory
tdir = tempname()
bdb_path = joinpath(tdir, "binDB")
mkpath(bdb_path)

# use ARTEMIS example genome
genome = joinpath(
    vcat(
        splitpath(dirname(pathof(ARTEMIS)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))

# finally, build a binDB
build_binDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0), 
    bdb_path)
```
"""
function build_binDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String,
    dictDB::String = "";
    probability_of_error::Float64 = 0.001,
    max_count::Int = 10)

    if motif.distance != 1
        @info "Distance and ambig_max enforced to 1."
        motif = setdist(motif, 1)
    end
    dbi = DBInfo(genomepath, name, motif)

    # first we measure how many unique guides there are
    max_count_type = smallestutype(unsigned(max_count))
    if dictDB == ""
        @info "Building Dictionary..."
        dict, ambig = build_guide_dict(dbi, max_count, UInt64)
    else 
        dict = load(joinpath(dictDB, "dictDB.bin"))
        dict = dict.dict
        ambig = dict.ambig
        ktype = valtype(dict)
        if max_count > typemax(ktype)
            throw("Dictionary supports only up to " * string(typemax(ktype)) * 
                " max_count and current max_count is " * string(max_count))
        end
    end
    @info "Building Motif templates..."
    mtp = build_motifTemplates(motif)
    max_count = convert(max_count_type, max_count)

    counts = IdDict{UInt8, Int}()
    for v in values(dict)
        if v >= max_count
            counts[max_count] = get(counts, max_count, 0) + 1
        else
            counts[v] = get(counts, v, 0) + 1
        end
    end

    # we sort size of the key from small to large
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
    len_noPAM = length_noPAM(dbi.motif) + dbi.motif.distance

    # fill up the bins with our guides
    for (guide, real_count) in dict
        if real_count > max_count
            real_count = max_count
        end
        idx = findfirst(x -> x == real_count, counts)
        push!(bins[idx], LongDNA{4}(guide, len_noPAM))
    end
    
    db = BinDB(bins, counts, mtp, dbi, ambig)
    save(db, joinpath(storagedir, "binDB.bin"))

    # use full db to estimate error rate!
    error_right = 0
    error_left = 0
    for (guide, real_count) in dict
        guide = LongDNA{4}(guide, len_noPAM)
        if iscertain(guide)
            if real_count >= max_count
                real_count = max_count
            end
            error_right += Int(estimate(db, guide, true) != real_count)
            error_left += Int(estimate(db, guide, false) != real_count)
        end
    end

    @info "Finished constructing binDB in " * storagedir * " consuming "  * 
        string(round((filesize(joinpath(storagedir, "binDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    @info "Estimated probability of miscounting an element in the bins is - Right: " * 
        string(round(error_right / length(dict); digits = 6)) * " Left: " * 
        string(round(error_left / length(dict); digits = 6))
    @info "Database is consuming: " * 
        string(round((filesize(joinpath(storagedir, "binDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


"""
```
search_binDB(
    storagedir::String,
    guides::Vector{LongDNA{4}},
    right::Bool)
```

Estimate off-target counts for `guides` using hashDB stored at `storagedir`.

Probabilistic filter offers a guarantee that it will always be correct when a sequence 
is in the set (no false negatives), but may overestimate that a sequence is in the set 
while it is not (false positive) with low probability. If both columns in the results are 0, 
it is guaranteed this gRNA has no off-targets in the genome!

Also, maximum count for each off-target in the database is capped,
to `max_count` specified during building of hashDB. This means 
that counts larger than `max_count` are no longer estimating correctly.
Its likely you would not care for those guides anyhow.

`right` argument specifies whether the databse should be checked in direction from 
unique off-targets which occur once to increasingly more occuring off-targets up until 
`max_count` is reached, which may result in assuming lower than real off-target counts 
(underestimate) for some of the sequences, however this approach will not reject any 
gRNAs that should not be rejected and is suitable for filtering of gRNAs we do not need. 
Left (`right = false`, or hight-counts to low-counts) approach is also supported, 
which can be used for ordering of gRNAs to the best of database ability. 
Left approach may overestimate counts for some gRNAs. When gRNA is reported as 
off-target free it is also guaranteed to be true in both cases (low-to-high and high-to-low).

# Examples
```julia-repl
# prepare libs
using ARTEMIS, BioSequences

# make a temporary directory
tdir = tempname()
hdb_path = joinpath(tdir, "hashDB")
mkpath(hdb_path)

# use ARTEMIS example genome
coh_path = splitpath(dirname(pathof(ARTEMIS)))[1:end-1]
genome = joinpath(vcat(coh_path, "test", "sample_data", "genome", "semirandom.fa"))

# build a hashDB
build_hashDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0), 
    hdb_path)

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(coh_path, "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)

# finally, get results!
hdb_res = search_hashDB(hdb_path, guides, false)
```
"""
function search_binDB(
    storagedir::String,
    guides::Vector{LongDNA{4}},
    right::Bool)

    if any(isambig.(guides))
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    bdb = load(joinpath(storagedir, "binDB.bin"))
    guides_ = copy(guides)
    len_noPAM_noEXT = length_noPAM(bdb.dbi.motif)
    len = len_noPAM_noEXT + bdb.dbi.motif.distance

    if any(len_noPAM_noEXT .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(bdb.dbi.motif))
    end
    # reverse guides so that PAM is always on the left
    if bdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), 2)
    for (i, s) in enumerate(guides_)

        pat = ARTEMIS.templates_to_sequences(s, bdb.mtp; dist = 1, reducible = false)
        d0 = Set(expand_path(pat[1], len))
        d1 = Set(mapreduce(x -> expand_path(x, len), vcat, pat[2:end]))
        d1 = setdiff(d1, d0) # remove d0 from d1
        d0 = collect(d0)
        d1 = collect(d1)

        res[i, 1] = Base.mapreduce(x -> estimate(bdb, x, right), +, d0)
        res[i, 2] = Base.mapreduce(x -> estimate(bdb, x, right), +, d1)

        if !isnothing(bdb.ambig)
            bits_mapped = Base.map(x -> findbits(x, bdb.ambig), d0)
            res[i, 1] += sum(reduce(.|, bits_mapped))
            bits_mapped = Base.map(x -> findbits(x, bdb.ambig), d1)
            res[i, 2] +=  sum(reduce(.|, bits_mapped))
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:1]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end