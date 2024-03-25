struct DictDB
    dict::IdDict
    mpt::PathTemplates
    dbi::DBInfo
    ambig::Union{AmbigPrefixHashDB, Nothing}
end


function build_guide_dict(dbi::DBInfo, max_count::Int, guide_type::Type{T}) where T <: Union{UInt64, UInt128}
    max_count_type = smallestutype(unsigned(max_count))
    guides = Vector{guide_type}()
    ambig = gatherofftargets!(guides, dbi)
    ambig = nothing # TODO
    # ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing
    guides = sort(guides)
    guides, counts = ranges(guides)
    counts = convert.(max_count_type, min.(length.(counts), max_count))
    return IdDict{guide_type, max_count_type}(zip(guides, counts)), ambig
end


"""
```
build_dictDB(
    name::String, 
    genomepath::String, 
    motif::Motif;
    storage_path::String = "")
```

Prepare dictDB index for future searches using `search_dictDB`.


# Arguments
`name` - Your preferred name for this index to ease future identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif defines what kind of gRNA to search for.

`storage_path`  - Path to the where index will be saved.


# Examples
```julia
# prepare libs
using CHOPOFF, BioSequences

# use CHOPOFF example genome
genome = joinpath(vcat(splitpath(dirname(pathof(CHOPOFF)))[1:end-1], 
    "test", "sample_data", "genome", "semirandom.fa"))

# build a hashDB!
db = build_dictDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0))
```
"""
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
    mpt = build_PathTemplates(motif)

    db = DictDB(dict, mpt, dbi, ambig)
    if storage_path != ""
        save(db, storage_path)
        @info "Finished constructing dictDB in " * storage_path
        @info "Database size is:" *
            "\n length -> " * string(length(db.dict)) *
            "\n consuming: " * string(round((filesize(storage_path) * 1e-6); digits = 3)) * 
            " mb of disk space."
    end
    return db
end


"""
```
search_dictDB(
    db::DictDB,
    guides::Vector{LongDNA{4}})
```

Summarize off-target counts for `guides` using `dictDB`.

This is simple dictionary storing all possible off-targets and their counts for a given `Motif`.
If you find no off-targets using this method it is guaranteed this gRNA has no off-targets in the genome!
Beware that the dictionary can be very big (e.g. human genome ~ 8Gb).

# Examples
```julia
# prepare libs
using CHOPOFF, BioSequences

# use CHOPOFF example genome
CHOPOFF_path = splitpath(dirname(pathof(CHOPOFF)))[1:end-1]
genome = joinpath(vcat(CHOPOFF_path, 
    "test", "sample_data", "genome", "semirandom.fa"))

# build a dictDB
db = build_dictDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0))

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(CHOPOFF_path, 
    "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)

# finally, get results!
res = search_dictDB(db, guides)
```
"""
function search_dictDB(
    db::DictDB,
    guides::Vector{LongDNA{4}})

    if any(isambig.(guides)) # TODO
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    guides_ = copy(guides)
    dist = db.mpt.motif.distance # use maximal distance as the performance is always bottlenecked by that
    # mpt = restrictDistance(db.mpt, dist)

    if any(length_noPAM(db.mpt.motif) .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(db.dbi.motif))
    end
    # reverse guides so that PAM is always on the left
    if db.mpt.motif.extends5
        guides_ = reverse.(guides_)
    end

    guides_uint2 = guide_to_template_format.(copy(guides_); alphabet = ALPHABET_TWOBIT)
    res = zeros(Int, length(guides_), dist + 1)
    for (i, s) in enumerate(guides_)

        pat = guides_uint2[i][db.mpt.paths]
        pat = map(x -> asUInt(UInt64, x), eachrow(pat))
        # further reduce non-unique seqeunces
        uniq = .!duplicated(pat)
        pat = pat[uniq]
        distances = db.mpt.distances[uniq]

        for di in 0:dist
            res[i, di + 1] = Base.mapreduce(x -> get(db.dict, x, 0), +, pat[distances .== di]; init = 0)

            #if !isnothing(db.ambig) # TODO?!
            #    res[i, di] += sum(Base.mapreduce(x -> findbits(x, db.ambig), .|, pat[this_di]))
            #end
        end
    end

    res = format_DF(res, dist, guides)
    return res
end