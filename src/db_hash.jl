struct HashDB{T<:Unsigned, K<:Union{UInt8, UInt16, UInt32}}
    dbi::DBInfo
    mtp::PathTemplates
    bins::Vector{BinaryFuseFilter{K}}
    counts::Vector{T}
    ambig::Union{AmbigIdx, Nothing}
end


function get_count_idx(bins::Vector{BinaryFuseFilter{K}}, guide::UInt64, right::Bool) where K<:Union{UInt8, UInt16, UInt32} 
    direction = right ? (1:length(bins)) : (length(bins):-1:1)
    for i in direction
        if guide in bins[i]
            return i
        end
    end
    return nothing
end


function guides_to_bins(
    guides::Vector{UInt64}, seed::UInt64, max_iterations::Int, max_count::Int; precision::DataType = UInt16)
    max_count_type = smallestutype(unsigned(max_count))
    max_count = convert(max_count_type, max_count)
    guides = sort(guides)
    guides, counts = ranges(guides)
    counts = convert.(max_count_type, min.(length.(counts), max_count))

    bins = Vector{BinaryFuseFilter{precision}}()
    bins_counts = Vector{max_count_type}()
    for count in 1:max_count
        idx = ThreadsX.findall(x -> x == count, counts)
        if length(idx) != 0
            push!(bins_counts, convert(max_count_type, count))
            push!(bins, BinaryFuseFilter{precision}(guides[idx]; seed = seed, max_iterations = max_iterations))
        end
    end
    # now order from smallest to largest
    order = sortperm(bins_counts)
    bins = bins[order]
    bins_counts = bins_counts[order]

    # use full db to estimate error rate!
    err_r = 0
    err_l = 0
    for (guide, real_count) in zip(guides, counts)
        est_count_right = get_count_idx(bins, guide, true)
        if isnothing(est_count_right) 
            throw("All guides should be inside!")
        else
            est_count_right = bins_counts[est_count_right]
        end
        est_count_left = get_count_idx(bins, guide, false)
        if isnothing(est_count_left) 
            throw("All guides should be inside!")
        else
            est_count_left = bins_counts[est_count_left]
        end
        err_r += Int(real_count != est_count_right)
        err_l += Int(real_count != est_count_left)
    end
    return bins, bins_counts, 
        length(err_l) / length(guides), 
        length(err_r) / length(guides)
end


"""
```
build_hashDB(
    name::String, 
    genomepath::String, 
    motif::Motif;
    storage_path::String = "",
    seed::UInt64 = UInt64(0x726b2b9d438b9d4d),
    max_iterations::Int = 10,
    max_count::Int = 10,
    precision::DataType = UInt16)
```

Prepare hashDB index for future searches using `search_hashDB`.


# Arguments
`name` - Your prefered name for this index to ease future identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif defines what kind of gRNA to search for.

`storage_path`  - Path to the where index will be saved.

`seed`  - Optional. Seed is used during hashing for randomization.

`max_iterations` - When finding hashing structure for binary fuse filter it might fail sometimes, 
                   we will retry `max_iterations` number of times though.
             
`max_count`  - Above this count we put all unique off-targets into one bin. 
               Put number here that is the minimum number of off-targets that you think is fine
               for the distance of 1.
                
`precision`- The higher the precision the larger the database, but also chances for error decrease dramatically.
             We support UInt8, UInt16, and UInt32.


# Examples
```julia-repl
# prepare libs
using ARTEMIS, BioSequences

# use ARTEMIS example genome
genome = joinpath(
    vcat(
        splitpath(dirname(pathof(ARTEMIS)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))

# finally, build a hashDB
build_hashDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0))
```
"""
function build_hashDB(
    name::String, 
    genomepath::String, 
    motif::Motif;
    storage_path::String = "",
    seed::UInt64 = UInt64(0x726b2b9d438b9d4d),
    max_iterations::Int = 10,
    max_count::Int = 10,
    precision::DataType = UInt16)

    dbi = DBInfo(genomepath, name, motif)
    if motif.distance != 1
        @info "Distance enforced to 1."
        motif = setdist(motif, 1)
    end
    @info "Building Motif templates..."
    mtp = build_PathTemplates(length_noPAM(motif), motif.distance)
    
    # gather all unique off-targets
    guides = Vector{UInt64}()
    ambig = gatherofftargets!(guides, dbi)
    ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing

    # guides are here of length 21
    bins, counts, err_left, err_right = 
        guides_to_bins(guides, seed, max_iterations, max_count; precision = precision)
    db = HashDB(dbi, mtp, bins, counts, ambig)

    if storage_path != ""
        save(db, storage_path)
        @info "Finished constructing hashDB in " * storage_path * " consuming "  * 
        string(round((filesize(joinpath(storage_path, "hashDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    end

    @info "Estimated probability of miscounting an elements in the bins is: " * 
        "\nRight: " * string(round(err_right; digits = 6)) *
        "\nLeft: " * string(round(err_left; digits = 6))
    return db
end


"Find which guides are matched with the `seq`."
function findbits(seq::LongDNA{4}, cols::Vector{BitVector})
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


"""
```
search_hashDB(
    db::HashDB,
    guides::Vector{LongDNA{4}},
    right::Bool)
```

Estimate off-target counts for `guides` using hashDB stored at `storage_dir`.

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

# use ARTEMIS example genome
ARTEMIS_path = splitpath(dirname(pathof(ARTEMIS)))[1:end-1]
genome = joinpath(vcat(ARTEMIS_path, "test", "sample_data", "genome", "semirandom.fa"))

# build a hashDB
db = build_hashDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0))

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(ARTEMIS_path, "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)

# finally, get results!
hdb_res = search_hashDB(db, guides, false)
```
"""
function search_hashDB(
    db::HashDB,
    guides::Vector{LongDNA{4}},
    right::Bool)

    if any(isambig.(guides)) # TODO I think we support it now - should not matter much
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    guides_ = copy(guides)
    len_noPAM_noEXT = length_noPAM(db.dbi.motif)
    len = len_noPAM_noEXT + db.dbi.motif.distance

    if any(len_noPAM_noEXT .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(db.dbi.motif))
    end
    # reverse guides so that PAM is always on the left
    if db.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = zeros(Int, length(guides_), 2)
    for (i, s) in enumerate(guides_)

        pat = templates_to_sequences_extended(s, db.mtp; dist = 1)
        d0 = collect(pat[1])
        d1 = collect(pat[2])

        for d0_s in convert.(UInt64, d0)
            idx = get_count_idx(db.bins, d0_s, right)
            if !isnothing(idx)
                res[i, 1] += db.counts[idx]
            end
        end

        for d1_s in convert.(UInt64, d1)
            idx = get_count_idx(db.bins, d1_s, right)
            if !isnothing(idx)
                res[i, 2] += db.counts[idx]
            end
        end

        if !isnothing(db.ambig)
            bits_mapped = Base.map(x -> findbits(x, db.ambig), d0)
            res[i, 1] += sum(reduce(.|, bits_mapped))
            bits_mapped = Base.map(x -> findbits(x, db.ambig), d1)
            res[i, 2] +=  sum(reduce(.|, bits_mapped))
        end
    end

    res = format_DF(res, 1, guides)
    return res
end