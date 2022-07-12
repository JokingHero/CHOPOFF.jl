struct HashDB{T<:Unsigned, K<:Union{UInt8, UInt16, UInt32}}
    dbi::DBInfo
    bins_d0::Vector{BinaryFuseFilter{K}}
    counts_d0::Vector{T}
    bins_d1::Vector{BinaryFuseFilter{K}}
    counts_d1::Vector{T}
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
    error_right = Vector{Int}()
    error_left = Vector{Int}()
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

        if  real_count != est_count_right
            push!(error_right, Int(est_count_right) - Int(real_count))
        end

        if  real_count != est_count_left
            push!(error_left, Int(est_count_left) - Int(real_count))
        end
    end
    error_rate_left = length(error_left) / length(guides)
    error_rate_right = length(error_right) / length(guides)
    return bins, bins_counts, error_rate_left, error_rate_right
end


"""
```
build_hashDB(
    name::String,
    genomepath::String, 
    motif::Motif,
    storagedir::String;
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

`storagedir`  - Folder path to the where index will be saved with name `hashDB.bin`.

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
# make a temporary directory
tdir = tempname()
hdb_path = joinpath(tdir, "hashDB")
mkpath(hdb_path)

# use CRISPRofftargetHunter example genome
genome = joinpath(
    vcat(
        splitpath(dirname(pathof(CRISPRofftargetHunter)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))

# finally, build a hashDB
build_hashDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0), 
    hdb_path)
```
"""
function build_hashDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String;
    seed::UInt64 = UInt64(0x726b2b9d438b9d4d),
    max_iterations::Int = 10,
    max_count::Int = 10,
    precision::DataType = UInt16)

    dbi = DBInfo(genomepath, name, motif)
    if motif.distance != 1
        @info "Distance enforced to 1."
        motif = setdist(motif, 1)
    end
    
    # gather all unique off-targets
    guides = Vector{UInt64}()
    ambig = gatherofftargets!(guides, dbi)
    ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing

    # guides are here of length 21
    bins_d1, counts_d1, error_d1_right, error_d1_left = 
        guides_to_bins(guides, seed, max_iterations, max_count; precision = precision)

    # now D0
    guides = guides .>> 2 # removes extension
    bins_d0, counts_d0, error_d0_right, error_d0_left = 
        guides_to_bins(guides, seed, max_iterations, max_count; precision = precision)

    db = HashDB(dbi, bins_d0, counts_d0, bins_d1, counts_d1, ambig)
    save(db, joinpath(storagedir, "hashDB.bin"))

    @info "Finished constructing hashDB in " * storagedir * " consuming "  * 
        string(round((filesize(joinpath(storagedir, "hashDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    @info "Estimated probability of miscounting an elements in the bins is: " * 
        "\nRight: D0 " * string(round(error_d0_right; digits = 6)) * 
        " D1 " * string(round(error_d1_right; digits = 6)) * 
        "\nLeft: D0 " * string(round(error_d0_left; digits = 6)) * 
        " D1 " * string(round(error_d1_left; digits = 6))
    return storagedir
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
using CRISPRofftargetHunter, BioSequences

# make a temporary directory
tdir = tempname()
hdb_path = joinpath(tdir, "hashDB")
mkpath(hdb_path)

# use CRISPRofftargetHunter example genome
coh_path = splitpath(dirname(pathof(CRISPRofftargetHunter)))[1:end-1]
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
function search_hashDB(
    storagedir::String,
    guides::Vector{LongDNA{4}},
    right::Bool)

    if any(isambig.(guides))
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    db = load(joinpath(storagedir, "hashDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if db.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), 2)
    for (i, s) in enumerate(guides_)
        # 0 distance
        res[i, 1] += isnothing(db.ambig) ? 0 : sum(findbits(s, db.ambig))
        idx = get_count_idx(db.bins_d0, convert(UInt64, s), right)
        if !isnothing(idx)
            res[i, 1] += db.counts_d0[idx]
        end
        
        norm, border = comb_of_d1(s)
        norm = collect(norm)
        for comb in convert.(UInt64, norm)
            idx = get_count_idx(db.bins_d0, comb, right)
            if !isnothing(idx)
                res[i, 2] += db.counts_d0[idx]
            end
        end
        
        for comb in convert.(UInt64, border)
            idx = get_count_idx(db.bins_d1, comb, right)
            if !isnothing(idx)
                idx0 = get_count_idx(db.bins_d0, comb >> 2, right)
                if !isnothing(idx0)
                    res[i, 2] += db.counts_d1[idx]
                end
            end
        end

        if !isnothing(db.ambig)
            bits_mapped = map(x -> findbits(x, db.ambig), vcat(norm, border))
            res[i, 2] += sum(reduce(.|, bits_mapped))
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:1]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, vcat(col_d, :guide))
    return res
end