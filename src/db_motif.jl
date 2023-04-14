struct MotifPos
    chroms::Vector{<: Unsigned}
    pos::Vector{<: Unsigned}
    isplus::BitVector
    sequences::Vector{LongDNA{4}} # large and potentially uneccessary suffixes
    ug::BitVector # bitvector encodes which guides positions share the same sequence!
    ug_count::Int
    bits::BitMatrix # columns are guides, rows are kmers
end


struct MotifDB
    dbi::DBInfo
    prefixes::Set{LongDNA{4}}
    kmer_size::Int
    kmers::Dict{LongDNA{4}, Int}
end


function Base.convert(::Type{BitVector}, x::Vector{LociRange})
    x_ = BitVector()
    bitflip = true
    for xi in x
        if bitflip
            append!(x_, BitVector(ones(length(xi))))
        else
            append!(x_, BitVector(zeros(length(xi))))
        end
        bitflip = !bitflip
    end
    return x_
end


function bits_to_counts(bits::BitVector, len::Int)
    x = zeros(Int, len)
    xi = 1
    bitflip = true
    i = 0
    for bi in bits
        if bi == bitflip
            i += 1
        else
            x[xi] = i
            bitflip = !bitflip
            xi += 1
            i = 1
        end
    end
    x[xi] = i
    return x
end


function get_ug_ranges(bv::BitVector, idx::Int)
    val = bv[idx]
    len = length(bv)
    start = idx
    while start != 1 && val == bv[start - 1]
        start -= 1
    end
    stop = idx
    while stop != len && val == bv[stop + 1]
        stop += 1
    end
    return start:stop
end


function gatherofftargets(
    dbi::DBInfo,
    chrom::K,
    chrom_name::String,
    is_antisense::Bool) where {K<:BioSequence}
 
    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    chrom_name_ = convert(dbi.gi.chrom_type, findfirst(isequal(chrom_name), dbi.gi.chrom))

    if length(dbi.motif) != 0
        guides_pos = findguides(dbi, chrom, is_antisense)
        guides = ThreadsX.map(x -> removepam(chrom[x], pam_loci), guides_pos)
        guides = add_extension(guides, guides_pos, dbi, chrom, is_antisense)
        guides, guides_pos = normalize_to_PAMseqEXT(guides, guides_pos, dbi, is_antisense)
        guides_pos = convert.(dbi.gi.pos_type, guides_pos)
        return (chrom_name_, guides, guides_pos, !is_antisense)
    end
    return nothing
end


"""
```
build_motifDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    prefix_len::Int = 7;
    skipmer_size::Int = Int(floor(length_noPAM(motif) / (motif.distance + 3))))
```

Prepare motifDB index for future searches using `search_motifDB`.

Will return a path to the database location, same as `storage_dir`.
When this database is used for the guide off-target scan it is similar 
to `search_linearDB`, however additional filter is applied on top of 
prefix filtering. Suffixes are used for next filter, similarly to 
pigeon hole principle - depending on the size of the skipmer `skipmer_size`.
For example, Cas9 off-target within distance 4 (d) might be 20bp long.
We skip `prefix_len` of 7, and are left with 13bp which can be split into 3 
skipmer (r) of size 4, 1bp will be left unused. However when searching within 
distance of 4 and for prefix where initial alignment was of distance 3 (m) and
adjustment parameter is 0 (a). We are obliged to find at least **r - (d - m + a)** 
which is **3 - (4 - 3 + 0) = 2** this many skipmers inside the off-targets. 

There exist also another approach which builds on the idea that it might be more efficient
to find at least two kmers of smaller size (named 01*0 seed) rather than one larger kmer 
(pigeon hole principle). You can use the `adjust` option for that during `search_linearDB` step.

Be sure to understand implications of using `motifDB` as using wrong parameters 
on `skipmer_size` might result in leaky filtering in relation to the assumed 
distance `dist` and adjustment `adjust` during search step in `search_motifDB`.

# Arguments

`name` - Your preferred name for this index for easier identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif defines what kind of gRNA to search for.

`storage_dir`  - Folder path to the where index will be saved and many prefix files.

`prefix_len`  - Size of the prefix by which off-targets are indexed. Prefix of 8 or larger will be the fastest,
                however it will also result in large number of files. 

`skipmer_size` - Size of the skipmer as described above. Be careful when setting this too large!

# Examples
```julia-repl
$(make_example_doc("motifDB"))
```
"""
function build_motifDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    prefix_len::Int = 7;
    skipmer_size::Int = Int(floor(length_noPAM(motif) / (motif.distance + 3))))

    dbi = ARTEMIS.DBInfo(genomepath, name, motif)

    # step 1
    @info "Step 1: Searching chromosomes."
    # For each chromosome parallelized we build database
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
    # Don't parallelize here as you can likely run out of memory (chromosomes are large)
    prefixes = Base.mapreduce(
        x -> do_linear_chrom(x, getchromseq(dbi.gi.is_fa, reader[x]), dbi, prefix_len, storage_dir), 
        union,
        dbi.gi.chrom)
    close(ref)
    prefixes = Set(prefixes)

    # step 2
    @info "Step 2: Constructing per prefix db."
    # Iterate over all prefixes and merge different chromosomes
    kmers = all_kmers(skipmer_size)
    kmers = Dict(zip(kmers, 1:length(kmers)))
    @showprogress 60 for prefix in prefixes # can be parallelized here ?! memory?!
        guides = Vector{LongDNA{4}}()
        loci = Vector{Loc}()
        for chrom in dbi.gi.chrom
            p = joinpath(storage_dir, string(prefix), string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                append!(guides, pdb.suffix)
                append!(loci, pdb.loci)
            end
        end
        rm(joinpath(storage_dir, string(prefix)), recursive = true)
        (guides, loci_range, loci) = unique_guides(guides, loci)
        guides_ = ThreadsX.map(x -> as_bitvector_of_kmers(x, kmers), guides) # only suffix part
        guides_ = hcat(guides_...) # as BitMatrix - columns are guides, rows are kmers

        loci_chrom = ThreadsX.map(x -> x.chrom, loci)
        loci_pos = ThreadsX.map(x -> x.pos, loci)
        loci_isplus = BitVector(ThreadsX.map(x -> x.isplus, loci))

        ug_count = length(loci_range)
        loci_range = convert(BitVector, loci_range)
        save(
            MotifPos(
                loci_chrom, loci_pos, loci_isplus, guides,
                loci_range, ug_count, guides_),
            joinpath(storage_dir, string(prefix) * ".bin"))
    end

    save(MotifDB(dbi, prefixes, skipmer_size, kmers), joinpath(storage_dir, "motifDB.bin"))
    @info "Finished constructing motifDB in " * storage_dir
    return storage_dir
end


function guide_to_bitvector(guide::Vector{LongDNA{4}}, bits::BitMatrix, kmers::Dict{LongDNA{4}, Int}, min_count::Int)
    in_guide = bits[kmers[guide[1]], :]
    for i in 2:length(guide)
        in_guide += bits[kmers[guide[i]], :]
    end
    return findall(in_guide .>= min_count)
end


function search_prefix(
    prefix::LongDNA{4},
    dist::Int,
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNA{4}},
    gskipmers::Vector{Vector{LongDNA{4}}},
    kmers::Dict{LongDNA{4}, Int},
    adjust::Int,
    storage_dir::String)

    # prefix alignment against all the guides
    prefix_len = length(prefix)
    suffix_len = length_noPAM(dbi.motif) + dbi.motif.distance - prefix_len
    prefix_aln = Base.map(g -> prefix_align(g, prefix, suffix_len, dist), guides)
    isfinal = Base.map(x -> x.isfinal, prefix_aln)

    if all(isfinal)
        return
    end

    detail_path = joinpath(detail, "detail_" * string(prefix) * ".csv")
    detail_file = open(detail_path, "w")

    # what is the alignment score so far - how many distance is accumulated
    row_min = Base.map(x -> minimum(x.v[prefix_len - 1, :]), prefix_aln)

    # if any of the guides requires further alignment
    sdb = load(joinpath(storage_dir, string(prefix) * ".bin"))
    sdb_counts = bits_to_counts(sdb.ug, sdb.ug_count)
    for i in 1:length(guides)
        if !isfinal[i]
            # at maximum can be length(gskipmers[i])
            # but that means
            # some distance can be eaten by prefix -row_min[i]
            # maximum allowed distance is dist
            min_count = length(gskipmers[i]) - (dist - row_min[i] + adjust)
            g_bits = guide_to_bitvector(gskipmers[i], sdb.bits, kmers, min_count)
            if length(g_bits) > 0
                offtargets = sdb.sequences[g_bits]
                for (j, suffix) in enumerate(offtargets)
                    suffix_aln = suffix_align(suffix, prefix_aln[i])

                    if suffix_aln.dist <= dist

                        if dbi.motif.extends5
                            guide_stranded = reverse(prefix_aln[i].guide)
                            aln_guide = reverse(suffix_aln.guide)
                            aln_ref = reverse(suffix_aln.ref)
                        else
                            guide_stranded = prefix_aln[i].guide
                            aln_guide = suffix_aln.guide
                            aln_ref = suffix_aln.ref
                        end
                        noloc = string(guide_stranded) * "," * aln_guide * "," * 
                                aln_ref * "," * string(suffix_aln.dist) * ","

                        for o in get_ug_ranges(sdb.ug, sum(sdb_counts[1:g_bits[j]]))
                            strand = "-"
                            if sdb.isplus[o]
                                strand = "+"
                            end
                            write(
                                detail_file, 
                                noloc * 
                                dbi.gi.chrom[sdb.chroms[o]] * "," * 
                                string(sdb.pos[o]) * "," * 
                                strand * "\n")
                        end
                    end
                end
            end
        end
    end

    close(detail_file)
    return 
end


"""
```
search_motifDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3,
    adjust::Int = 0)
```

Find all off-targets for `guides` within distance of `distance` using motifDB located at `storage_dir`.

Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

```
5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM
    
3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide
```

# Arguments

`storage_dir` - Directory containing motifDB.

`guides` - a vector of gRNAs without PAM.

`output_file` - File to which write detailed results. This search will create intermediate 
files which will have same name as output_file, but with a sequence prefix. Final file
will contain all those intermediate files, and other files with sequence prefix will be deleted.

`distance` - Defines maximum levenshtein distance (insertions, deletions, mismatches) for 
which off-targets are considered.  

`adjust` - This will be crucial parameter for tightening second layer of filtering after, 
the initial prefix alignment. For example, Cas9 off-target within distance 4 (d) might be 20bp long.
We skip `prefix_len` of 7, and are left with 13bp which can be split into 3 skipmers (r) of size 4, 
1bp will be left unused. However when searching within distance of 4 and for prefix where initial 
alignment was of distance 3 (m) and adjustment parameter is 0 (a). We are obliged to find at least 
`r - (d - m + a)` which is `3 - (4 - 3 + 0) = 2` this many skipmers inside the off-targets. 


# Examples
```julia-repl
$(make_example_doc("motifDB"))
```
"""
function search_motifDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3,
    adjust::Int = 0)

    sdb = load(joinpath(storage_dir, "motifDB.bin"))
    if distance > sdb.dbi.motif.distance
        error("Maximum distance is " * string(sdb.dbi.motif.distance))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    prefix_len = length(first(sdb.prefixes))
    gskipmers = ThreadsX.map(x -> collect(Set(as_skipkmers(x[(prefix_len + 1):end], sdb.kmer_size))), guides_)
    ThreadsX.map(p -> search_prefix(
        p, distance, sdb.dbi, dirname(output_file), guides_, gskipmers, sdb.kmers, adjust, storage_dir), sdb.prefixes)
    
    cleanup_detail(output_file)
    return
end