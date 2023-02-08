"
Final SuffixDB unit that contains all guides from
all chromosomes that start with the `prefix` and their locations.
"
struct SuffixDB
    prefix::LongDNA{4}
    suffix::Vector{LongDNA{4}}
    suffix_loci_idx::Vector{LociRange}
    loci::Vector{Loc}
end


function unique_guides(
    guides::Vector{LongDNA{4}}, 
    loci::Vector{Loc})
    order = sortperm(guides)
    guides = guides[order]
    loci = loci[order]
    (guides, loci_range) = ranges(guides)
    return (guides, loci_range, loci)
end


struct LinearDB
    dbi::DBInfo
    prefixes::Set{LongDNA{4}}
end


"""
```
build_linearDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    prefix_len::Int = 7)
```

Prepare linearDB index for future searches using `search_linearDB`.

Will return a path to the database location, same as `storage_dir`.
When this database is used for the guide off-target scan it is similar 
to linear in performance, hence the name. There is an optimization that 
if the alignment becomes impossible against
the prefix we don't search the off-targets grouped inside the prefix.
Therefore, it is advantageous to select much larger prefix than maximum 
search distance, however in that case number of files also grows. For example,
if interested with searches within distance 4, preferably use prefix length of 
7 or 8.

# Arguments
`name` - Your preferred name for this index for easier identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif defines what kind of gRNA to search for.

`storage_dir`  - Folder path to the where index will be saved with name `linearDB.bin` and many prefix files.

`prefix_len`  - Size of the prefix by which off-targets are indexed. Prefix of 8 or larger will be the fastest,
                however it will also result in large number of files.

# Examples
```julia-repl
$(make_example_doc())
```
"""
function build_linearDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    prefix_len::Int = 7)

    if prefix_len <= motif.distance
        throw("prefix_len $prefix_len is <= " * string(motif.distance))
    end

    dbi = DBInfo(genomepath, name, motif)

    # step 1
    @info "Step 1: Searching chromosomes."
    # For each chromsome paralelized we build database
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
    # Don't paralelize here as you can likely run out of memory (chromosomes are large)
    prefixes = Base.map(x -> do_linear_chrom(x, getchromseq(dbi.gi.is_fa, reader[x]), dbi, prefix_len, storage_dir), dbi.gi.chrom)
    close(ref)

    prefixes = Set(vcat(prefixes...))

    # step 2
    @info "Step 2: Constructing per prefix db."
    # Iterate over all prefixes and merge different chromosomes
    ThreadsX.map(prefixes) do prefix
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
        sdb = SuffixDB(prefix, guides, loci_range, loci)
        save(sdb, joinpath(storage_dir, string(prefix) * ".bin"))
    end

    linDB = LinearDB(dbi, prefixes)
    save(linDB, joinpath(storage_dir, "linearDB.bin"))
    @info "Finished constructing linearDB in " * storage_dir
    return storage_dir
end


function search_prefix(
    prefix::LongDNA{4},
    dist::Int,
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNA{4}},
    storage_dir::String)

    # prefix alignment against all the guides
    suffix_len = length_noPAM(dbi.motif) + dbi.motif.distance - length(prefix)
    prefix_aln = Base.map(g -> prefix_align(g, prefix, suffix_len, dist), guides)
    isfinal = Base.map(x -> x.isfinal, prefix_aln)

    if all(isfinal)
        return
    end

    detail_path = joinpath(detail, "detail_" * string(prefix) * ".csv")
    detail_file = open(detail_path, "w")

    # if any of the guides requires further alignment 
    # load the SuffixDB and iterate
    sdb = load(joinpath(storage_dir, string(prefix) * ".bin"))
    for i in 1:length(guides)
        if !isfinal[i]
            for (j, suffix) in enumerate(sdb.suffix)
                suffix_aln = suffix_align(suffix, prefix_aln[i])
                if suffix_aln.dist <= dist
                    sl_idx = sdb.suffix_loci_idx[j]
                    offtargets = sdb.loci[sl_idx.start:sl_idx.stop]
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
                    for offt in offtargets
                        write(detail_file, noloc * decode(offt, dbi) * "\n")
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
search_linearDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 4)
```

Find all off-targets for `guides` within distance of `dist` using linearDB located at `storage_dir`.

Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

```
5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM
    
3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide
```

# Arguments

`output_file` - Path and name for the output file, this will be comma separated table, therefore `.csv` extension is preferred. 
This search will create intermediate files which will have same name as `output_file`, but with a sequence prefix. Final file
will contain all those intermediate files.

`distance` - Defines maximum levenshtein distance (insertions, deletions, mismatches) for 
which off-targets are considered.

# Examples
```julia-repl
$(make_example_doc())
```
"""
function search_linearDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 4)

    ldb = load(joinpath(storage_dir, "linearDB.bin"))
    prefixes = collect(ldb.prefixes)
    if distance > length(first(prefixes)) || distance > ldb.dbi.motif.distance
        error("For this database maximum distance is " * 
              string(min(ldb.dbi.motif.distance, length(first(prefixes)))))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if ldb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    ThreadsX.map(p -> search_prefix(p, distance, ldb.dbi, dirname(output_file), guides_, storage_dir), prefixes)
    
    cleanup_detail(output_file)
    return
end