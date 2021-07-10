"
Final SuffixDB unit that contains all guides from
all chromosomes that start with the `prefix` and their locations.
"
struct SuffixDB
    prefix::LongSequence{DNAAlphabet{4}}
    suffix::Vector{LongSequence{DNAAlphabet{4}}}
    suffix_loci_idx::Vector{LociRange}
    loci::Vector{Loc}
end


function unique_guides(
    guides::Vector{LongSequence{DNAAlphabet{4}}}, 
    loci::Vector{Loc})
    order = sortperm(guides)
    guides = guides[order]
    loci = loci[order]
    (guides, loci_range) = ranges(guides)
    return (guides, loci_range, loci)
end


struct LinearDB
    dbi::DBInfo
    prefixes::Set{LongSequence{DNAAlphabet{4}}}
end


"""
`build_linearDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String,
    prefix_len::Int = 7)`

Build a DB of offtargets for the given `motif`,
DB groups off-targets by their prefixes.

Will return a path to the database location, same as `storagedir`.
When this database is used for the guide off-target scan it is similar 
to linear in performance, hence the name.

There is an optimization that if the alignment becomes imposible against
the prefix we don't search the off-targets grouped inside the prefix.
Therefore it is advantageous to select larger prefix than maximum 
search distance, however in that case number of files also grows.
"""
function build_linearDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String,
    prefix_len::Int = 7)

    if prefix_len <= motif.distance
        throw("prefix_len $prefix_len is <= " * string(motif.distance))
    end

    dbi = DBInfo(genomepath, name, motif)

    # step 1
    @info "Step 1: Searching chromosomes."
    # For each chromsome paralelized we build database
    ref = open(dbi.filepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = dbi.filepath * ".fai") : TwoBit.Reader(ref)
    # Don't paralelize here as you can likely run out of memory (chromosomes are large)
    prefixes = Base.map(x -> do_linear_chrom(x, getchromseq(dbi.is_fa, reader[x]), dbi, prefix_len, storagedir), dbi.chrom)
    close(ref)

    prefixes = Set(vcat(prefixes...))

    # step 2
    @info "Step 2: Constructing per prefix db."
    # Iterate over all prefixes and merge different chromosomes
    ThreadsX.map(prefixes) do prefix
        guides = Vector{LongDNASeq}()
        loci = Vector{Loc}()
        for chrom in dbi.chrom
            p = joinpath(storagedir, string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                append!(guides, pdb.suffix)
                append!(loci, pdb.loci)
                rm(p)
            end
        end
        (guides, loci_range, loci) = unique_guides(guides, loci)
        sdb = SuffixDB(prefix, guides, loci_range, loci)
        save(sdb, joinpath(storagedir, string(prefix) * ".bin"))
    end

    linDB = LinearDB(dbi, prefixes)
    save(linDB, joinpath(storagedir, "linearDB.bin"))
    @info "Finished constructing linearDB in " * storagedir
    return storagedir
end


function search_prefix(
    prefix::LongSequence{DNAAlphabet{4}},
    dist::Int,
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNASeq},
    storagedir::String)

    res = zeros(Int, length(guides), dist + 1)
    if detail != ""
        detail_path = joinpath(detail, "detail_" * string(prefix) * ".csv")
        detail_file = open(detail_path, "w")
    end

    # prefix alignment against all the guides
    suffix_len = length_noPAM(dbi.motif) + dbi.motif.distance - length(prefix)
    prefix_aln = Base.map(g -> prefix_align(g, prefix, suffix_len, dist), guides)
    isfinal = Base.map(x -> x.isfinal, prefix_aln)

    if all(isfinal)
        return res
    end

    # if any of the guides requires further alignment 
    # load the SuffixDB and iterate
    sdb = load(joinpath(storagedir, string(prefix) * ".bin"))
    for i in 1:length(guides)
        if !isfinal[i]
            for (j, suffix) in enumerate(sdb.suffix)
                suffix_aln = suffix_align(suffix, prefix_aln[i])
                if suffix_aln.dist <= dist
                    sl_idx = sdb.suffix_loci_idx[j]
                    res[i, suffix_aln.dist + 1] += length(sl_idx)
                    if detail != ""
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
    end

    if detail != ""
        close(detail_file)
    end
    return res
end


"""
`search_linearDB(storagedir::String, guides::Vector{LongDNASeq}, dist::Int = 4; detail::String = "")`

Will search the previously build database for the off-targets of the `guides`. 
Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

`
5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM

3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide
`

# Arguments

`dist` - defines maximum levenshtein distance (insertions, deletions, mismatches) for
which off-targets are considered.  
`detail` - path and name for the output file. This search will create intermediate 
files which will have same name as detail, but with a sequence prefix. Final file
will contain all those intermediate files. Leave `detail` empty if you are only 
interested in off-target counts returned by the linearDB.  
"""
function search_linearDB(storagedir::String, guides::Vector{LongDNASeq}, dist::Int = 4; detail::String = "")
    ldb = load(joinpath(storagedir, "linearDB.bin"))
    prefixes = collect(ldb.prefixes)
    if dist > length(first(prefixes)) || dist > ldb.dbi.motif.distance
        error("For this database maximum distance is " * 
              string(min(ldb.dbi.motif.distance, length(first(prefixes)))))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if ldb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    #res = zeros(Int, length(guides_), dist + 1)
    res = ThreadsX.mapreduce(p -> search_prefix(p, dist, ldb.dbi, dirname(detail), guides_, storagedir), +, prefixes)
    #for p in prefixes
    #    res += search_prefix(p, dist, ldb.dbi, dirname(detail), guides_, storagedir)
    #end
    
    if detail != ""
        # clean up detail files into one file
        open(detail, "w") do detail_file
            write(detail_file, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")
            for prefix in filter(x -> occursin("detail_", x), readdir(dirname(detail)))
                prefix_file = joinpath(dirname(detail), prefix)
                open(prefix_file, "r") do prefix_file
                    for ln in eachline(prefix_file)
                        write(detail_file, ln * "\n")
                    end
                end
                rm(prefix_file)
            end
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end