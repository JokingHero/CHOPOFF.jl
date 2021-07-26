
"
Final CSuffixDB unit that contains all guides from
all chromosomes that start with the `prefix` and their locations.
"
struct CSuffixDB
    prefix::LongDNASeq
    #vcf expected max 5M

    # these are 1 to 1, guide to Loci
    ambig_suffixes::Vector{LongDNASeq}
    ambig_loci::Vector{Loc}
    
    # the ones that are repeated very often - around 10M max
    overflow_suffixes::Vector{LongDNASeq}
    overflow_loci_range::Vector{LociRange}
    overflow_loci::Vector{Loc}

    suffix1::LongDNASeq
    n_changes::Vector{UInt8}
    changes::Vector{UInt8}
    suffix_loci_count::Vector{UInt8}
    loci::Vector{Loc}
end


struct CompactDB
    dbi::DBInfo
    prefixes::Set{LongDNASeq}
end


const SEQ_DIFF = vec(collect(Iterators.product(1:30, [DNA_A, DNA_C, DNA_G, DNA_T])))
@inline function encode_seq_diff(pos::Int64, base::T) where T <: DNA
    offset = 0
    if base == DNA_C
        offset = 30
    elseif base == DNA_G
        offset = 60
    elseif base == DNA_T
        offset = 90
    end
    return UInt8(pos + offset)
end


@inline function next_suffix!(suffix::LongDNASeq, changes::Vector{UInt8})
    aln_start_idx = 1
    if length(changes) > 0
        aln_start_idx = first(SEQ_DIFF[changes[1]])
    end
    for c in changes
        i, base = SEQ_DIFF[c]
        suffix[i] = base
    end
    return aln_start_idx
end


function guide_diff(p::LongDNASeq, n::LongDNASeq)
    diff = Vector{UInt8}()
    @inbounds for (i, g) in enumerate(p)
        if !ismatch(g, n[i])
            push!(diff, encode_seq_diff(i, n[i]))
        end
    end
    return diff
end


function guide_changes(guides::Vector{LongDNASeq})
    n_changes = zeros(UInt8, length(guides))
    changes = Vector{UInt8}()
    for i in 2:length(guides)
        i_changes = guide_diff(guides[i-1], guides[i])
        n_changes[i] = length(i_changes)
        append!(changes, i_changes)
    end
    return (n_changes, changes)
end


"""
`build_compactDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String,
    prefix_len::Int = 7;
    bed_filter = "",
    vcf = "")`

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
function build_compactDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String,
    prefix_len::Int = 7;
    bed_filter = "",
    vcf = "") # TODO

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
    @info "Step 2: Constructing per prefix compactDB."
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
        # ambig suffixes
        ambig_idx = ThreadsX.findall(x -> n_ambigous(x) > 0, guides)
        ambig_suffixes = guides[ambig_idx]
        ambig_loci = loci[ambig_idx]
        deleteat!(guides, ambig_idx)
        deleteat!(loci, ambig_idx)
        # repeated guides with count >= 255
        (guides, loci_range, loci) = unique_guides(guides, loci)
        overflow_idx = ThreadsX.findall(x -> length(x) >= 255)
        overflow_suffixes = guides[overflow_idx]
        deleteat!(guides, overflow_idx)
        overflow_loci_range = loci_range[overflow_idx]
        deleteat!(loci_range, overflow_idx)
        unrolled_overflow_loci = vcat(overflow_loci_range...)
        overflow_loci = loci[unrolled_overflow_loci]
        deleteat!(loci, unrolled_overflow_loci)
        overflow_loci_range = recalculate_ranges(overflow_loci_range)
        
        # compressed guides depend on the previous guide
        order = order_by_hamming_and_prefix(guides)
        guides = guides[order]
        loci_range = recalculate_ranges(loci_range)
        loci_range = loci_range[order]
        loci = loci[vcat(loci_range...)]
        suffix_loci_count = UInt8.(length.(loci_range))
        (n_changes, changes) = guide_changes(guides)
        csdb = CSuffixDB(
            prefix, ambig_suffixes, ambig_loci,
            overflow_suffixes, overflow_loci_range, overflow_loci,
            guides[1], n_changes, changes, suffix_loci_count, loci)
        save(csdb, joinpath(storagedir, string(prefix) * ".bin"))
    end

    linDB = LinearDB(dbi, prefixes)
    save(linDB, joinpath(storagedir, "linearDB.bin"))
    @info "Finished constructing linearDB in " * storagedir
    return storagedir
end


function write_detail(
    detail_file::IOStream, offtargets::Vector{Loc}, dbi::DBInfo, 
    guide_stranded::LongDNASeq, aln_guide::LongDNASeq, 
    aln_ref::LongDNASeq, aln_dist::Int)
    if dbi.motif.extends5
        guide_stranded = reverse(guide_stranded)
        aln_guide = reverse(aln_guide)
        aln_ref = reverse(aln_ref)
    end
    noloc = string(guide_stranded) * "," * aln_guide * "," * 
            aln_ref * "," * string(aln_dist) * ","
    for offt in offtargets
        write(detail_file, noloc * decode(offt, dbi) * "\n")
    end
    return nothing
end


function search_compact_prefix!(
    res::Matrix{Int64},
    is_early_stopped::Vector{Bool},
    prefix::LongDNASeq,
    dist::Int,
    early_stop::Vector{Int},
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNASeq},
    storagedir::String)

    if detail != ""
        detail_path = joinpath(detail, "detail_" * string(prefix) * ".csv")
        detail_file = open(detail_path, "w")
    end

    # prefix alignment against all the guides
    suffix_len = length_noPAM(dbi.motif) + dbi.motif.distance - length(prefix)
    prefix_aln = Base.map(g -> prefix_align(g, prefix, suffix_len, dist), guides)
    isfinal = Base.map(x -> x.isfinal, prefix_aln)

    if all(isfinal)
        return nothing
    end

    # if any of the guides requires further alignment 
    # load the SuffixDB and iterate
    db = load(joinpath(storagedir, string(prefix) * ".bin"))
    # first check highly repeated guides
    for (j, suffix) in enumerate(db.overflow_suffixes)
        for i in 1:length(guides)
            if !isfinal[i] && !is_early_stopped[i]
                suffix_aln = suffix_align(suffix, prefix_aln[i])
                if suffix_aln.dist <= dist
                    sl_idx = db.overflow_suffix_loci_range[j]
                    res[i, suffix_aln.dist + 1] += length(sl_idx)
                    if res[i, suffix_aln.dist + 1] >= early_stop[suffix_aln.dist + 1]
                        is_early_stopped[i] <- true
                    end
                    if detail != ""
                        write_detail(
                            detail_file, db.overflow_loci[sl_idx], dbi, 
                            prefix_aln[i].guide, suffix_aln.guide, 
                            suffix_aln.ref, suffix_aln.dist)
                    end
                end
            end
        end
    end

    # now check the rest of the guides
    suffix = db.suffix1
    idx_changes_start = 1
    idx_loci_start = 1
    for (j, n) in enumerate(db.n_changes)
        aln_start_idx = next_suffix!(suffix, db.changes[idx_changes_start:(idx_changes_start + n - 1)])
        idx_changes_start += n
        
        for i in 1:length(guides)
            if !isfinal[i] && !is_early_stopped[i]
                suffix_aln = suffix_align!(suffix, prefix_aln[i], aln_start_idx)
                if suffix_aln.dist <= dist
                    res[i, suffix_aln.dist + 1] += db.suffix_loci_count[j]
                    if res[i, suffix_aln.dist + 1] >= early_stop[suffix_aln.dist + 1]
                        is_early_stopped[i] <- true
                    end
                    if detail != ""
                        write_detail(
                            detail_file, 
                            db.loci[idx_loci_start:(idx_loci_start + db.suffix_loci_count[j] - 1)], 
                            dbi, prefix_aln[i].guide, suffix_aln.guide,
                            suffix_aln.ref, suffix_aln.dist)
                    end
                end
            end
        end

        idx_loci_start += db.suffix_loci_count[j]
    end

    # TODO vcf guides

    if detail != ""
        close(detail_file)
    end
    return nothing
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
function search_compactDB(
    storagedir::String, 
    guides::Vector{LongDNASeq}, 
    dist::Int = 4,
    early_stop::Vector{Int} = zeros(Int, dist); 
    detail::String = "")

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

    res = zeros(Int, length(guides_), dist + 1)
    is_early_stopped = ones(Bool, length(guides))
    for p in prefixes
        search_compact_prefix!(
            res, is_early_stopped, p, dist, early_stop, 
            ldb.dbi, dirname(detail), guides_, storagedir)
    end
    
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
    res.early_stopped = is_early_stopped
    sort!(res, [order(Symbol("early_stopped"), rev = true), col_d])
    return res
end