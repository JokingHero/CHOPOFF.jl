struct LociRange
    start::UInt32
    stop::UInt32
end

import Base.length
@inline length(x::LociRange) = x.stop - x.start + 1

"
Temporary structure that is build per chromosome, per prefix.
Slow to save/load.
"
struct PrefixDB
    prefix::LongSequence{DNAAlphabet{4}}
    suffix::Dict{LongSequence{DNAAlphabet{4}}, Vector{Loc}}
end


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


function to_suffix(prefix::LongSequence{DNAAlphabet{4}}, d::Dict)
    guides = collect(keys(d))
    loci = collect(values(d))
    lengths = length.(loci)
    cs = cumsum(lengths)
    loci = reduce(vcat, loci)
    stops = cs
    starts = stops .- lengths .+ 1
    loci_range = LociRange.(starts, stops)
    return SuffixDB(prefix, guides, loci_range, loci)
end


struct LinearDB
    dbi::DBInfo
    prefixes::Set{LongSequence{DNAAlphabet{4}}}
end


" 
Extract 3' extension of the guide:
TTN ... EXT
CCN ... EXT
"
function getExt3(chrom::K, chrom_max::Int, ext_start::Int, dist::Int) where K <: BioSequence
    ext_end = ext_start + dist
    if ext_start > chrom_max
        ext = LongDNASeq(repeat("-", dist))
    elseif ext_end > chrom_max
        ext = chrom[ext_start:chrom_max] 
        ext = ext * LongDNASeq(repeat("-", dist - length(ext)))
    else
        ext = chrom[ext_start:ext_end]
    end
    return ext
end

" 
Extract 5' extension of the guide:
EXT ... NAA
EXT ... NGG
"
function getExt5(chrom::K, ext_end::Int, dist::Int)  where K <: BioSequence
    ext_start = ext_end - dist + 1
    if ext_end < 1
        ext = LongDNASeq(repeat("-", dist))
    elseif ext_start < 1
        ext = chrom[1:ext_end]
        ext = LongDNASeq(repeat("-", dist - length(ext))) * ext
    else
        ext = chrom[ext_start:ext_end]
    end
    return ext
end


"
Prefix is the length of our common part between guides!
Guide reversal happens here too, so that alignment is more
convienient.
"
function pushguides!(
    output::Vector{PrefixDB},
    dbi::DBInfo,
    chrom::K,
    chrom_name::String,
    reverse_comp::Bool,
    prefix_len::Int) where {K<:BioSequence}

    query = reverse_comp ? dbi.motif.rve : dbi.motif.fwd 
    pam_loci = reverse_comp ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd

    g_len = length(query) - length(pam_loci)
    suffix_len = g_len - prefix_len

    if length(query) != 0
        chrom_max = lastindex(chrom)

        for x in findall(query, chrom)
            guide = LongDNASeq(chrom[x])
            guide = removepam(guide, pam_loci)
            # add extension
            if dbi.motif.extends5 && reverse_comp
                # CCN ... EXT
                guide = guide * getExt3(chrom, chrom_max, last(x) + 1, dbi.motif.distance)
                guide = complement(guide)
                pos_ = first(x)
                # becomes GGN ... EXT
            elseif dbi.motif.extends5 && !reverse_comp 
                # EXT ... NGG
                guide = getExt5(chrom, first(x) - 1, dbi.motif.distance) * guide
                guide = reverse(guide)
                pos_ = last(x)
                # becomes GGN ... EXT
            elseif !dbi.motif.extends5 && reverse_comp 
                # EXT ... NAA
                guide = getExt5(chrom, first(x) - 1, dbi.motif.distance) * guide
                guide = reverse_complement(guide)
                pos_ = last(x)
                # becomes TTA ... EXT
            else #!dbi.motif.extends5 && !reverse_comp
                # TTN ... EXT
                guide = guide * getExt3(chrom, chrom_max, last(x) + 1, dbi.motif.distance)
                pos_ = first(x)
            end

            gprefix = guide[1:prefix_len]
            gsuffix = guide[prefix_len+1:end]
            chrom_name_ = convert(dbi.chrom_type, findfirst(isequal(chrom_name), dbi.chrom))
            pos_ = convert(dbi.pos_type, pos_)
            loc = Loc(chrom_name_, pos_, !reverse_comp)

            prefix_idx = findfirst(x -> gprefix == x.prefix, output)

            if prefix_idx !== nothing
                if haskey(output[prefix_idx].suffix, gsuffix)
                    push!(output[prefix_idx].suffix[gsuffix], loc)
                else
                    output[prefix_idx].suffix[gsuffix] = [loc]
                end
            else
                push!(output, PrefixDB(gprefix, Dict(gsuffix => [loc])))
            end
        end
    end
    return output
end


function do_linear_chrom(chrom_name::String, chrom::K, dbi::DBInfo, prefix_len::Int, storagedir::String) where K<:BioSequence
    @info "Working on $chrom_name"
    output = Vector{PrefixDB}()
    pushguides!(output, dbi, chrom, chrom_name, false, prefix_len)
    pushguides!(output, dbi, chrom, chrom_name, true, prefix_len)
    # save small files as it is to large to store them in memory
    for pdb in output
        save(pdb, joinpath(storagedir, string(pdb.prefix) * "_" * chrom_name * ".bin"))
    end
    return [x.prefix for x in output]
end


function getseq(isfa, record)
    return isfa ? FASTA.sequence(record) : TwoBit.sequence(record)
end


"
Build a DB of offtargets for the given `motif`,
DB groups off-targets by their prefixes.

Will return a path to the database location, same as `storagedir`.
When this database is used for the guide off-target scan it is similar 
to linear in performance, hence the name.

There is an optimization that if the alignment becomes imposible against
the prefix we don't search the off-targets grouped inside the prefix.
Therefore it is advantageous to select larger prefix than maximum 
search distance, however in that case number of files also grows.
"
function buildlinearDB(
    name::String,
    genomepath::String,
    motif::Motif,
    prefix_len::Int,
    storagedir::String)

    dbi = DBInfo(genomepath, name, motif)

    # step 1
    @info "Step 1: Searching chromosomes."
    # For each chromsome paralelized we build database
    ref = open(dbi.filepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = dbi.filepath * ".fai") : TwoBit.Reader(ref)
    prefixes = Base.map(x -> do_linear_chrom(x, getseq(dbi.is_fa, reader[x]), dbi, prefix_len, storagedir), dbi.chrom)
    close(ref)

    prefixes = Set(vcat(prefixes...))

    # step 2
    @info "Step 2: Constructing per prefix db."
    # Iterate over all prefixes and merge different chromosomes
    for prefix in prefixes
        merged = Dict()
        for chrom in dbi.chrom
            p = joinpath(storagedir, string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                merge!(vcat, merged, pdb.suffix)
                rm(p)
            end
        end
        sdb = to_suffix(prefix, merged)
        save(sdb, joinpath(storagedir, string(prefix) * ".bin"))
    end

    linDB = LinearDB(dbi, prefixes)
    save(linDB, joinpath(storagedir, "linearDB.bin"))
    @info "Finished building linearDB in " * storagedir
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
    # CGCATG-CCATCACAGCAAGG - guide
    # CGCATGACCATCAtAGCAAtG AGG - ref
    #=if string(prefix) == "GTAACGA"
        @info "Working on $prefix"
        println(decode(Loc{UInt8,UInt32}(0x08, 0x0000440c, true), dbi))
    else
        return res
    end =#

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
    #@info "$sdb"

    for (i, g) in enumerate(guides)
        if !isfinal[i]
            for (j, suffix) in enumerate(sdb.suffix)
                suffix_aln = suffix_align(suffix, prefix_aln[i])
                if suffix_aln.dist <= dist
                    sl_idx = sdb.suffix_loci_idx[j]
                    res[i, suffix_aln.dist + 1] += length(sl_idx)
                    if detail != ""
                        offtargets = sdb.loci[sl_idx.start:sl_idx.stop]
                        #@info "$offtargets"
                        #show(sdb.loci)
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
                            #@info "$offt"
                            write(detail_file, noloc * decode(offt, dbi) * "\n")
                        end
                    end
                end
            end
        end
    end

    if detail != ""
        close(detail_file)
        @info "closing $prefix"
    end
    return res
end


"
Will search the previously build database for the off-targets of the `guides`. 
Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM

3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide

`dist` - defines maximum levenshtein distance (insertions, deletions, mismatches) for
which off-targets are considered.
`detail` - path and name for the output file. This search will create intermediate 
files which will have same name as detail, but with a sequence prefix. Final file
will contain all those intermediate files. Leave `detail` empty if you are only 
interested in off-target counts returned by the searchDB.
"
function searchlinearDB(storagedir::String, dist::Int, guides::Vector{LongDNASeq}; detail::String = "")
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
    for p in prefixes
        res += search_prefix(p, dist, ldb.dbi, dirname(detail), guides_, storagedir)
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

    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end