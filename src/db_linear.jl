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
    prefix::DNAMer
    suffix::Dict{DNAMer, Vector{Loc}}
end


"
Final SuffixDB unit that contains all guides from
all chromosomes that start with the `prefix` and their locations.
"
struct SuffixDB
    prefix::DNAMer
    suffix::Vector{DNAMer}
    suffix_loci_idx::Vector{LociRange}
    loci::Vector{Loc}
end


function to_suffix(prefix::DNAMer, d::Dict)
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
    prefixes::Set{DNAMer}
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

    g_len = length(query) - length(collect(vcat(pam_loci...)))
    suffix_len = g_len - prefix_len

    if length(query) != 0
        for x in findall(query, chrom)
            guide = LongDNASeq(chrom[x])
            guide = removepam(guide, pam_loci)
            pos_ = dbi.motif.extends5 ? last(x) : first(x)
            if reverse_comp
                guide = reverse_complement(guide)
                pos_ = dbi.motif.extends5 ? first(x) : last(x)
            end
            # all the guides are in dist*N + 20*N + *NGG
            # we want PAM-guide style for alignment
            if dbi.motif.extends5
                guide = reverse(guide)
            end

            gprefix = DNAMer{prefix_len}(guide[1:prefix_len])
            gsuffix = DNAMer{suffix_len}(guide[prefix_len+1:end])
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
                merge!(merged, pdb.suffix)
                rm(p)
            end
        end
        sdb = to_suffix(prefix, merged)
        save(sdb, joinpath(storagedir, string(prefix) * ".bin"))
    end

    linDB = LinearDB(dbi, prefixes)
    save(linDB, joinpath(storagedir, "linearDB.bin"))
    @info "Finished building linearDB in " * storagedir
    return nothing
end


function search_prefix(
    prefix::DNAMer,
    dist::Int,
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNASeq},
    storagedir::String)

    #@info "Working on $prefix"
    res = zeros(Int, length(guides), dist + 1)

    # step 1 produce prefix alignment against all the guides
    prefix_aln = Base.map(g -> prefix_levenshtein(g, prefix, dist), guides)
    isfinal = Base.map(x -> x.isfinal, prefix_aln)
    if prefix == mer"CCTGAAG"
        println(prefix_aln) 
        println("G: ", guides)
        println("Pref: ", prefix)
    end

    if all(isfinal)
        return res
    end

    # step 2, if any of the guides requires further alignment 
    #         load the SuffixDB and iterate
    sdb = load(joinpath(storagedir, string(prefix) * ".bin"))

    for (i, g) in enumerate(guides)
        if !isfinal[i]
            for (j, suffix) in enumerate(sdb.suffix)
                if suffix == mer"ACTCTTTCGACCGCCCC"
                    println("g: ", g)
                    println("suf: ", suffix)
                    println("pa: ", prefix_aln[i])
                end
                suffix_aln = suffix_levenshtein(g, suffix, prefix_aln[i], dist)
                if suffix_aln <= dist
                    sl_idx = sdb.suffix_loci_idx[j]
                    res[i, suffix_aln + 1] += length(sl_idx)
                    #println("Loc:", sdb.loci[sl_idx.start:sl_idx.stop])
                    #println("sa: ", suffix_aln) 
                    #if detail != ""
                        # write to file
                        # dbi to decode(Loc)
                    #end
                end
            end
        end
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
    if dist > length(first(prefixes))
        error("For this database maximum distance is " * string(length(first(prefixes))))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if ldb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = Base.map(x -> search_prefix(x, dist, ldb.dbi, detail, guides_, storagedir), prefixes)
    res = reduce(+, res)
    
    if detail != ""
        # clean up detail files into one file
    end

    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end