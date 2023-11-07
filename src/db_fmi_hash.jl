
function search_chrom2(
    chrom::String, 
    detail::String, 
    guides::Vector{LongDNA{4}},
    motif::Motif,
    fwd_offt::Vector{Vector{Path}},
    storage_dir::String)

    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    #seq = getchromseq(gi.is_fa, reader[chrom])
    fmi = load(joinpath(storage_dir, chrom * ".bin"))
    
    detail_path = joinpath(detail, "detail_" * string(chrom) * ".csv")
    detail_file = open(detail_path, "w")

    # wroking on this guide and his all possible off-targets
    for (i, fwd_offt_i) in enumerate(fwd_offt) # for each guide

        for offt in fwd_offt_i # for each potential OT

            fwd_iter = locate(offt.seq, fmi)
            if motif.extends5 
                fwd_iter = fwd_iter .+ length(offt.seq) .- 1
            end
            for pos in fwd_iter # TODO this might be slow too
                line = string(guides[i]) * "," * "no_alignment" * "," * 
                    string(offt.seq) * "," * string(offt.dist) * "," *
                    chrom * "," * string(pos) * "," * "+" * "\n"
                write(detail_file, line)
            end
                
            rve_iter = locate(reverse_complement(offt.seq), fmi)
            if !motif.extends5
                rve_iter = rve_iter .+ length(offt.seq) .- 1
            end
            for pos in rve_iter
                line = string(guides[i]) * "," * "no_alignment" * "," * 
                    string(offt.seq) * "," * string(offt.dist) * "," *
                    chrom * "," * string(pos) * "," * "-" * "\n"
                write(detail_file, line)
            end
        end
    end
    close(detail_file)
    #close(ref)
    return 
end



"""
```
search_fmiDB(
    guides::Vector{LongDNA{4}}, 
    mpt::PathTemplates, 
    motif::Motif, 
    fmidbdir::String,
    output_file::String; 
    distance::Int = 2)
```

Search FM-index for off-targets using hash based brute-force enumeration method.

**Experimental! Proof-of-concept!**

This method uses `PathTemplates` build on top of `Motif` to enumerate all
possible off-target sequences, next these sequences are filtered for presence inside the 
genome using hash based methods. Finally, leftover sequences are found in the genome
using FM-index.

# Arguments
`guides` - a vector of gRNAs without PAM.

`mpt` - PathTemplates object that contains abstraction for all possible alignments

`motif` - Motif defines what kind of gRNA to search for. Has to be compatible with `mpt`.

`fmidbdir`   - Path to the folder where FM-index was build using `build_fmi`.

`bff` - Hashed version of all offtargets - Binary Fuse Filter. Has to be compatible with the motif.

`output_file`  - Where output will be saved.

`distance`  - Search distance.


# Examples
```julia-repl
# prepare libs
using ARTEMIS, BioSequences

distance = 3

# make a temporary directory
tdir = tempname()
fmi_dir = joinpath(tdir, "fmi")
mkpath(fmi_dir)

# use ARTEMIS example genome
artemis_path = splitpath(dirname(pathof(ARTEMIS)))[1:end-1]
genome = joinpath(vcat(artemis_path, 
    "test", "sample_data", "genome", "semirandom.fa"))
# build FM-index
build_fmiDB(genome, fmi_dir)
motif = Motif("Cas9"; distance = distance)
mpt = build_PathTemplates(motif)

# build hashDB
db = build_hashDB("samirandom", genome, motif)

# prepare output folder
res_dir = joinpath(tdir, "results")
mkpath(res_dir)

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(artemis_path, 
    "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)
    
# finally, make results!
res_path = joinpath(res_dir, "results.csv")
search_fmiDB_hash(guides, mpt, motif, fmi_dir, db, res_path, distance)

# load results
using DataFrames, CSV
res = DataFrame(CSV.File(res_path))

# filter results by close proximity
res = filter_overlapping(res, 20 + distance)

# summarize results into a table of counts by distance
summary = summarize_offtargets(res, distance)
```
"""
function search_fmiDB_hash(
    guides::Vector{LongDNA{4}}, mpt::PathTemplates, motif::Motif, fmidbdir::String, db::HashDB,
    output_file::String, distance::Int; right::Bool = true)

    if any(isambig.(guides)) # TODO I think we support it now - should not matter much
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    guides_ = copy(guides)
    len_noPAM_noEXT = length_noPAM(motif)
    len = len_noPAM_noEXT + motif.distance

    if any(len_noPAM_noEXT .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(db.dbi.motif))
    end
    # reverse guides so that PAM is always on the left
    if db.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    fwd_offt = Vector{Vector{Path}}()
    for (i, s) in enumerate(guides_)
        
        # TODO - we can have sequences of specific length, but not produce them for larger distances, e.g. seaerrching for d2 with max dist 3 - not make dist 3 patterns
        pat = ARTEMIS.templates_to_sequences(s, mpt; dist = distance) # still super slow for larger distances

        pat_in_genome = falses(lastindex(pat))
        for (o, ot) in enumerate(pat)
            subs = expand_ambiguous(LongDNA{4}(ot.seq) * repeat(dna"N", len - length(ot.seq)))
            for sub in subs 
                if !isnothing(ARTEMIS.get_count_idx(db.bins, convert(UInt64, sub), right))
                    pat_in_genome[o] = true
                    continue # we skip the checks as one of the subsequences was found in the genome
                end
            end
        end

        @info "Fraction of seqeunces that passed " * string(ceil(sum(pat_in_genome)/lastindex(pat_in_genome); digits = 4))
        push!(fwd_offt, pat[pat_in_genome]) # TODO probably replace with preallocated 
    end

    # we input guides that are in forward search configuration
    ThreadsX.map(ch -> search_chrom2(ch, dirname(output_file), guides_, motif, fwd_offt, fmidbdir), gi.chrom)
    cleanup_detail(output_file)
    return 
end