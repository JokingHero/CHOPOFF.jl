function search_chrom(
    chrom::String, 
    detail::String, 
    guides::Vector{LongDNA{4}},
    motif::Motif,
    mpt::PathTemplates,
    storage_dir::String)

    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    #seq = getchromseq(gi.is_fa, reader[chrom])
    fmi = load(joinpath(storage_dir, chrom * ".bin"))
    detail_path = joinpath(detail, "detail_" * string(chrom) * ".csv")
    detail_file = open(detail_path, "w")

    guides_uint64 = guide_to_template_format.(copy(guides))
    guides_uint64_rc = guide_to_template_format.(reverse_complement.(copy(guides)))
    guides_fmi = guide_to_template_format.(copy(guides_); alphabet = ALPHABET_UINT8)
    guides_fmi_rc = guide_to_template_format.(reverse_complement.(copy(guides_)); alphabet = ALPHABET_UINT8)

    # wroking on this guide and his all possible off-targets
    for (i, fwd_offt_i) in enumerate(fwd_offt) # for each guide
        fwd_pos_filter = Set{Int64}([])
        rve_pos_filter = Set{Int64}([])

        curr_dist = 0
        for offt in fwd_offt_i # for each potential OT
            # keep track of overlapping positions
            if offt.dist > curr_dist
                fwd_pos_filter = union(fwd_pos_filter, fwd_pos_filter .+ 1, fwd_pos_filter .-1)
                rve_pos_filter = union(rve_pos_filter, rve_pos_filter .+ 1, rve_pos_filter .-1)
                curr_dist = offt.dist
            end

            fwd_iter = locate(offt.seq, fmi)
            if motif.extends5 
                fwd_iter = fwd_iter .+ length(offt.seq) .- 1
            end

            for pos in fwd_iter
                if !(pos in fwd_pos_filter)
                    push!(fwd_pos_filter, pos)
                    
                    line = string(guides[i]) * "," * "no_alignment" * "," * 
                        string(offt.seq) * "," * string(offt.dist) * "," *
                        chrom * "," * string(pos) * "," * "+" * "\n"
                    write(detail_file, line)
                end
            end
                
            rve_iter = locate(reverse_complement(offt.seq), fmi)
            if !motif.extends5
                rve_iter = rve_iter .+ length(offt.seq) .- 1
            end
            for pos in rve_iter
                if !(pos in rve_pos_filter)
                    push!(rve_pos_filter, pos)
                    line = string(guides[i]) * "," * "no_alignment" * "," * 
                        string(offt.seq) * "," * string(offt.dist) * "," *
                        chrom * "," * string(pos) * "," * "-" * "\n"
                    write(detail_file, line)
                end
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

Search FM-index for off-targets using brute-force enumeration method.

**Experimental! Proof-of-concept!**

This method uses `PathTemplates` build on top of `Motif` to enumerate all
possible off-target sequences, next, these sequences are found in the genome
using FM-index. This method is impractically slow above distance of 2.

# Arguments
`guides` - a vector of gRNAs without PAM.

`mpt` - PathTemplates object that contains abstraction for all possible alignments

`motif` - Motif defines what kind of gRNA to search for. Has to be compatible with `mpt`.

`fmidbdir`   - Path to the folder where FM-index was build using `build_fmi`.

`output_file`  - Where output will be saved.

`distance`  - Search distance, maximum of 2 is practical.


# Examples
```julia-repl
# prepare libs
using ARTEMIS, BioSequences

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
motif = Motif("Cas9"; distance = 1)
mpt = build_PathTemplates(motif)

# prepare output folder
res_dir = joinpath(tdir, "results")
mkpath(res_dir)

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(artemis_path, 
    "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)
    
# finally, make results!
res_path = joinpath(res_dir, "results.csv")
search_fmiDB(guides, mpt, motif, fmi_dir, res_path; distance = 1)

# load results
using DataFrames, CSV
res = DataFrame(CSV.File(res_path))

# filter results by close proximity
res = filter_overlapping(res, 23)

# summarize results into a table of counts by distance
summary = summarize_offtargets(res, 1)
```
"""
function search_fmiDB(
    guides::Vector{LongDNA{4}}, mpt::PathTemplates, motif::Motif, fmidbdir::String,
    output_file::String; distance::Int = 2)

    if distance > 2
        throw("Max distance is 2, anything more than that is impractically slow.")
    end

    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    guides_ = copy(guides)


    if length(guide) != mpt.len
        throw("Wrong guide length.")
    end
    mpt = restrictDistance(mpt, distance)

    # we input guides that are in forward search configuration
    ThreadsX.map(ch -> search_chrom(ch, dirname(output_file), guides_, motif, mpt, fmidbdir), gi.chrom)
    cleanup_detail(output_file)
    return 
end