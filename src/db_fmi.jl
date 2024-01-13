function search_chrom(
    chrom::String, 
    detail::String, 
    guides::Vector{LongDNA{4}},
    mpt::PathTemplates,
    storage_dir::String)

    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    #seq = getchromseq(gi.is_fa, reader[chrom])
    fmi = load(joinpath(storage_dir, chrom * ".bin"))
    detail_path = joinpath(detail, "detail_" * string(chrom) * ".csv")
    detail_file = open(detail_path, "w")

    guides_fmi = guide_to_template_format.(copy(guides); alphabet = ALPHABET_UINT8)
    guides_fmi_rc = guide_to_template_format.(copy(guides), true; alphabet = ALPHABET_UINT8)

    # wroking on this guide and his all possible off-targets
    for (i, g) in enumerate(guides) # for each guide
        if mpt.motif.extends5
            guide_stranded = reverse(g)
        else
            guide_stranded = g
        end

        ot = guides_fmi[i][mpt.paths] # GGN + 20N + extension
        ot_rc = guides_fmi_rc[i][mpt.paths] # CCN + 20N + extension
        if mpt.motif.extends5
            reverse!(ot; dims = 2) # extension + 20N + NGG
        else # Cpf1
            reverse!(ot_rc; dims = 2) # extension + 20N + NAAA
        end
        
        ot_len = size(ot)[2]

        fwd_pos_filter = Set{Int64}([])
        rve_pos_filter = Set{Int64}([])

        curr_dist = 0
        for i in 1:size(ot)[1] # for each potential OT
            ot_i = @view ot[i, :]
            ot_rc_i = @view ot_rc[i, :]
            # keep track of overlapping positions
            if mpt.distances[i] > curr_dist
                fwd_pos_filter = union(fwd_pos_filter, fwd_pos_filter .+ 1, fwd_pos_filter .-1)
                rve_pos_filter = union(rve_pos_filter, rve_pos_filter .+ 1, rve_pos_filter .-1)
                curr_dist = mpt.distances[i]
            end

            fwd_iter = locate(ot_i, fmi)
            if mpt.motif.extends5 
                fwd_iter = fwd_iter .+ ot_len .- 1
            end

            for pos in fwd_iter
                if !(pos in fwd_pos_filter)
                    push!(fwd_pos_filter, pos)
                    line = string(guide_stranded) * "," * "no_alignment" * "," * 
                        string(LongDNA{4}(reinterpret(DNA, ot_i))) * "," * string(mpt.distances[i]) * "," *
                        chrom * "," * string(pos) * "," * "+" * "\n"
                    write(detail_file, line)
                end
            end
                
            rve_iter = locate(ot_rc_i, fmi)
            if !mpt.motif.extends5
                rve_iter = rve_iter .+ ot_len .- 1
            end
            for pos in rve_iter
                if !(pos in rve_pos_filter)
                    push!(rve_pos_filter, pos)
                    line = string(guide_stranded) * "," * "no_alignment" * "," * 
                    string(LongDNA{4}(reinterpret(DNA, ot_rc_i))) * "," * string(mpt.distances[i]) * "," *
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
mpt = build_PathTemplates(motif; withPAM = true) # its important to add PAM here!

# prepare output folder
res_dir = joinpath(tdir, "results")
mkpath(res_dir)

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(artemis_path, 
    "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)
    
# finally, make results!
res_path = joinpath(res_dir, "results.csv")
search_fmiDB(guides, mpt, fmi_dir, res_path; distance = 1)

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
    guides::Vector{LongDNA{4}}, mpt::PathTemplates, fmidbdir::String,
    output_file::String; distance::Int = 2)

    if distance > 2
        throw("Max distance is 2, anything more than that is impractically slow.")
    end

    if any(length.(guides) .!= length_noPAM(mpt.motif))
        throw("Wrong guide length.")
    end

    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if mpt.motif.extends5
        guides_ = reverse.(guides_)
    end

    mpt = restrictDistance(mpt, distance)
    ThreadsX.map(ch -> search_chrom(ch, dirname(output_file), guides_, mpt, fmidbdir), gi.chrom)
    cleanup_detail(output_file)
    return 
end