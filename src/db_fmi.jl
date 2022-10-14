# TODO PathTemplate should contain alignment
# TODO alignment has to be correctly written to the detail
function search_chrom(
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
    for (i, fwd_offt_i) in enumerate(fwd_offt)
        fwd_pos_filter = Set{Int64}([])
        rve_pos_filter = Set{Int64}([])

        curr_dist = 0
        for offt in fwd_offt_i
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
                rve_iter = rve_iter .- length(offt.seq) .+ 1
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


function search_fmiDB(
    guides::Vector{LongDNA{4}}, mpt::PathTemplates, motif::Motif, fmidbdir::String,
    detail::String; distance::Int = 2)

    if distance > 2
        throw("Max distance is 2, anything more than that is impractically slow.")
    end

    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    guides_ = copy(guides)

    # we input guides that are in forward search configuration
    fwd_offt = map(x -> templates_to_sequences(x, mpt, motif; dist = distance), guides_)
    ThreadsX.map(ch -> search_chrom(ch, dirname(detail), guides_, motif, fwd_offt, fmidbdir), gi.chrom)
    cleanup_detail(detail)
    return detail
end