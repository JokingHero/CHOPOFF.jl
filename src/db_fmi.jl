# TODO PathTemplate contains alignment
function search_chrom(
    chrom::String, 
    detail::String, 
    guides::Vector{LongDNA{4}},
    fwd_offt::Vector{Vector{Path}},
    storagedir::String)

    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    #seq = ARTEMIS.getchromseq(gi.is_fa, reader[chrom])
    
    res = zeros(Int, length(guides), fwd_offt[end][end].dist + 1)
    fmi = load(joinpath(storagedir, chrom * ".bin"))
    
    if detail != ""
        detail_path = joinpath(detail, "detail_" * string(chrom) * ".csv")
        detail_file = open(detail_path, "w")
    end

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

            fwd_iter = ARTEMIS.locate(offt.seq, fmi)
            for pos in fwd_iter
                if !(pos in fwd_pos_filter)
                    push!(fwd_pos_filter, pos)
                    res[i, offt.dist + 1] += 1
                    if detail != ""
                        line = string(guides[i]) * "," * "no_alignment" * "," * 
                            string(offt.seq) * "," * string(offt.dist) * "," *
                            chrom * "," * string(pos) * "," * "+" * "\n"
                        write(detail_file, line)
                    end
                end
            end
                
            rve_iter = ARTEMIS.locate(reverse_complement(offt.seq), fmi)
            for pos in rve_iter
                if !(pos in rve_pos_filter)
                    push!(rve_pos_filter, pos)
                    res[i, offt.dist + 1] += 1
                    if detail != ""
                        line = string(guides_[i]) * "," * "no_alignment" * "," * 
                            string(offt.seq) * "," * string(offt.dist) * "," *
                            chrom * "," * string(pos) * "," * "-" * "\n"
                        write(detail_file, line)
                    end
                end
            end
        end
    end
    if detail != ""
        close(detail_file)
    end
    #close(ref)
    return res
end


function search_fmiDB_patterns(
    fmidbdir::String, mpt::PathTemplates,
    guides::Vector{LongDNA{4}}; detail::String = "", distance::Int = 2)

    if distance > 2
        throw("Max distance is 2, anything more than that is impractically slow.")
    end

    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    guides_ = copy(guides)

    # we input guides that are in forward search configuration
    fwd_offt = map(x -> templates_to_sequences(x, mpt, motif; dist = distance), guides_)
    res = ThreadsX.mapreduce(ch -> search_chrom(ch, dirname(detail), guides_, fwd_offt, storagedir), +, gi.chrom)

    if detail != ""
        # clean up detail files into one file
        cleanup_detail(detail)
    end

    res = format_DF(res, distance, guides)
    return res
end


# NOT FINISHED?! - untested
function as_partial_alignments(s::String, motif::Motif, len::Int = 10, alphabet::Vector{Char} = ['A', 'C', 'T', 'G'])
    s = s[1:(len + motif.distance)]
    if motif.distance == 0
        return [s]
    end
    comb = comb_of_d1(s, alphabet)
    for i in 1:(motif.distance-1)
        comb = foldxt(union, Map(x -> comb_of_d1(x, alphabet)), comb)
    end

    comb = collect(Set(map(x -> x[1:len], collect(comb))))

    if motif.extends5
        pam = reverse(motif.fwd[motif.pam_loci_fwd])
        comb .= pam .* comb
        comb_rev = complement.(comb)
        comb = reverse.(comb)
    else
        pam = motif.fwd[motif.pam_loci_fwd]
        comb .= pam .* comb
        comb_rev = reverse_complement.(comb)
    end
    comb = vcat(expand_ambiguous.(comb)...)
    comb_rev = vcat(expand_ambiguous.(comb_rev)...)
    return (comb, comb_rev)
end


# NOT FINISHED!!!
function search_fmiDB_raw(
    fmidbdir::String, genomepath::String, motif::Motif,
    guides::Vector{LongDNA{4}}; detail::String = "", prefix_length::Int = 10)

    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if motif.extends5
        guides_ = reverse.(guides_)
    end

    partials = map(x -> as_partial_alignments(String(x), motif, prefix_length), guides_)
    
    res = zeros(Int, length(guides_), motif.distance + 1)
    ref = open(genomepath, "r")
    reader = gi.is_fa ? FASTA.Reader(ref, index = genomepath * ".fai") : TwoBit.Reader(ref)

    for (ic, chrom) in enumerate(gi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        seq = getchromseq(gi.is_fa, reader[chrom])
    
        for (idx, gp) in enumerate(partials)
            p, p_rev = gp
            pos = mapreduce(x -> locateall(x, fmi), union, p)
            dists = ThreadsX.map(x -> levenshtein(
                guides_[idx], 
                seq[(pos - prefix_length - motif.distance - 1):(pos + prefix_length)], 
                motif.distance), 
                pos)

            for d in 0:dist
                res[idx, d + 1] += sum(dists .== d)
            end

            
            pos = mapreduce(x -> locateall(x, fmi), union, p_rev)
            dists = ThreadsX.map(x -> levenshtein(
                guides_[idx], 
                seq[(pos - prefix_length - motif.distance - 1):(pos + prefix_length)], 
                motif.distance), 
                pos)

            for d in 0:dist
                res[idx, d + 1] += sum(dists .== d)
            end
        end
    end
    close(ref)
    
    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end


# https://mikael-salson.univ-lille.fr/talks/VST_010seed.pdf
# TODO
# this takes only closest PAM upstream of the two seq - suboptimal? optimal?
function search_pamDB(
    fmidbdir::String, genomepath::String, pamdbpath::String,
    guides::Vector{LongDNA{4}}; detail::String = "", distance::Int = 3)

    pamDB = ARTEMIS.load(pamdbpath)
    if distance > pamDB.motif.distance
        throw("DIstance is too large for selected PathTemplates. Max distance is" * 
            string(pamDB.motif.distance))
    end

    # 01*0 seed always exists when pattern is partitioned in at least distance + 2 parts
    adj_seed_size = Int(floor(length(guides[1])/(distance + 2)))
    parts_number = Int(floor(length(guides[1])/adj_seed_size))

    gi = ARTEMIS.load(joinpath(fmidbdir, "genomeInfo.bin"))
    is_fa = ARTEMIS.is_fasta(genomepath)

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    # TODO CPf1 style...
    if pamDB.motif.extends5
        guides_ = reverse.(guides_)
    end
    # this can't be unique and skipmers are ambig handleable
    guides_skipmers = Base.map(x -> unique(ARTEMIS.as_skipkmers(x, adj_seed_size)), guides_)

    res = zeros(Int, length(guides_), distance + 1)
    ref = open(genomepath, "r")
    reader = is_fa ? FASTA.Reader(ref, index = genomepath * ".fai") : TwoBit.Reader(ref)

    offt_len = ARTEMIS.length_noPAM(pamDB.motif) + distance
    # 1       23
    # hit---hitPAM
    # |     |  
    scan_dist = offt_len - adj_seed_size
    # worst case
    # hithit---PAM
    #       |  | = 23 - adj_seed_size * 2 = 15
    # hit-hit--PAM
    # |   |
    # diff
    #        | | 
    # 23 - diff - adj_seed_size 
    # adj_seed_size is the difference here 
    # hit---hitPAM
    #          | = 0
    for (ic, chrom) in enumerate(gi.chrom)
        fmi = ARTEMIS.load(joinpath(fmidbdir, chrom * ".bin"))
        seq = ARTEMIS.getchromseq(is_fa, reader[chrom])
        pam_loc_fwd = pamDB.pam_loc_fwd[ic]
        pam_loc_fwd_len = length(pam_loc_fwd)
        pam_loc_rve = pamDB.pam_loc_rve[ic]
            
        for (i, gs) in enumerate(guides_skipmers)
            # FORWARD
            gs_fwd = reverse.(gs)
            # we lose knowledge on order and which one was which so this is no good!!!
            hits = mapreduce(x -> ARTEMIS.locateall(x, fmi), vcat, gs_fwd)
            hits = unique(hits)
            sort!(hits)
            diff_hits = diff(hits)
            two_are_close = diff_hits .<= scan_dist
            diff_hits = diff_hits[two_are_close]
            pushfirst!(two_are_close, 0) # we take the most right element
            hits = hits[two_are_close]
            # hitsPAM
            # |-> |
            hits .= hits .+ adj_seed_size
            check_offtarget_pam = zeros(Int, length(hits))
            check_offtarget = falses(length(hits))
            for (ih, h) in enumerate(hits)
                idx = searchsortedfirst(pam_loc_fwd, h) # find first >= h
                if idx <= pam_loc_fwd_len
                    # |    |-PAM
                    # 14 462 469
                    if (pam_loc_fwd[idx] - h) <= (scan_dist - diff_hits[ih])
                        check_offtarget[ih] = true
                        check_offtarget_pam[ih] = pam_loc_fwd[idx]
                    end
                end
            end
            check_offtarget_pam = unique(check_offtarget_pam[check_offtarget])
            check_offtarget_pam .= check_offtarget_pam .- 1
            offt_seq = map(x -> seq[(x - offt_len + 1):x], check_offtarget_pam)
            if pamDB.motif.extends5
                offt_seq = reverse.(offt_seq)
            end
            align_score = map(x -> ARTEMIS.levenshtein(guides_[i], x, distance), offt_seq)
            for d in 0:distance
                res[i, d + 1] += sum(align_score .== d)
            end

            # REVERSE
            gs_rve = complement.(gs)
            hits = mapreduce(x -> ARTEMIS.locateall(x, fmi), vcat, gs_rve)
            hits = unique(hits)
            sort!(hits)
            diff_hits = diff(hits)
            two_are_close = diff_hits .<= scan_dist
            diff_hits = diff_hits[two_are_close]
            push!(two_are_close, 0) # we take the most left element
            hits = hits[two_are_close]
            hits .= hits .- 1

            check_offtarget_pam = zeros(Int, length(hits))
            check_offtarget = falses(length(hits))
            for (ih, h) in enumerate(hits)
                idx = searchsortedlast(pam_loc_rve, h) # find last value <= h
                if idx != 0 
                    if (h - pam_loc_rve[idx]) <= (scan_dist - diff_hits[ih])
                        check_offtarget[ih] = true
                        check_offtarget_pam[ih] = pam_loc_rve[idx]
                    end
                end
            end
            check_offtarget_pam = unique(check_offtarget_pam[check_offtarget])
            check_offtarget_pam .= check_offtarget_pam .+ 1
            offt_seq = map(x -> seq[x:(x + offt_len - 1)], check_offtarget_pam)
            align_score = map(x -> ARTEMIS.levenshtein(guides_[i], complement(x), distance), offt_seq)
            for d in 0:distance
                res[i, d + 1] += sum(align_score .== d)
            end
        end
    end
    close(ref)

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:distance]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
end