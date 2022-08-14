

function build_fmiDB(
    genomepath::String,
    storagedir::String)

    gi = GenomeInfo(genomepath)
    ref = open(gi.filepath, "r")
    reader = gi.is_fa ? FASTA.Reader(ref, index = gi.filepath * ".fai") : TwoBit.Reader(ref)
    @showprogress 60 for chrom in gi.chrom
        fmi = FMIndex(getchromseq(gi.is_fa, reader[chrom]), 16, r = 32)
        p = joinpath(storagedir, chrom * ".bin")
        save(fmi, p)
    end
    close(ref)
    save(gi, joinpath(storagedir, "genomeInfo.bin"))
    @info "FMindex is build."
    return storagedir
end


# assumes pattern_pos is sorted
function is_in_range(
    all_pos::Vector{UInt32}, # contains all pos
    pos::BitVector, # should be as long as all_pos - indicates which pos should be evaluated
    pattern_pos::Vector{Int}, # positions of the kmer
    reverse_comp::Bool,
    extends5::Bool,
    dist::Int, # maximal distance we consider here
    shift::Int, # shift from the begining 
    pam_length::Int,
    kmer_length::Int)
    configuration = (extends5 & reverse_comp) | (!extends5 & !reverse_comp)
    pattern_last_insert_idx = length(pattern_pos) + 1
    pos = copy(pos)
    for (i, b) in enumerate(pos)
        if b # true means we should evaluate the potential
            # check if is close enough
            is_close = false
            if configuration
                # CCN 1 2 3 4 5 6 7 8 9 ... EXT
                # p
                #     p + pam_length
                #     p + shift - 1
                #     T G A C # first >= 
                all_pos_i = all_pos[i] + pam_length + shift - 1
                idx = searchsortedfirst(pattern_pos, all_pos_i)
                if idx != pattern_last_insert_idx
                    is_close = (pattern_pos[idx] - all_pos_i) <= dist
                end
            else
                # EXT ... 1 2 3 4 5 6 7 8 9 N G G
                #                               p
                #                         p - pam_length (3)
                #                         p - shift (1) + 1
                #                   p - kmer_length (4) + 1 # we get first position of the kmer too
                all_pos_i = all_pos[i] - pam_length - shift - kmer_length + 2
                # first value that pattern_pos[idx] >= all_pos[i]
                idx = searchsortedfirst(pattern_pos, all_pos_i)
                if idx != pattern_last_insert_idx
                    is_close = (pattern_pos[idx] - all_pos_i) <= dist
                end
            end
            pos[i] = is_close
        end
    end
    return pos
end


# will update long too (expand part)
function shrink_and_expand!(long::BitVector, counts::Vector{Int})
    short = BitVector(zeros(length(counts)))
    start = 1
    for (i, c) in enumerate(counts)
        stop = start + c - 1
        short[i] = any(long[start:stop])
        if short[i]
            long[start:stop] = ones(Bool, c)
        end
        start = stop + 1
    end
    return short
end


"
Funny enough this can easily support ambig guides!

THIS FUNCTION IS BUGGED somewhere :( GL figuring this out
"
function search_fmiDB(
    fmidbdir::String, motifdbdir::String,
    guides::Vector{LongDNA{4}}, dist::Int = 4; detail::String = "")

    mp = load(joinpath(motifdbdir, "motifDB.bin"))
    mp_counts = bits_to_counts(mp.ug, mp.ug_count)

    if dist > mp.dbi.motif.distance
        error("Maximum distance is " * string(mp.dbi.motif.distance))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if mp.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    pam_length = length(mp.dbi.motif) - length_noPAM(mp.dbi.motif) 
    kmer_length = minkmersize(length_noPAM(mp.dbi.motif), dist)
    gkmers = as_kmers.(guides_, kmer_length)
    if detail != ""
        # clean up detail files into one file
        detail_file = open(detail, "w")
        write(detail_file, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")
    end

    res = zeros(Int, length(guides_), dist + 1)
    # early stopping + write detail
    for (ic, chrom) in enumerate(mp.dbi.gi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        pos_chrom = mp.chroms .== ic
        pos_chrom_fwd = pos_chrom .& mp.isplus # 1 have to be chcked
        pos_chrom_rev = pos_chrom .& (.!mp.isplus)
        for (idx, gk) in enumerate(gkmers)
            # keep track of all pos that have been checked already
            # exclude these from further analyses!
            was_checked = BitVector(zeros(length(mp.pos)))
            for (i , kmer) in enumerate(gk)
                if mp.dbi.motif.extends5
                    kmer_rev = complement(kmer)
                    kmer = reverse(kmer)
                else
                    kmer_rev = copy(kmer)
                    kmer = reverse_complement(kmer)
                end
                kmer_pos = sort(locateall(kmer, fmi))
                pos_chrom_fwd_in_range = is_in_range(
                    mp.pos, (.!was_checked) .& pos_chrom_fwd, 
                    kmer_pos, false, mp.dbi.motif.extends5, dist, i, pam_length, kmer_length)

                pos_seqnames = shrink_and_expand!(pos_chrom_fwd_in_range, mp_counts)
                was_checked .|= pos_chrom_fwd_in_range

                dists = ThreadsX.map(x -> levenshtein(guides_[idx], x, dist), mp.sequences[pos_seqnames])
                mp_counts_g = mp_counts[pos_seqnames]
                for d in 0:dist
                    res[idx, d + 1] += sum(mp_counts_g[dists .== d])
                end

                # reverse now
                kmer_pos = sort(locateall(kmer_rev, fmi))
                pos_chrom_rve_in_range = is_in_range(
                    mp.pos, (.!was_checked) .& pos_chrom_rev, 
                    kmer_pos, true, mp.dbi.motif.extends5, dist, i, pam_length, kmer_length)

                pos_seqnames = shrink_and_expand!(pos_chrom_rve_in_range, mp_counts)
                was_checked .|= pos_chrom_rve_in_range

                dists = ThreadsX.map(x -> levenshtein(guides_[idx], x, dist), mp.sequences[pos_seqnames])
                mp_counts_g = mp_counts[pos_seqnames]
                for d in 0:dist
                    res[idx, d + 1] += sum(mp_counts_g[dists .== d])
                end
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


# Working - but pattern, reduceability is to be improved - currently it is not precise!
function search_fmiDB_patterns(
    fmidbdir::String, genomepath::String, mpt::MotifPathTemplates,
    guides::Vector{LongDNA{4}}; detail::String = "", distance::Int = 3)

    if distance > mpt.motif.distance
        throw("DIstance is too large for selected MotifPathTemplates. Max distance is" * 
            string(mpt.motif.distance))
    end

    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    pam = mpt.motif.fwd[mpt.motif.pam_loci_fwd]
    if mpt.motif.extends5
        guides_ = reverse.(guides_)
        pam = reverse(pam)
    end

    pam = expand_ambiguous(pam)
    no_pam = map(x -> templates_to_sequences(x, mpt; dist = distance), guides_)
    res = zeros(Int, length(guides_), mpt.motif.distance + 1)
    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index = genomepath * ".fai") : TwoBit.Reader(ref)

    for (ic, chrom) in enumerate(gi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        #seq = ARTEMIS.getchromseq(gi.is_fa, reader[chrom])
        
        for pam_i in pam
            for (i, t) in enumerate(no_pam)
                for path in t
                    sequence = reverse(pam_i * path.seq)
                    count = ARTEMIS.count(sequence, fmi)
                    count_rve = ARTEMIS.count(reverse_complement(sequence), fmi)
                    res[i, path.dist + 1] += (count + count_rve)
                    if path.reducible != 0
                        res[i, t[path.reducible].dist + 1] -= (count + count_rve)
                    end
                end
            end
        end
    end
    #close(ref)
        
    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:distance]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end