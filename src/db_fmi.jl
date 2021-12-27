struct MotifPos
    dbi::DBInfo
    kmer_size::Int # 0 indicates there were no kmer storage
    chroms::Vector{<: Unsigned}
    pos::Vector{<: Unsigned}
    isplus::BitVector
    sequences::Vector{LongDNASeq} # large and potentially uneccessary
    ug::BitVector # bitvector encodes which guides positions share the same sequence!
    ug_count::Int
end


function Base.convert(::Type{BitVector}, x::Vector{LociRange})
    x_ = BitVector()
    bitflip = true
    for xi in x
        if bitflip
            append!(x_, BitVector(ones(length(xi))))
        else
            append!(x_, BitVector(zeros(length(xi))))
        end
        bitflip = !bitflip
    end
    return x_
end


function bits_to_counts(bits::BitVector, len::Int)
    x = zeros(Int, len)
    xi = 1
    bitflip = true
    i = 0
    for bi in bits
        if bi == bitflip
            i += 1
        else
            x[xi] = i
            bitflip = !bitflip
            xi += 1
            i = 1
        end
    end
    x[xi] = i
    return x
end


function get_ug_ranges(bv::BitVector, idx::Int)
    val = bv[idx]
    len = length(bv)
    start = idx
    while start != 1 && val == bv[start - 1]
        start -= 1
    end
    stop = idx
    while stop != len && val == bv[stop + 1]
        stop += 1
    end
    return (start, stop)
end


function gatherofftargets(
    dbi::DBInfo,
    chrom::K,
    chrom_name::String,
    reverse_comp::Bool) where {K<:BioSequence}
 
    pam_loci = reverse_comp ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    chrom_name_ = convert(dbi.chrom_type, findfirst(isequal(chrom_name), dbi.chrom))

    if length(dbi.motif) != 0
        guides_pos = findguides(dbi, chrom, reverse_comp)
        guides = ThreadsX.map(x -> removepam(chrom[x], pam_loci), guides_pos)
        guides = add_extension(guides, guides_pos, dbi, chrom, reverse_comp)
        guides, guides_pos = normalize_to_PAMseqEXT(guides, guides_pos, dbi, reverse_comp)
        guides_pos = convert.(dbi.pos_type, guides_pos)
        return (chrom_name_, guides, guides_pos, !reverse_comp)
    end
    return nothing
end


function build_motifDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String;
    store_kmers = true)

    dbi = CRISPRofftargetHunter.DBInfo(genomepath, name, motif)
    prefix_len = 4

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
    guides_all = Vector{LongDNASeq}()
    gbits_all = Vector{BitVector}()
    loci_range_all = Vector{LociRange}()
    loci_chrom_all = Vector{dbi.chrom_type}()
    loci_pos_all = Vector{dbi.pos_type}()
    loci_isplus_all = BitVector()
    kmer_size = minkmersize(length_noPAM(motif), motif.distance)
    kmers = all_kmers(kmer_size)
    kmers = IdDict(zip(kmers, 1:length(kmers)))
    @showprogress 60 for prefix in prefixes
        guides = Vector{LongDNASeq}()
        loci = Vector{Loc}()
        for chrom in dbi.chrom
            p = joinpath(storagedir, string(prefix), string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                append!(guides, pdb.suffix)
                append!(loci, pdb.loci)
            end
        end
        rm(joinpath(storagedir, string(prefix)), recursive = true)
        (guides, loci_range, loci) = unique_guides(guides, loci)
        guides = Base.map(x -> prefix * x, guides)
        if store_kmers
            guides_ = ThreadsX.map(x -> as_bitvector_of_kmers(x, kmers), guides)
            append!(gbits_all, guides_)
        else
            kmer_size = 0  
        end
        append!(guides_all, guides)
        append!(loci_range_all, loci_range)
        append!(loci_chrom_all, map(x -> x.chrom, loci))
        append!(loci_pos_all, map(x -> x.pos, loci))
        append!(loci_isplus_all, BitVector(map(x -> x.isplus, loci)))
    end
    ug_count = length(loci_range_all)
    loci_range_all = convert(BitVector, loci_range_all)
    motifpos = MotifPos(
        dbi, kmer_size, loci_chrom_all, loci_pos_all, 
        loci_isplus_all, guides_all, loci_range_all, ug_count)
    motifpos_path = joinpath(storagedir, "motifDB.bin")
    save(motifpos, motifpos_path)

    if store_kmers
        # now kmers for each guide - skip FMindex stuff - pay with disk space
        @info "Writing down kmers."
        @showprogress 60 for (kmer, pos) in kmers
            kmer_bit = zeros(length(gbits_all))
            for (idx, g) in enumerate(gbits_all)
                kmer_bit[idx] = g[pos]
            end
            save(BitVector(kmer_bit), joinpath(storagedir, string(kmer) * ".bin"))
        end
    end

    @info "Finished constructing motifDB in " * storagedir
    return storagedir
end


function build_fmiDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String)

    dbi = DBInfo(genomepath, name, motif)
    ref = open(dbi.filepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = dbi.filepath * ".fai") : TwoBit.Reader(ref)
    @showprogress 60 for chrom in dbi.chrom
        fmi = FMIndex(getchromseq(dbi.is_fa, reader[chrom]), 16, r = 32)
        p = joinpath(storagedir, chrom * ".bin")
        save(fmi, p)
    end
    close(ref)
    save(dbi, joinpath(storagedir, "dbi.bin"))
    @info "FMindex is build."
    return storagedir
end


function guide_to_bitvector(guide::LongDNASeq, path::String, kmer_size::Int)
    guide = collect(Set(as_kmers(guide, kmer_size)))
    guide_bits = load(joinpath(path, string(guide[1]) * ".bin"))
    for i in 2:length(guide)
        bits = load(joinpath(path, string(guide[i]) * ".bin"))
        guide_bits .|= bits
    end
    return guide_bits
end


function search_motifDB(
    storagedir::String, 
    guides::Vector{LongDNASeq}, 
    dist::Int = 4; 
    detail::String = "")

    mp = load(joinpath(storagedir, "motifDB.bin"))
    mp_counts = bits_to_counts(mp.ug, mp.ug_count)

    if dist > mp.dbi.motif.distance
        error("Maximum distance is " * string(mp.dbi.motif.distance))
    end

    # we work on each chrom separately
    # each guide separately?

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if mp.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = zeros(Int, length(guides_), dist + 1)
    for (idx, g) in enumerate(guides_) # if we add detail and early stopping - paralelize on the guides
        g_bits = guide_to_bitvector(g, storagedir, mp.kmer_size)
        offtargets = mp.sequences[g_bits]
        mp_counts_g = mp_counts[g_bits]
        dists = ThreadsX.map(x -> levenshtein(g, x, dist), offtargets)
        for d in 0:dist
            res[idx, d + 1] += sum(mp_counts_g[dists .== d])
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
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
    guides::Vector{LongDNASeq}, dist::Int = 4; detail::String = "")

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
    for (ic, chrom) in enumerate(mp.dbi.chrom)
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


function as_partial_alignments(s::String, motif::Motif, len::Int = 10)
    s = s[1:(len + d)]
    if d == 0
        return [s]
    end
    comb = comb_of_d1(s, alphabet)
    for i in 1:(d-1)
        comb = foldxt(union, Map(x -> comb_of_d1(x, alphabet)), comb)
    end

    comb = collect(Set(map(x -> x[1:len], comb)))

    if motif.extends5
        partials = map(x -> dna"GGN" * x, comb)
        partials_rev = complement(partials)
        partials = reverse(partials)
    else
        partials = map(x -> dna"TTTN" * x, comb)
        partials_rev = reverse_complement(partials)
    end
    partials = expand_ambiguous(partials)
    partials_rev = expand_ambiguous(partials_rev)
    return (partials, partials_rev)
end

function search_fmiDB_raw(
    fmidbdir::String, genomepath::String, motif::Motif,
    guides::Vector{LongDNASeq}; detail::String = "", prefix_length::Int = 10)

    dbi = load(joinpath(storagedir, "dbi.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if motif.extends5
        guides_ = reverse.(guides_)
    end

    partials = as_partial_alignments.(String.(guides_), motif, prefix_length)
    
    res = zeros(Int, length(guides_), motif.distance + 1)
    ref = open(genomepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = genomepath * ".fai") : TwoBit.Reader(ref)

    for (ic, chrom) in enumerate(dbi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        seq = getchromseq(dbi.is_fa, reader[chrom])
    
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