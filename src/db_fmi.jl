struct SeedDB
    dbi::DBInfo
    prefixes::Set{LongDNASeq}
    kmer_size::Int
    kmers::Dict{LongDNASeq, Int}
end


struct MotifPos
    chroms::Vector{<: Unsigned}
    pos::Vector{<: Unsigned}
    isplus::BitVector
    sequences::Vector{UInt128} # large and potentially uneccessary suffixes
    ug::BitVector # bitvector encodes which guides positions share the same sequence!
    ug_count::Int
    bits::BitMatrix # columns are guides, rows are kmers
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
    chrom_name_ = convert(dbi.gi.chrom_type, findfirst(isequal(chrom_name), dbi.gi.chrom))

    if length(dbi.motif) != 0
        guides_pos = findguides(dbi, chrom, reverse_comp)
        guides = ThreadsX.map(x -> removepam(chrom[x], pam_loci), guides_pos)
        guides = add_extension(guides, guides_pos, dbi, chrom, reverse_comp)
        guides, guides_pos = normalize_to_PAMseqEXT(guides, guides_pos, dbi, reverse_comp)
        guides_pos = convert.(dbi.gi.pos_type, guides_pos)
        return (chrom_name_, guides, guides_pos, !reverse_comp)
    end
    return nothing
end


function build_motifDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String,
    prefix_len = 7)

    dbi = CRISPRofftargetHunter.DBInfo(genomepath, name, motif)

    # step 1
    @info "Step 1: Searching chromosomes."
    # For each chromsome paralelized we build database
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
    # Don't paralelize here as you can likely run out of memory (chromosomes are large)
    prefixes = Base.map(x -> do_linear_chrom(x, getchromseq(dbi.gi.is_fa, reader[x]), dbi, prefix_len, storagedir), dbi.gi.chrom)
    close(ref)

    prefixes = Set(vcat(prefixes...))

    # step 2
    @info "Step 2: Constructing per prefix db."
    # Iterate over all prefixes and merge different chromosomes
    kmer_size = Int(floor(length_noPAM(motif) / (motif.distance + 2))) # lossless seed 01*0
    kmers = LongDNASeq.(all_kmers(kmer_size))
    kmers = Dict(zip(kmers, 1:length(kmers)))
    @showprogress 60 for prefix in prefixes # can be paralelized here ?! memory?!
        guides = Vector{LongDNASeq}()
        loci = Vector{Loc}()
        for chrom in dbi.gi.chrom
            p = joinpath(storagedir, string(prefix), string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                append!(guides, pdb.suffix)
                append!(loci, pdb.loci)
            end
        end
        rm(joinpath(storagedir, string(prefix)), recursive = true)
        (guides, loci_range, loci) = unique_guides(guides, loci)
        guides_ = ThreadsX.map(x -> as_bitvector_of_kmers(prefix * x, kmers), guides)
        guides_ = hcat(guides_...) # as BitMatrix - columns are guides, rows are kmers

        loci_chrom = ThreadsX.map(x -> x.chrom, loci)
        loci_pos = ThreadsX.map(x -> x.pos, loci)
        loci_isplus = BitVector(ThreadsX.map(x -> x.isplus, loci))

        ug_count = length(loci_range)
        loci_range = convert(BitVector, loci_range)
        guides = ThreadsX.map(x -> convert(UInt128, x), guides)
        save(
            MotifPos(
                loci_chrom, loci_pos, loci_isplus, guides,
                loci_range, ug_count, guides_),
            joinpath(storagedir, string(prefix) * ".bin"))
    end

    save(SeedDB(dbi, prefixes, kmer_size, kmers), joinpath(storagedir, "seedDB.bin"))
    @info "Finished constructing seedDB in " * storagedir
    return storagedir
end


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





function guide_to_bitvector(guide::Vector{LongDNASeq}, bits::BitMatrix, kmers::Dict{LongDNASeq, Int})
    in_guide = bits[kmers[guide[1]], :]
    for i in 2:length(guide)
        in_guide += bits[kmers[guide[i]], :]
    end
    return in_guide .> 1 # lossless seed requiries at least two skipmers to be found
end


function search_prefix(
    prefix::LongDNASeq,
    dist::Int,
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNASeq},
    gskipmers::Vector{Vector{LongDNASeq}},
    kmers::Dict{LongDNASeq, Int},
    storagedir::String)

    res = zeros(Int, length(guides), dist + 1)
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
    sdb = load(joinpath(storagedir, string(prefix) * ".bin"))
    sdb_counts = bits_to_counts(sdb.ug, sdb.ug_count)
    for i in 1:length(guides)
        if !isfinal[i]
            g_bits = guide_to_bitvector(gskipmers[i], sdb.bits, kmers)
            if any(g_bits)
                offtargets = map(x -> LongDNASeq(x, suffix_len), sdb.sequences[g_bits])
                sdb_counts_g = sdb_counts[g_bits]
                for (j, suffix) in enumerate(offtargets)
                    suffix_aln = suffix_align(suffix, prefix_aln[i])
                    if suffix_aln.dist <= dist
                        res[i, suffix_aln.dist + 1] += sdb_counts_g[j]
                        if detail != ""
                            #=
                            offtargets = sdb.loci[sl_idx.start:sl_idx.stop]
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
                                write(detail_file, noloc * decode(offt, dbi) * "\n")
                            end
                            =#
                        end
                    end
                end
            end
        end
    end

    if detail != ""
        close(detail_file)
    end
    return res
end


function search_motifDB(
    storagedir::String, 
    guides::Vector{LongDNASeq}, 
    dist::Int = 4; 
    detail::String = "")

    sdb = load(joinpath(storagedir, "seedDB.bin"))
    if dist > sdb.dbi.motif.distance
        error("Maximum distance is " * string(sdb.dbi.motif.distance))
    end

    # we work on each chrom separately
    # each guide separately?

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    gskipmers = ThreadsX.map(x -> collect(Set(as_skipkmers(x, sdb.kmer_size))), guides_)
    #res = zeros(Int, length(guides_), dist + 1)
    res = ThreadsX.mapreduce(p -> search_prefix(
        p, dist, sdb.dbi, dirname(detail), guides_, gskipmers, sdb.kmers, storagedir), +, sdb.prefixes)
    #for p in prefixes
    #    res += search_prefix(p, dist, ldb.dbi, dirname(detail), guides_, storagedir)
    #end
    
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
    guides::Vector{LongDNASeq}; detail::String = "", prefix_length::Int = 10)

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
    guides::Vector{LongDNASeq}; detail::String = "", distance::Int = 3)

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
        #seq = CRISPRofftargetHunter.getchromseq(gi.is_fa, reader[chrom])
        
        for pam_i in pam
            for (i, t) in enumerate(no_pam)
                for path in t
                    sequence = reverse(pam_i * path.seq)
                    count = CRISPRofftargetHunter.count(sequence, fmi)
                    count_rve = CRISPRofftargetHunter.count(reverse_complement(sequence), fmi)
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


function search_fmiDB_patterns_cashed(
    fmidbdir::String, genomepath::String, mpt::MotifPathTemplates,
    guides::Vector{LongDNASeq}; detail::String = "", distance::Int = 3)

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
    # we need to add GGN in front and reverse for it to be forward sequence to search on the genome
    # for reverse we need to complement

    res = zeros(Int, length(guides_), mpt.motif.distance + 1)
    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index = genomepath * ".fai") : TwoBit.Reader(ref)

    for (ic, chrom) in enumerate(gi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        cashe = Dict{LongDNASeq, UnitRange{Int}}()
        #seq = CRISPRofftargetHunter.getchromseq(gi.is_fa, reader[chrom])
        
        for pam_i in pam
            for (i, t) in enumerate(no_pam)
                for path in t
                    sequence = reverse(pam_i * path.seq)
                    count = count_with_cashe!(sequence, fmi, cashe)
                    count_rve = count_with_cashe!(reverse_complement(sequence), fmi, cashe)
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