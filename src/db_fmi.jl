struct MotifPos
    dbi::DBInfo
    kmer_size::Int
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
            p = joinpath(storagedir, string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                append!(guides, pdb.suffix)
                append!(loci, pdb.loci)
                rm(p)
            end
        end
        (guides, loci_range, loci) = unique_guides(guides, loci)
        if store_kmers
            guides = Base.map(x -> prefix * x, guides)
            guides_ = ThreadsX.map(x -> as_bitvector_of_kmers(x, kmers), guides)
            append!(gbits_all, guides_)
        else
            guides = nothing  
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

    @info "Finished constructing motifDB in " * motifpos_path
    return motifpos_path
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
        fmi = FMIndex(getchromseq(dbi.is_fa, reader[x]), 16, r = 32)
        p = joinpath(storagedir, chrom * ".bin")
        save(fmi, p)
    end
    close(ref)
    @info "FMindex is build."
    return storagedir
end


function guide_to_bitvector(guide::LongDNASeq, path::String, kmer_size::Int)
    guide = collect(as_kmers(guide, kmer_size))
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


function search_fmiDB(
    fmidir::String, motifdbpath::String, 
    guides::Vector{LongDNASeq}, dist::Int = 4; detail::String = "")

    mp = load(motifdbpath)

    if dist > mp.dbi.motif.distance
        error("Maximum distance is " * string(mp.dbi.motif.distance))
    end

    # we work on each chrom separately
    # each guide separately?
    
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

    #res = zeros(Int, length(guides_), dist + 1)
    res = ThreadsX.mapreduce(p -> search_prefix(p, dist, ldb.dbi, dirname(detail), guides_, storagedir), +, prefixes)
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