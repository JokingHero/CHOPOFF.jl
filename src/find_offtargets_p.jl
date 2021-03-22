
function prepare_for_align(
    guide_ref::LongDNASeq,
    pam_loci::Vector{UnitRange{<:Integer}},
    reverse_comp::Bool)

    guide_ref = removepam(guide_ref, pam_loci)
    if reverse_comp
        guide_ref = reverse_complement(guide_ref)
    end
    guide_ref
end


function align_to_ref(
    guide_ref::LongDNASeq,
    guides::Vector{LongDNASeq},
    max_dist::Int)
    output = zeros(Int, length(guides), max_dist + 1)
    @inbounds for (i, gi) in enumerate(guides)
        dist = levenshtein(reverse(gi), reverse(guide_ref), max_dist)
        if dist <= max_dist
            output[i, dist + 1] += 1
        end
    end
    return output
end


function countofftargets_p(
    chrom::K,
    query::LongDNASeq,
    pam_loci::Vector{UnitRange{<:Integer}},
    max_dist::Int,
    reverse_comp::Bool,
    guides::Vector{LongDNASeq}) where {K<:BioSequence}

    output = zeros(Int, length(guides), max_dist + 1)

    if length(query) != 0
        guide_ref = findall(query, chrom)
        guide_ref = [prepare_for_align(chrom[x], pam_loci, reverse_comp) for x in guide_ref]
        output = ThreadsX.mapreduce(x -> align_to_ref(x, guides, max_dist), +, guide_ref)
    end

    return output
end


function findofftargets_p_chrom(
    genome::String,
    motif::Motif, # has to be motive including the max_dist extensions on the 3'
    max_dist::Int,
    g::Array{LongDNASeq, 1}) # vector of guides in 1-20 (NGG) config)

    ext = extension(genome)
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        error(
        "Wrong extension of the genome.",
        "We can parse only .fa/.fna/.fasta or .2bit references.")
    end

    ref = open(genome, "r")
    reader = is_fa ? FASTA.Reader(ref) : TwoBit.Reader(ref)
    res_vec = ThreadsX.mapi(x -> do_chrom(x, motif, max_dist, g), reader)
    close(ref)

    res = reduce(+, res_vec)
    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:max_dist]
    rename!(res, col_d)
    res.guide = [string(gi) for gi in g]
    sort!(res, col_d)
    return res
end


function do_chrom_p(record, motif, max_dist, g)
    res = zeros(Int, length(g), max_dist + 1)
    chrom = FASTA.sequence(record)

    if FASTA.hasidentifier(record) # TODO 2bit?!
        println("Working on: ", FASTA.identifier(record))
    end

    res += countofftargets_p(chrom, motif.fwd, motif.pam_loci_fwd, max_dist, false, g)
    res += countofftargets_p(chrom, motif.rve, motif.pam_loci_rve, max_dist, true, g)

    return res
end


function findofftargets_p_refg(
    genome::String,
    motif::Motif, # has to be motive including the max_dist extensions on the 3'
    max_dist::Int,
    g::Array{LongDNASeq, 1}) # vector of guides in 1-20 (NGG) config)

    ext = extension(genome)
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        error(
        "Wrong extension of the genome.",
        "We can parse only .fa/.fna/.fasta or .2bit references.")
    end

    ref = open(genome, "r")
    reader = is_fa ? FASTA.Reader(ref) : TwoBit.Reader(ref)
    res_vec = Base.map(x -> do_chrom_p(x, motif, max_dist, g), reader)
    close(ref)

    res = reduce(+, res_vec)
    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:max_dist]
    rename!(res, col_d)
    res.guide = [string(gi) for gi in g]
    sort!(res, col_d)
    return res
end
