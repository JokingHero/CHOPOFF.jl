"
Find all instances of pat inside seq, disregard seq ambiguous bases.
Restrict seq to subset of start:stop.
"
# TODO we might want to not disregard ambigous bases when they are on
# the extension! - as this could potentially lose couple overlapping
# endings NNNNN off-targets
function Base.findall(pat::BioSequence, seq::BioSequence,
    start::Integer = 1, stop::Integer = lastindex(seq))

    res = Vector{UnitRange{Int64}}()
    m = length(pat)
    n = length(seq)
    stop_ = min(stop, n) - m
    s::Int = max(start - 1, 0)

    if m == 0  # empty query
        return nothing
    end

    while s ≤ stop_
        if isambiguous(seq[s+m])
            s += m
        elseif iscompatible(pat[m], seq[s+m])
            i = m - 1
            while i > 0
                if isambiguous(seq[s+i]) || !iscompatible(pat[i], seq[s+i])
                    break
                end
                i -= 1
            end
            if i == 0
                push!(res, (s+1):(s+length(pat))) #found
            end
            s += 1
        else
            s += 1
        end
    end

    return res  # not found
end


function pushguides!(
    output::T,
    dbi::DBInfo,
    chrom::K,
    reverse_comp::Bool) where {T<:Union{CountMinSketch, HyperLogLog, Vector{String}}, K<:BioSequence}

    query = reverse_comp ? dbi.motif.fwd : dbi.motif.rve
    pam_loci = reverse_comp ? dbi.motif.pam_loci_fwd : dbi.motif.pam_loci_rve

    if length(query) != 0
        for x in findall(query, chrom)
            guide = LongDNASeq(chrom[x])
            guide = removepam(guide, pam_loci)
            if reverse_comp
                guide = reverse_complement(guide)
            end

            push!(output, string(guide))
        end
    end
    return output
end


function pushguides!(
    output::IdDict{String, Int},
    dbi::DBInfo,
    chrom::K,
    reverse_comp::Bool) where {K<:BioSequence}

    query = reverse_comp ? dbi.motif.fwd : dbi.motif.rve
    pam_loci = reverse_comp ? dbi.motif.pam_loci_fwd : dbi.motif.pam_loci_rve

    g_len = length(query) - length(collect(vcat(pam_loci...)))
    kmer_size = minkmersize(g_len, dbi.motif.distance)

    if length(query) != 0
        for x in findall(query, chrom)
            guide = LongDNASeq(chrom[x])
            guide = removepam(guide, pam_loci)
            if reverse_comp
                guide = reverse_complement(guide)
            end

            # kmer treatment - just unique count
            kmers = Vector{String}()
            for i in 1:g_len-kmer_size
                kmer = string(guide[i:(i + kmer_size - 1)])
                if !(kmer in kmers)
                    push!(kmers, kmer)
                    output[kmer] = get(output, kmer, 0) + 1
                end
            end
        end
    end
    return output
end


function gatherofftargets!(
    output::Any,
    dbi::DBInfo;
    prefix_len::Union{Int, Nothing} = nothing)

    ref = open(dbi.filepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = dbi.filepath * ".fai") : TwoBit.Reader(ref)

    for chrom_name in dbi.chrom
        record = reader[chrom_name] # this is possible only with index!
        @info "Working on $chrom_name"
        chrom = dbi.is_fa ? FASTA.sequence(record) : TwoBit.sequence(record)

        if isnothing(prefix_len)
            pushguides!(output, dbi, chrom, false)
            pushguides!(output, dbi, chrom, true)
        else
            pushguides!(output, dbi, chrom, chrom_name, false, prefix_len)
            pushguides!(output, dbi, chrom, chrom_name, true, prefix_len)
        end
    end

    close(ref)
    return
end


# SHENANIGANS

function countofftargets(
    chrom::K,
    query::LongDNASeq,
    pam_loci::Vector{UnitRange{<:Integer}},
    max_dist::Int,
    reverse_comp::Bool,
    guides::Vector{LongDNASeq}) where {K<:BioSequence}

    output = zeros(Int, length(guides), max_dist + 1)

    if length(query) != 0
        for x in findall(query, chrom)
            guide_ref = LongDNASeq(chrom[x])
            guide_ref = removepam(guide_ref, pam_loci)
            if reverse_comp
                guide_ref = reverse_complement(guide_ref)
            end

            for (i, gi) in enumerate(guides)
                dist = levenshtein(reverse(gi), reverse(guide_ref), max_dist)
                if dist <= max_dist
                    output[i, dist + 1] += 1
                end
            end
        end
    end
    return output
end


function do_chrom(record, motif, max_dist, g)
    res = zeros(Int, length(g), max_dist + 1)
    chrom = FASTA.sequence(record)

    if FASTA.hasidentifier(record) # TODO 2bit?!
        @info "Working on: " * FASTA.identifier(record)
    end

    res += countofftargets(chrom, motif.fwd, motif.pam_loci_fwd, max_dist, false, g)
    res += countofftargets(chrom, motif.rve, motif.pam_loci_rve, max_dist, true, g)

    return res
end


function getdist(guide::String, ref::String)
    if levenshtein_bp(guide, ref, 4)
        return levenshtein(LongDNASeq(guide), LongDNASeq(ref), 4)
    else
        return 5
    end
end

# I checked whether kmers strategy for prefiltering of offtargets
# is faster but around 0.7 of time is wasted and gains
# are small, and even negative especially when set construction is
# needed
function iterate_over_offtargets(guides, all_guides, oft_file)
    guides_str = [string(g) for g in guides]

    res = zeros(Int, length(guides), 5)
    f = open(all_guides)
    out = open(oft_file, "w")
    i = 0
    w = 0
    for ln in eachline(f)
        if ln == "guide,location"
            continue
        end
        i += 1
        w += 1
        if w == 1000000
            println("Iter: ", i)
            w = 0
        end
        ln = split(ln, ",")
        ref = String(ln[1][4:end])
        loci = ln[2]

        dist = ThreadsX.map(g -> getdist(g, ref), guides_str)
        for j in 1:length(dist)
            if dist[j] <= 4
                res[j, dist[j] + 1] += 1
                write(out, guides_str[j] * "," * ref * "," * string(dist[j]) * "," * loci * "\n")
            end
        end
    end
    close(f)
    close(out)
    return res
end
