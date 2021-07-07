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


function strandedguide(guide, reverse_comp::Bool, extends5::Bool)
    if extends5 && reverse_comp
        # CCN ... EXT
        guide = complement(guide)
        # becomes GGN ... EXT
    elseif extends5 && !reverse_comp 
        # EXT ... NGG
        guide = reverse(guide)
        # becomes GGN ... EXT
    elseif !extends5 && reverse_comp 
        # EXT ... NAA
        guide = reverse_complement(guide)
        # becomes TTA ... EXT
    end
    return guide
end


add_guides!(vec::Vector{String}, guides::Vector{LongDNASeq}) = append!(vec, string.(guides))
add_guides!(vec::Vector{LongDNASeq}, guides::Vector{LongDNASeq}) = append!(vec, guides)
add_guides!(vec::Vector{UInt64}, guides::Vector{LongDNASeq}) = append!(vec, unsigned.(DNAMer.(guides)))

"
Will push guides found by the dbi.motif into the output as strings.
Guides are as is for forward strand and reverse complemented when on reverse strand.
Guides do not contain PAM sequence here.
"
function pushguides!(
    output::T,
    dbi::DBInfo,
    chrom::K,
    reverse_comp::Bool) where {
        T<:Union{
            IdDict, CountMinSketch, HyperLogLog, 
            Vector{String}, Vector{LongDNASeq}, Vector{UInt64}}, 
        K<:BioSequence}
    
    query = reverse_comp ? dbi.motif.rve : dbi.motif.fwd
    pam_loci = reverse_comp ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd

    if length(query) != 0
        guides = ThreadsX.map(findall(query, chrom)) do x 
            guide = LongDNASeq(chrom[x])
            guide = removepam(guide, pam_loci)
            # we don't need extensions here
            guide = strandedguide(guide, reverse_comp, dbi.motif.extends5)
            guide
        end
        add_guides!(output, guides)
    end
    return output
end


function gatherofftargets!(
    output::Any,
    dbi::DBInfo)

    ref = open(dbi.filepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = dbi.filepath * ".fai") : TwoBit.Reader(ref)

    for chrom_name in dbi.chrom
        record = reader[chrom_name] # this is possible only with index!
        @info "Working on $chrom_name"
        chrom = dbi.is_fa ? FASTA.sequence(record) : TwoBit.sequence(record)
        pushguides!(output, dbi, chrom, false)
        pushguides!(output, dbi, chrom, true)
    end

    close(ref)
    return
end
