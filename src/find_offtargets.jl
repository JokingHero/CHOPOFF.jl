"
Find start and ends of the telomeres if present in the sequence, 
return start:stop of nontelomeric sequence.

e.g. NNNACTGNNN
start = 4
stop = 7
"
function locate_telomeres(seq::BioSequence, telomere_symbol = DNA_N)
    stop = length(seq)
    start = 1
    @inbounds while  start <= stop && isequal(telomere_symbol, seq[start])
        start += 1
    end

    @inbounds while stop > 0 && isequal(telomere_symbol, seq[stop])
        stop -= 1
    end

    return (start, stop)
end


"
Find all instances of pat inside seq, uses iscompatible for pattern matching.
Restrict seq to subset of start:stop and allow maximum `ambig_max` ambiguous bases.
"
function Base.findall(pat::T, seq::T,
    start::Int = 1, stop::Int = lastindex(seq); ambig_max::Int = length(pat)) where T <: BioSequence

    res = Vector{UnitRange{Int64}}()
    m = length(pat)
    n = length(seq)
    stop_ = min(stop, n) - m
    s::Int = max(start - 1, 0)

    if m == 0  # empty query
        return nothing
    end

    ambig = 0
    @inbounds while s ≤ stop_
        if iscompatible(pat[m], seq[s+m])
            i = m - 1
            if isambiguous(seq[s+m])
                ambig = 1
            else
                ambig = 0
            end
            @inbounds while i > 0
                if !iscompatible(pat[i], seq[s+i])
                    break
                end
                if isambiguous(seq[s+i])
                    ambig += 1
                end
                i -= 1
            end
            if i == 0 && ambig <= ambig_max
                push!(res, (s+1):(s+length(pat))) #found
            end
            s += 1
        else
            s += 1
        end
    end

    return res  # not found
end


function strandedguide(guide::LongDNASeq, reverse_comp::Bool, extends5::Bool)
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
add_guides!(vec::Vector{DNAMer{20}}, guides::Vector{LongDNASeq}) = append!(vec, DNAMer{20}.(guides))
add_guides!(vec::Vector{UInt128}, guides::Vector{UInt128}) = append!(vec, guides)
add_guides!(vec::Vector{UInt64}, guides::Vector{UInt64}) = append!(vec, guides)

"
Will push guides found by the dbi.motif into the output as strings.
Guides are as is for forward strand and reverse complemented when on reverse strand.
Guides do not contain PAM sequence here.
"
function pushguides!(
    output::T,
    ambig::Vector{LongDNASeq},
    dbi::DBInfo,
    chrom::K,
    reverse_comp::Bool) where {
        T<:Union{
            Vector{String}, Vector{DNAMer{20}},
            Vector{UInt128}, Vector{UInt64}}, 
        K<:BioSequence}

    query = reverse_comp ? dbi.motif.rve : dbi.motif.fwd
    pam_loci = reverse_comp ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    as_UInt128 = output isa Vector{UInt128}
    as_UInt64 = output isa Vector{UInt64}

    if length(query) != 0
        (seq_start, seq_stop) = locate_telomeres(chrom)
        if reverse_comp ⊻ dbi.motif.extends5
            seq_start -= dbi.motif.distance
            seq_start = max(seq_start, 1)
        else
            seq_stop += dbi.motif.distance
            seq_stop = min(seq_stop, lastindex(chrom))
        end
        guides = ThreadsX.map(findall(query, chrom, seq_start, seq_stop; ambig_max = dbi.motif.ambig_max)) do x 
            guide = LongDNASeq(chrom[x])
            guide = removepam(guide, pam_loci)
            guide = strandedguide(guide, reverse_comp, dbi.motif.extends5)
            return guide
        end

        if dbi.motif.ambig_max > 0
            idx = ThreadsX.map(x -> n_ambiguous(x) > 0, guides)
            push!(ambig, guides[idx])
            guides = guides[.!idx]
        end

        if as_UInt128
            guides = convert.(UInt128, guides)
        end

        if as_UInt64
            guides = convert.(UInt64, guides)
        end

        add_guides!(output, guides)
    end
    return
end


function gatherofftargets!(
    output::Any,
    dbi::DBInfo)

    ref = open(dbi.filepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = dbi.filepath * ".fai") : TwoBit.Reader(ref)
    ambig = Vector{LongDNASeq}()
    for chrom_name in dbi.chrom
        record = reader[chrom_name] # this is possible only with index!
        @info "Working on $chrom_name"
        chrom = dbi.is_fa ? FASTA.sequence(record) : TwoBit.sequence(record)
        pushguides!(output, ambig, dbi, chrom, false)
        pushguides!(output, ambig, dbi, chrom, true)
    end

    close(ref)
    return ambig
end
