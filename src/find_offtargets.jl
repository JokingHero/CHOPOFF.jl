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


function strandedguide(guide::LongDNA{4}, is_antisense::Bool, extends5::Bool)
    if extends5 && is_antisense
        # CCN ... EXT
        guide = complement(guide)
        # becomes GGN ... EXT
    elseif extends5 && !is_antisense 
        # EXT ... NGG
        guide = reverse(guide)
        # becomes GGN ... EXT
    elseif !extends5 && is_antisense 
        # EXT ... NAA
        guide = reverse_complement(guide)
        # becomes TTA ... EXT
    end
    return guide
end


add_guides!(vec::Vector{String}, guides::Vector{LongDNA{4}}) = append!(vec, string.(guides))
add_guides!(vec::Vector{UInt128}, guides::Vector{UInt128}) = append!(vec, guides)
add_guides!(vec::Vector{UInt64}, guides::Vector{UInt64}) = append!(vec, guides)


"
Will push guides found by the dbi.motif into the output as strings.
Guides are as is for forward strand and reverse complemented when on reverse strand.
Guides do not contain PAM sequence here.

remove_pam - whether PAM sequence should be removed
normalzie - whether all guides should be flipped into PAMseqEXT e.g. GGn-20N-3bp
"
function pushguides!(
    output::T,
    ambig::Vector{LongDNA{4}},
    dbi::DBInfo,
    chrom::K,
    is_antisense::Bool;
    remove_pam::Bool = true,
    normalize::Bool = true) where {
        T<:Union{
            Vector{String},
            Vector{UInt128}, 
            Vector{UInt64}}, 
        K<:BioSequence}

    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    as_UInt128 = output isa Vector{UInt128}
    as_UInt64 = output isa Vector{UInt64}

    if length(dbi.motif) != 0
        guides_pos = findguides(dbi, chrom, is_antisense)
        guides = ThreadsX.map(guides_pos) do x 
            guide = LongDNA{4}(chrom[x])
            if remove_pam 
                guide = removepam(guide, pam_loci)
            end
            return guide
        end

        if dbi.motif.distance > 0
            guides = add_extension(guides, guides_pos, dbi, chrom, is_antisense)
        end

        if normalize
            guides, guides_pos = normalize_to_PAMseqEXT(guides, guides_pos, dbi, is_antisense)
        end
        guides_pos = nothing # not needed for anything and uses a lot of space
        
        idx = ThreadsX.map(isambig, guides)
        if sum(idx) != 0
            append!(ambig, guides[idx])
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


"""
```
gatherofftargets!(
    output::T,
    dbi::DBInfo) where {T<:Union{Vector{String}, Vector{UInt64}, Vector{UInt128}}}
```

Gathers all off-targets that conform to the given `dbi` Motif.

This function appends to the `output` during the run, however it will also return all ambiguous 
guides in return object. We can use UInt64 and UInt128 to compress space that the gRNAs use. When using 
large genomes or non-specific PAMs you might run out of memory when using this function.

# Examples
```julia-repl
# use ARTEMIS example genome
genome = joinpath(
    vcat(
        splitpath(dirname(pathof(ARTEMIS)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))
# construct example DBInfo
dbi = DBInfo(genome, "Cas9_semirandom_noVCF", Motif("Cas9"))
# finally gather all off-targets
guides = Vector{String}()
ambig = gatherofftargets!(guides, dbi)

# here in the format of UInt64 encoding
guides2 = Vector{UInt64}()
ambig2 = gatherofftargets!(guides2, dbi)
guide_with_extension_len = length_noPAM(dbi.motif) + dbi.motif.distance

# transform UInt64 to LongDNA and String
guides2 = String.(LongDNA{4}.(guides2, guide_with_extension_len))
```
"""
function gatherofftargets!(
    output::T,
    dbi::DBInfo) where {T<:Union{Vector{String}, Vector{UInt64}, Vector{UInt128}}}

    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
    ambig = Vector{LongDNA{4}}()
    for chrom_name in dbi.gi.chrom
        record = reader[chrom_name] # this is possible only with index!
        @info "Working on $chrom_name"
        chrom = dbi.gi.is_fa ? FASTA.sequence(LongDNA{4}, record) : TwoBit.sequence(LongDNA{4}, record)
        pushguides!(output, ambig, dbi, chrom, false)
        pushguides!(output, ambig, dbi, chrom, true)
    end

    close(ref)
    return ambig
end


function gatherofftargets(
    seq::LongDNA{4},
    dbi::DBInfo)

    guides_pos_fwd = findguides(dbi, seq, false)
    guides_fwd = ThreadsX.map(guides_pos_fwd) do x 
        guide = LongDNA{4}(seq[x])
        guide = removepam(guide, dbi.motif.pam_loci_fwd)
        return guide
    end
    guides_fwd = add_extension(guides_fwd, guides_pos_fwd, dbi, seq, false)
    guides_fwd, guides_pos = normalize_to_PAMseqEXT(guides_fwd, guides_pos_fwd, dbi, false)
    

    guides_pos_rve = findguides(dbi, seq, true)
    guides_rve = ThreadsX.map(guides_pos_rve) do x 
        guide = LongDNA{4}(seq[x])
        guide = removepam(guide, dbi.motif.pam_loci_rve)
        return guide
    end
    guides_rve = add_extension(guides_rve, guides_pos_rve, dbi, seq, true)
    guides_rve, guides_pos = normalize_to_PAMseqEXT(guides_rve, guides_pos_rve, dbi, true)

    return vcat(guides_fwd, guides_rve), vcat(guides_pos_fwd, guides_pos_rve)
end
