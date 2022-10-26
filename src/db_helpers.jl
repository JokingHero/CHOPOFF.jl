const LociRange = UnitRange{UInt32}


function recalculate_ranges(x::Vector{LociRange})
    x = cumsum(length.(x))
    return LociRange.(vcat(1, x[1:end-1] .+ 1), x)
end


"
Input has to be sorted beforehand!!!
"
function ranges(guides::Vector{T}) where T <: Any
    suffixes = Vector{T}()
    loci_range = Vector{LociRange}()
    start = 1
    stop = 1
    current_guide = guides[1]
    for i in 1:length(guides)
        if guides[i] == current_guide
            stop = i
        else
            push!(loci_range, LociRange(start, stop))
            push!(suffixes, current_guide)
            current_guide = guides[i]
            start = i
            stop = i
        end
    end
    push!(loci_range, LociRange(start, stop))
    push!(suffixes, current_guide)
    return (suffixes, loci_range)
end


"
Temporary structure that is build per chromosome, per prefix.
Slow to save/load.
"
struct PrefixDB
    prefix::LongDNA{4}
    suffix::Vector{LongDNA{4}}
    loci::Vector{Loc}
end


" 
Extract 3' extension of the guide:
TTN ... EXT
CCN ... EXT
"
function getExt3(chrom::K, chrom_max::Int, ext_start::Int, dist::Int) where K <: BioSequence
    ext_end = ext_start + dist - 1
    if ext_start > chrom_max
        ext = LongDNA{4}(repeat("-", dist))
    elseif ext_end > chrom_max
        ext = chrom[ext_start:chrom_max] 
        ext = ext * LongDNA{4}(repeat("-", dist - length(ext)))
    else
        ext = chrom[ext_start:ext_end]
    end
    return ext
end

" 
Extract 5' extension of the guide:
EXT ... NAA
EXT ... NGG
"
function getExt5(chrom::K, ext_end::Int, dist::Int) where K <: BioSequence
    ext_start = ext_end - dist + 1
    if ext_end < 1
        ext = LongDNA{4}(repeat("-", dist))
    elseif ext_start < 1
        ext = chrom[1:ext_end]
        ext = LongDNA{4}(repeat("-", dist - length(ext))) * ext
    else
        ext = chrom[ext_start:ext_end]
    end
    return ext
end


function add_extension(
    guides::Vector{LongDNA{4}},
    guides_pos::Vector{UnitRange{Int64}},
    dbi::DBInfo,
    chrom::K,
    is_antisense::Bool) where K <: BioSequence

    chrom_max = lastindex(chrom)
    if dbi.motif.extends5 && is_antisense
        ext = last.(guides_pos) .+ 1
        ext = ThreadsX.map(x -> getExt3(chrom, chrom_max, x, dbi.motif.distance), ext)
        guides .= guides .* ext
    elseif dbi.motif.extends5 && !is_antisense
        ext = first.(guides_pos) .- 1
        ext = ThreadsX.map(x -> getExt5(chrom, x, dbi.motif.distance), ext)
        guides .= ext .* guides
    elseif !dbi.motif.extends5 && is_antisense
        ext = first.(guides_pos) .- 1
        ext = ThreadsX.map(x -> getExt5(chrom, x, dbi.motif.distance), ext)
        guides .= ext .* guides
    else #!dbi.motif.extends5 && !is_antisense
        ext = last.(guides_pos) .+ 1
        ext = ThreadsX.map(x -> getExt3(chrom, chrom_max, x, dbi.motif.distance), ext)
        guides .= guides .* ext
    end
    return guides
end


function normalize_to_PAMseqEXT(
    guides::Vector{LongDNA{4}},
    guides_pos::Vector{UnitRange{Int64}},
    dbi::DBInfo,
    is_antisense::Bool)

    if dbi.motif.extends5 && is_antisense
        # CCN ... EXT
        guides .= complement.(guides)
        guides_pos = first.(guides_pos)
        # becomes GGN ... EXT
    elseif dbi.motif.extends5 && !is_antisense 
        # EXT ... NGG
        guides .= reverse.(guides)
        guides_pos = last.(guides_pos)
        # becomes GGN ... EXT
    elseif !dbi.motif.extends5 && is_antisense 
        # EXT ... NAA
        guides .= reverse_complement.(guides)
        guides_pos = last.(guides_pos)
        # becomes TTA ... EXT
    else #!dbi.motif.extends5 && !is_antisense
        # TTN ... EXT
        guides_pos = first.(guides_pos)
    end
    return guides, guides_pos
end


function findguides(
    dbi::DBInfo,
    chrom::K,
    is_antisense::Bool) where {K<:BioSequence}

    query = is_antisense ? dbi.motif.rve : dbi.motif.fwd
    (seq_start, seq_stop) = locate_telomeres(chrom)
    if is_antisense ⊻ dbi.motif.extends5
        seq_start -= dbi.motif.distance
        seq_start = max(seq_start, 1)
    else
        seq_stop += dbi.motif.distance
        seq_stop = min(seq_stop, length(chrom))
    end
    return findall(query, chrom, seq_start, seq_stop; ambig_max = dbi.motif.ambig_max)
end



"
Prefix is the length of our common part between guides!
Guide reversal happens here too, so that alignment is more
convienient.
"
function gatherofftargets(
    dbi::DBInfo,
    chrom::K,
    chrom_name::String,
    is_antisense::Bool,
    prefix_len::Int) where {K<:BioSequence}
 
    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    chrom_name_ = convert(dbi.gi.chrom_type, findfirst(isequal(chrom_name), dbi.gi.chrom))

    if length(dbi.motif) != 0
        guides_pos = findguides(dbi, chrom, is_antisense)
        guides = ThreadsX.map(x -> removepam(chrom[x], pam_loci), guides_pos)
        if dbi.motif.distance > 0
            guides = add_extension(guides, guides_pos, dbi, chrom, is_antisense)
        end
        guides, guides_pos = normalize_to_PAMseqEXT(guides, guides_pos, dbi, is_antisense)
        guides_pos = convert.(dbi.gi.pos_type, guides_pos)
        loc = Loc.(chrom_name_, guides_pos, !is_antisense)

        gprefix = getindex.(guides, [1:prefix_len])
        gsuffix = getindex.(guides, [(prefix_len + 1):length(guides[1])])
        uprefixes = unique(gprefix)

        output = ThreadsX.map(uprefixes) do prefix
            prefix_guides = findall(x -> x == prefix, gprefix)
            return PrefixDB(prefix, gsuffix[prefix_guides], loc[prefix_guides])
        end
    end
    return output
end


function do_linear_chrom(chrom_name::String, chrom::K, dbi::DBInfo, prefix_len::Int, storage_dir::String) where K<:BioSequence
    @info "Working on $chrom_name"
    output_fwd = gatherofftargets(dbi, chrom, chrom_name, false, prefix_len)
    output_rve = gatherofftargets(dbi, chrom, chrom_name, true, prefix_len)
    fwd = ThreadsX.collect(x.prefix for x in output_fwd)
    rve = ThreadsX.collect(x.prefix for x in output_rve)
    prefixes = unique(vcat(fwd, rve))
    ThreadsX.map(prefixes) do prefix
        fwd_idx = findfirst(x -> x == prefix, fwd)
        rve_idx = findfirst(x -> x == prefix, rve)
        if fwd_idx === nothing
            suffixes = output_rve[rve_idx].suffix
            loci = output_rve[rve_idx].loci
        elseif rve_idx === nothing
            suffixes = output_fwd[fwd_idx].suffix
            loci = output_fwd[fwd_idx].loci
        else
            suffixes = vcat(output_fwd[fwd_idx].suffix, output_rve[rve_idx].suffix)
            loci = vcat(output_fwd[fwd_idx].loci, output_rve[rve_idx].loci)
        end
        pdb = PrefixDB(prefix, suffixes, loci)
        pdir = joinpath(storage_dir, string(prefix))
        if !isdir(pdir)
            mkdir(pdir)
        end
        save(pdb, joinpath(pdir, string(pdb.prefix) * "_" * chrom_name * ".bin"))
    end
    return prefixes
end


function getchromseq(isfa, record)
    return isfa ? FASTA.sequence(LongDNA{4}, record) : TwoBit.sequence(LongDNA{4}, record)
end