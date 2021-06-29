struct LociRange
    start::UInt32
    stop::UInt32
end

import Base.length
@inline length(x::LociRange) = x.stop - x.start + 1

"
Temporary structure that is build per chromosome, per prefix.
Slow to save/load.
"
struct PrefixDB
    prefix::LongSequence{DNAAlphabet{4}}
    suffix::Vector{LongSequence{DNAAlphabet{4}}}
    loci::Vector{Loc}
end


" 
Extract 3' extension of the guide:
TTN ... EXT
CCN ... EXT
"
function getExt3(chrom::K, chrom_max::Int, ext_start::Int, dist::Int) where K <: BioSequence
    ext_end = ext_start + dist
    if ext_start > chrom_max
        ext = LongDNASeq(repeat("-", dist))
    elseif ext_end > chrom_max
        ext = chrom[ext_start:chrom_max] 
        ext = ext * LongDNASeq(repeat("-", dist - length(ext)))
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
function getExt5(chrom::K, ext_end::Int, dist::Int)  where K <: BioSequence
    ext_start = ext_end - dist + 1
    if ext_end < 1
        ext = LongDNASeq(repeat("-", dist))
    elseif ext_start < 1
        ext = chrom[1:ext_end]
        ext = LongDNASeq(repeat("-", dist - length(ext))) * ext
    else
        ext = chrom[ext_start:ext_end]
    end
    return ext
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
    reverse_comp::Bool,
    prefix_len::Int) where {K<:BioSequence}

    query = reverse_comp ? dbi.motif.rve : dbi.motif.fwd 
    pam_loci = reverse_comp ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    chrom_name_ = convert(dbi.chrom_type, findfirst(isequal(chrom_name), dbi.chrom))

    if length(query) != 0
        chrom_max = lastindex(chrom)
        guides_pos = findall(query, chrom)
        guides = Base.map(x -> removepam(chrom[x], pam_loci), guides_pos)

        # add extension
        if dbi.motif.extends5 && reverse_comp
            # CCN ... EXT
            ext = last.(guides_pos) .+ 1
            ext = Base.map(x -> getExt3(chrom, chrom_max, x, dbi.motif.distance), ext)
            guides .= complement.(guides .* ext)
            guides_pos = first.(guides_pos)
            # becomes GGN ... EXT
        elseif dbi.motif.extends5 && !reverse_comp 
            # EXT ... NGG
            ext = first.(guides_pos) .- 1
            ext = Base.map(x -> getExt5(chrom, x, dbi.motif.distance), ext)
            guides .= reverse.(ext .* guides)
            guides_pos = last.(guides_pos)
            # becomes GGN ... EXT
        elseif !dbi.motif.extends5 && reverse_comp 
            # EXT ... NAA
            ext = first.(guides_pos) .- 1
            ext = Base.map(x -> getExt5(chrom, x, dbi.motif.distance), ext)
            guides .= reverse_complement.(ext .* guides)
            guides_pos = last.(guides_pos)
            # becomes TTA ... EXT
        else #!dbi.motif.extends5 && !reverse_comp
            # TTN ... EXT
            ext = last.(guides_pos) .+ 1
            ext = Base.map(x -> getExt3(chrom, chrom_max, x, dbi.motif.distance), ext)
            guides = guides .* ext
            guides_pos = first.(guides_pos)
        end
        guides_pos = convert.(dbi.pos_type, guides_pos)
        loc = Loc.(chrom_name_, guides_pos, !reverse_comp)

        gprefix = getindex.(guides, [1:prefix_len])
        gsuffix = getindex.(guides, [prefix_len+1:length(guides[1])])
        uprefixes = unique(gprefix)

        output = ThreadsX.map(uprefixes) do prefix
            prefix_guides = findall(x -> x == prefix, gprefix)
            return PrefixDB(prefix, gsuffix[prefix_guides], loc[prefix_guides])
        end
    end
    return output
end


function do_linear_chrom(chrom_name::String, chrom::K, dbi::DBInfo, prefix_len::Int, storagedir::String) where K<:BioSequence
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
        save(pdb, joinpath(storagedir, string(pdb.prefix) * "_" * chrom_name * ".bin"))
    end
    return prefixes
end


function getchromseq(isfa, record)
    return isfa ? FASTA.sequence(record) : TwoBit.sequence(record)
end