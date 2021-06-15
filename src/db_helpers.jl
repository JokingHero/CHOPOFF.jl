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
    suffix::Dict{LongSequence{DNAAlphabet{4}}, Vector{Loc}}
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
function pushguides!(
    output::Vector{PrefixDB},
    dbi::DBInfo,
    chrom::K,
    chrom_name::String,
    reverse_comp::Bool,
    prefix_len::Int) where {K<:BioSequence}

    query = reverse_comp ? dbi.motif.rve : dbi.motif.fwd 
    pam_loci = reverse_comp ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd

    # This can't be paralelized with Threads and lock
    if length(query) != 0
        chrom_max = lastindex(chrom)
        for x in findall(query, chrom)
            guide = LongDNASeq(chrom[x])
            guide = removepam(guide, pam_loci)
            # add extension
            if dbi.motif.extends5 && reverse_comp
                # CCN ... EXT
                guide = guide * getExt3(chrom, chrom_max, last(x) + 1, dbi.motif.distance)
                guide = complement(guide)
                pos_ = first(x)
                # becomes GGN ... EXT
            elseif dbi.motif.extends5 && !reverse_comp 
                # EXT ... NGG
                guide = getExt5(chrom, first(x) - 1, dbi.motif.distance) * guide
                guide = reverse(guide)
                pos_ = last(x)
                # becomes GGN ... EXT
            elseif !dbi.motif.extends5 && reverse_comp 
                # EXT ... NAA
                guide = getExt5(chrom, first(x) - 1, dbi.motif.distance) * guide
                guide = reverse_complement(guide)
                pos_ = last(x)
                # becomes TTA ... EXT
            else #!dbi.motif.extends5 && !reverse_comp
                # TTN ... EXT
                guide = guide * getExt3(chrom, chrom_max, last(x) + 1, dbi.motif.distance)
                pos_ = first(x)
            end

            gprefix = guide[1:prefix_len]
            gsuffix = guide[prefix_len+1:end]
            chrom_name_ = convert(dbi.chrom_type, findfirst(isequal(chrom_name), dbi.chrom))
            pos_ = convert(dbi.pos_type, pos_)
            loc = Loc(chrom_name_, pos_, !reverse_comp)
            prefix_idx = findfirst(x -> gprefix == x.prefix, output)

            if prefix_idx !== nothing
                if haskey(output[prefix_idx].suffix, gsuffix)
                    push!(output[prefix_idx].suffix[gsuffix], loc)
                else
                    output[prefix_idx].suffix[gsuffix] = [loc]
                end
            else
                push!(output, PrefixDB(gprefix, Dict(gsuffix => [loc])))
            end
        end
    end
    return output
end


function do_linear_chrom(chrom_name::String, chrom::K, dbi::DBInfo, prefix_len::Int, storagedir::String) where K<:BioSequence
    @info "Working on $chrom_name"
    output = Vector{PrefixDB}()
    pushguides!(output, dbi, chrom, chrom_name, false, prefix_len)
    pushguides!(output, dbi, chrom, chrom_name, true, prefix_len)
    # save small files as they are too large to store them in memory
    for pdb in output
        save(pdb, joinpath(storagedir, string(pdb.prefix) * "_" * chrom_name * ".bin"))
    end
    return ThreadsX.collect(x.prefix for x in output)
end


function getchromseq(isfa, record)
    return isfa ? FASTA.sequence(record) : TwoBit.sequence(record)
end