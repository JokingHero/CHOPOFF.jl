# My container that stores ambiguous guides (of same length)
# Its purpose is to determine efficiently how many exact matches
# of some guide are inside.
struct AmbigIdx
    ambig::SMatrix
    annot::Union{Vector{String}, Nothing}
end


function AmbigIdx(guides::Vector{<: BioSequence}, annot::Union{Vector{String}, Nothing})
    K = length(guides)
    T = length(guides[1])
    order = sortperm(guides)
    guides = vcat(collect.(guides[order])...)
    guides = SMatrix{T, K}(guides)
    if !isnothing(annot)
        return AmbigIdx(guides, annot[order])
    else
        return AmbigIdx(guides, annot)
    end
end


import Base.length
function length(idx::AmbigIdx)
    return Size(idx.ambig)[2]
end


"
Find for an input sequence how many sequences inside idx are compatible.
Returns a vector of bits, for subsetting on idx::annot.
"
function findbits(guide::LongDNA{4}, idx::AmbigIdx)
    bits = trues(length(idx))
    for (i, base) in enumerate(guide)
        bits[bits] = iscompatible.(idx.ambig[i, bits], base)
        if sum(bits) == 0
            return bits
        end
    end
    return bits
end