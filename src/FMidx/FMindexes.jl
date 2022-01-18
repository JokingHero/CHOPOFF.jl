__precompile__()

module FMIndexes

using SuffixArrays
#using WaveletMatrices # this is not on registry too ffs!
include("WaveletMatrices.jl")
using BioSequences
using .WaveletMatrices
using IndexableBitVectors
using Humanize
using Mmap

include("saca.jl")


"""
Index for full-text search.

Type parameters:

* `w`: the number of bits required to encode the alphabet
* `T`: the type to represent positions of a sequence
"""
struct FMIndex{w,T}
    bwt::WaveletMatrix{w,UInt8,SucVector}
    sentinel::Int
    samples::Vector{T}
    sampled::SucVector
    count::Vector{Int}
end

function FMIndex(seq, sa, σ, r)
    wm = WaveletMatrix(make_bwt(seq, sa), log2(Int, σ))
    # sample suffix array
    samples, sampled = sample_sa(sa, r)
    sentinel = something(findfirst(isequal(0), sa), 0) + 1
    # count characters
    count = count_bytes(seq, σ)
    count[1] = 1  # sentinel '$' is smaller than any character
    cumsum!(count, count)
    return FMIndex(wm, sentinel, samples, SucVector(sampled), count)
end

"""
    FMIndex(seq, σ=256; r=32, program=:SuffixArrays, mmap::Bool=false, opts...)

Build an FM-Index from a sequence `seq`.
The sequence must support `convert(UInt8, seq[i])` for each character and the
alphabet size should be less than or equal to 256. The second parameter, `σ`, is
the alphabet size. The third parameter, `r`, is the interval of sampling values
from a suffix array. If you set it large, you can save the memory footprint but
it requires more time to locate the position.
"""
function FMIndex(seq, σ=256; r=32, program=:SuffixArrays, mmap::Bool=false, opts...)
    T = index_type(length(seq))
    opts = Dict(opts)
    @assert !haskey(opts, :σ) "σ should be passed as the second argument"
    if program === :SuffixArrays
        @assert 1 ≤ σ ≤ typemax(UInt8) + 1
        sa = make_sa(T, seq, σ, mmap)
    elseif program === :psascan || program === :pSAscan
        @assert 1 ≤ σ ≤ typemax(UInt8)
        psascan = get(opts, :psascan, "psascan")
        workdir = get(opts, :workdir, pwd())
        sa = make_sa_pscan(T, seq, psascan, workdir, mmap)
    else
        error("unknown program name: $program")
    end
    return FMIndex(seq, sa, σ, r)
end

"""
    FMIndex(text; opts...)

Build an FM-Index from an ASCII text.
"""
function FMIndex(text::Union{String, SubString{String}}; opts...)
    return FMIndex(codeunits(text), 128; opts...)
end

Base.length(index::FMIndex) = length(index.bwt)

function Base.show(io::IO, fmindex::FMIndex{w,T}) where {w,T}
    println(io, summary(fmindex), ':')
    totalsize = (
        sizeof(fmindex.bwt) +
        sizeof(fmindex.samples) +
        sizeof(fmindex.sampled) +
        sizeof(fmindex.count)
    )
    print("     length: ", length(fmindex), '\n')
    print("  data size: ", Humanize.datasize(totalsize, style=:bin))
end

"""
Restore the original text from the index.
"""
function restore(index::FMIndex)
    n = length(index)
    text = Vector{UInt8}(undef, n)
    p = index.sentinel
    while n > 0
        p = lfmap(index, p)
        text[n] = index.bwt[p ≥ index.sentinel ? p - 1 : p]
        n -= 1
    end
    return text
end

function log2(::Type{Int}, x::Integer)
    return 64 - leading_zeros(convert(UInt64, x-1))
end

function count_bytes(seq, σ)
    count = zeros(Int, σ + 1)
    for i in 1:length(seq)
        count[convert(UInt8,seq[i])+2] += 1
    end
    resize!(count, σ)
    return count
end

# LF-mapping
function lfmap(index::FMIndex, i)
    if i == index.sentinel
        return 1
    elseif i > index.sentinel
        i -= 1
    end
    char = index.bwt[i]
    @inbounds return index.count[char+1] + rank(char, index.bwt, i)
end

function sa_range(query, index::FMIndex)
    sa_range(query, index::FMIndex, 1:(length(index)+1))
end

function sa_range(query, index::FMIndex, init_range::UnitRange{Int})
    sp, ep = init_range.start, init_range.stop
    # backward search
    i = length(query)
    while sp ≤ ep && i ≥ 1
        char = convert(UInt8, query[i])
        c = index.count[char+1]
        sp = c + rank(char, index.bwt, (sp > index.sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + rank(char, index.bwt, (ep > index.sentinel ? ep - 1 : ep))
        i -= 1
    end
    return sp:ep
end

@inline function sa_range(char::UInt8, index::FMIndex, range::UnitRange{Int})
    sp, ep = range.start, range.stop
    c = index.count[char+1]
    sp = c + rank(char, index.bwt, (sp > index.sentinel ? sp - 1 : sp) - 1) + 1
    ep = c + rank(char, index.bwt, (ep > index.sentinel ? ep - 1 : ep))
    return sp:ep
end

function sa_range!(query::LongDNASeq, index::FMIndex, cashe::Dict{LongDNASeq, UnitRange{Int}})
    sp, ep = 1, (length(index) + 1)
    cashe_val = sp:ep
    i = length(query)
    while !isnothing(cashe_val) && i ≥ 1
        cashe_val = get(cashe, LongDNASeq(query[i:end]), nothing)
        if !isnothing(cashe_val)
            sp = cashe_val.start
            ep = cashe_val.stop
            i -= 1
        end
    end

    # backward search
    while sp ≤ ep && i ≥ 1
        char = convert(UInt8, query[i])
        c = index.count[char+1]
        sp = c + rank(char, index.bwt, (sp > index.sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + rank(char, index.bwt, (ep > index.sentinel ? ep - 1 : ep))
        cashe[LongDNASeq(query[i:end])] = sp:ep
        i -= 1
    end
    return sp:ep
end

function sa_value(i::Int, index::FMIndex)
    if i == 1
        # point to the sentinel '$'
        return length(index) + 1
    end
    d = 0
    @inbounds while !index.sampled[i-1]
        i = lfmap(index, i)
        d += 1
    end
    return index.samples[rank1(index.sampled, i - 1)] + d
end

sa_value(i::Integer, index::FMIndex) = sa_value(Int(i), index)

"""
Count the number of occurrences of the given query.
"""
function Base.count(query, index::FMIndex)
    return length(sa_range(query, index))
end

"""
Count the number of occurrences of the given query.
Use cashe. Update cashe with each iteration.
"""
function count_with_cashe!(query::LongDNASeq, index::FMIndex, cashe::Dict{LongDNASeq, UnitRange{Int}})
    return length(sa_range!(query, index, cashe))
end

struct LocationIterator{w,T}
    range::UnitRange{Int}
    index::FMIndex{w,T}
end

Base.length(iter::LocationIterator) = length(iter.range)

function Base.iterate(iter::LocationIterator, i::Int=1)
    if i > length(iter)
        return nothing
    end
    return sa_value(iter.range[i], iter.index) + 1, i + 1
end

"""
Locate the positions of occurrences of the given query.
This method returns an iterator of positions:

    for pos in locate(query, index)
        # ...
    end
"""
function locate(query, index::FMIndex)
    return LocationIterator(sa_range(query, index), index)
end

"""
Locate the positions of all occurrences of the given query.
"""
function locateall(query, index::FMIndex)
    iter = locate(query, index)
    locs = Vector{Int}(undef, length(iter))
    for (i, loc) in enumerate(iter)
        locs[i] = loc
    end
    return locs
end


export
    FMIndex,
    restore,
    count,
    count_with_cashe!,
    locate,
    locateall

end # module