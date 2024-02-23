" Don't overflow the typemax. "
safeadd(x::T, y::T) where {T} = ifelse(x + y ≥ x, x + y, typemax(T))


"""
`getseq(n = 20, letters = ['A', 'C', 'G', 'T'])`

Randomize sequence of length `n` from `letters`.
"""
function getseq(n = 20, letters = ['A', 'C', 'G', 'T'])
    return LongDNA{4}(randstring(letters, n))
end


function base_to_idx(letter::Char)
    if letter == 'A'
        return 1
    elseif letter == 'C'
        return 2
    elseif letter == 'T'
        return 3
    elseif letter == 'G'
        return 4
    else
        throw("Not ACTG character.")
    end
end


"
Get file extension from the string path `s`.
"
function extension(s::String)
    extension = match(r"\.[A-Za-z0-9]+$", s)
    if extension !== nothing
        return extension.match
    else
        return ""
    end
end


"
Returns smallest possible Unsigned type that can contain
given `max_value`.
"
function smallestutype(max_value::Unsigned)
    if typemax(UInt8) >= max_value
        return UInt8
    elseif typemax(UInt16) >= max_value
        return UInt16
    elseif typemax(UInt32) >= max_value
        return UInt32
    elseif typemax(UInt64) >= max_value
        return UInt64
    elseif typemax(UInt128) >= max_value
        return UInt128
    else
        throw("Too big unsigned value to fit in our types.")
    end
end


"""
`minkmersize(len::Int = 20, d::Int = 4)`

Pigeon hole principle: minimum
k-mer size that is required for two strings of
size `len` to be aligned within distance of `d`.

# Examples
```jldoctest
julia> minkmersize(20, 3)
5

julia> minkmersize(20, 6)
2
```
"""
function minkmersize(len::Int = 20, d::Int = 4)
    return Int(floor(len / (d + 1)))
end


"
    Will try to find value in `x` that will allow for almost equal
    split of values into buckets/leafs.
"
function balance(x::Vector{Int})
    if isempty(x)
        return nothing
    end
    uniq = unique(x)
    sort!(uniq)
    counts = [count(y -> y == i, x) for i in uniq]
    balance = argmin(abs.([sum(counts[1:i]) - sum(counts[i:end]) for i = 1:length(counts)]))
    return uniq[argmin(abs.(uniq .- uniq[balance]))]
end


"
Instead of this there is also BigMER!

Transform DNA to UInt128, useful for hashing, and space saving.
We zero bits that are beyond sequence length. No clue why these bits are set.
"
function Base.convert(::Type{UInt128}, x::LongDNA{4})
    gap_idx = findfirst(DNA_Gap, x)
    len = length(x)
    x = x.data
    if !isnothing(gap_idx)
        throw("DNA with gaps can't be converted!")
    end
    if length(x) == 2
        b = 4 * (16 - (len - 16))
        return parse(UInt128, bitstring((x[2] << b) >> b) * bitstring(x[1]); base = 2)
    elseif length(x) == 1
        b = 4 * (16 - len)
        return parse(UInt128, bitstring(UInt64(0)) * bitstring((x[1] << b) >> b); base = 2)
    else
        throw("Too long DNA sequence, it has to be less than 31 of length.")
    end
end


# Slow!
@inline function Base.convert(::Type{LongDNA{4}}, x::UInt128)
    x = bitstring(x)
    x_seq = LongDNA{4}("")
    @inbounds for i in reverse(1:4:length(x))
        xi = reinterpret(DNA, parse(UInt8, x[i:i + 3]; base = 2))
        if xi != DNA_Gap
            push!(x_seq, xi)
        else
            break # GAP is not allowed
        end
    end
    return x_seq
end


function Base.convert(::Type{UInt64}, x::LongDNA{4})
    # BioSequences.encoded_data(DNAMer(x))
    x = LongDNA{2}(x)
    if (length(x.data) > 1) 
        throw("Sequence too long to save as UInt64.")
    end

    y = zero(UInt64)
    for c in x
        nt = convert(DNA, c)
        if isambiguous(nt)
            throw(ArgumentError("cannot create a mer with ambiguous nucleotides"))
        elseif isgap(nt)
            throw(ArgumentError("cannot create a mer with gaps"))
        end
        y = (y << 2) | UInt64(BioSequences.twobitnucs[reinterpret(UInt8, nt) + 0x01])
    end

    mask = (one(UInt64) << (2 * length(x))) - 1
    return reinterpret(UInt64, y & mask) # encoded_data
end


function Base.convert(::Type{UInt32}, x::LongDNA{4})
    x = LongDNA{2}(x)
    if (length(x) > 16) 
        throw("Sequence too long to save as UInt32.")
    end

    y = zero(UInt32)
    for c in x
        nt = convert(DNA, c)
        if isambiguous(nt)
            throw(ArgumentError("cannot create a mer with ambiguous nucleotides"))
        elseif isgap(nt)
            throw(ArgumentError("cannot create a mer with gaps"))
        end
        y = (y << 2) | UInt32(BioSequences.twobitnucs[reinterpret(UInt8, nt) + 0x01])
    end

    mask = (one(UInt32) << (2 * length(x))) - one(UInt32)
    return reinterpret(UInt32, y & mask) # encoded_data
end


# this is needed inside saca.jl
# TODO make sure the queried sequences are from the same type
# or conform to the same encoding as the reference
# aka. I am pretty sure it won't work {2} vs {4}
@inline function Base.convert(::Type{UInt8}, x::DNA)
    return reinterpret(UInt8, x)
end

@inline function Base.convert(::Type{String}, x::BioSequence)
    return String(x)
end


@inline function BioSequences.LongDNA{4}(x::UInt128, len::Int)
    data = [convert(UInt64, (x << 64) >> 64), convert(UInt64, x >> 64)]
    return LongDNA{4}(data, UInt64(len))
end


@inline function BioSequences.LongDNA{4}(x::Union{UInt64, UInt32}, len::Int)
    y = []
    for i in 1:len
        push!(y, reinterpret(DNA, 0x01 << ((x >> 2(len - i)) & 0b11)))
    end
    
    # LongDNA{4}(DNAMer{len}(x))
    #return LongDNA{4}([nt for nt in x])
    #return LongDNA{4}(LongDNA{2}([x], UInt64(len)))
    return LongDNA{4}(y)
end


"This can be used to asses whether query is inside the guide."
@inline function is_equal(guide::UInt64, query::UInt64, shift::Int64)
    return count_ones(~((guide >> shift) ⊻ query)) == 64
end


import BioSymbols.iscertain
function BioSymbols.iscertain(x::LongDNA{4})
    return all(iscertain.(x))
end


function isambig(x::LongDNA{4})
    return !iscertain(x)
end


# assumes elements are sorted A->C->G->T
const TO_AMBIGUOUS = IdDict(
    "ACGT" => DNA_N,

    "CGT" => DNA_B,
    "AGT" => DNA_D,
    "ACG" => DNA_V,
    "ACT" => DNA_H,

    "GT" => DNA_K,
    "AG" => DNA_R,
    "AC" => DNA_M,
    "AT" => DNA_W,
    "CG" => DNA_S,
    "CT" => DNA_Y,

    "A" => DNA_A,
    "C" => DNA_C,
    "G" => DNA_G,
    "T" => DNA_T,

    "N" => DNA_N,
    )


const FROM_AMBIGUOUS = IdDict(
    DNA_N => [DNA_A, DNA_C, DNA_T, DNA_G],

    DNA_B => [       DNA_C, DNA_T, DNA_G],
    DNA_D => [DNA_A,        DNA_T, DNA_G],
    DNA_V => [DNA_A, DNA_C,        DNA_G],
    DNA_H => [DNA_A, DNA_C, DNA_T,      ],

    DNA_K => [              DNA_T, DNA_G],
    DNA_R => [DNA_A,               DNA_G],
    DNA_M => [DNA_A, DNA_C,             ],
    DNA_W => [DNA_A,        DNA_T,      ],
    DNA_S => [       DNA_C,        DNA_G],
    DNA_Y => [       DNA_C, DNA_T,      ],
    )


"""
`expand_ambiguous(x::LongDNA{4})`

Turn ambiguous bases e.g. N, into all possible combinations.

# Examples
```jldoctest
julia> expand_ambiguous(dna"AN")
4-element Vector{LongSequence{DNAAlphabet{4}}}:
 AA
 AC
 AT
 AG
 ```
"""
function expand_ambiguous(x::LongDNA{4})
    amb_dna = Vector{Vector{DNA}}()
    amb_idx = Vector{Int64}()
    for (i, dna) in enumerate(x)
        if isambiguous(dna)
            push!(amb_idx, i)
            push!(amb_dna, FROM_AMBIGUOUS[dna])
        end
    end
    iter = Iterators.product(amb_dna...)
    res = [copy(x) for i in 1:length(iter)]
    i = 1
    for comb in iter
        for (idx, dna) in zip(amb_idx, comb)
            res[i][idx] = dna
        end
        i += 1
    end
    return res
end


"""
`all_kmers(size = 4; alphabet = [DNA_A, DNA_C, DNA_G, DNA_T]`

Make a list of all possible kmers with given`size` using bases in the `alphabet`. 

# Examples
```jldoctest
julia> all_kmers(2; alphabet = [DNA_A, DNA_N])
4-element Vector{LongSequence{DNAAlphabet{4}}}:
 AA
 AN
 NA
 NN
```
"""
function all_kmers(size = 4; alphabet = [DNA_A, DNA_C, DNA_G, DNA_T])
    iters = Iterators.product(ntuple(_ -> alphabet, size)...)
    iters = LongDNA{4}.(collect(iters)[:])
    return sort(iters)
end


function as_bitvector_of_kmers(x::LongDNA{4}, kmers::Dict{LongDNA{4}, Int})
    BioSequences.ungap!(x)
    bits = zeros(length(kmers))
    kmer_size = length(first(first(kmers)))
    for i in 1:(length(x) - kmer_size + 1)
        xi = x[i:(i + kmer_size - 1)]
        if iscertain(xi)
            bits[kmers[xi]] = 1
        else # replace all ambigs with nonambigs
            xexp = expand_ambiguous(xi)
            for xi in xexp
                bits[kmers[xi]] = 1
            end
        end
    end
    return BitVector(bits)
end



"""
`as_kmers(x::LongDNA{4}, kmer_size::Int)`

Transforms `x` into vector of kmers of size `kmer_size`. 
All ambiguous bases will be expanded.

# Examples
```jldoctest
julia> as_kmers(dna"ACTGG", 4)
2-element Vector{LongSequence{DNAAlphabet{4}}}:
 ACTG
 CTGG
```
"""
function as_kmers(x::LongDNA{4}, kmer_size::Int)
    kmers = Vector{LongDNA{4}}()
    for i in 1:(length(x) - kmer_size + 1)
        xi = x[i:(i + kmer_size - 1)]
        if iscertain(xi)
            push!(kmers, xi)
        else # replace all ambigs with nonambigs
            xexp = expand_ambiguous(xi)
            for xi in xexp
                push!(kmers, xi)
            end
        end
    end
    return kmers
end


"""
`as_skipkmers(x::LongDNA{4}, kmer_size::Int)`

Transforms `x` into vector of skip-kmers of size `kmer_size`. 
All ambiguous bases will be expanded. 
Leftover-bases are ignored!

# Examples
```jldoctest
julia> as_skipkmers(dna"ACTGG", 2)
2-element Vector{LongSequence{DNAAlphabet{4}}}:
 AC
 TG
```
"""
function as_skipkmers(x::LongDNA{4}, kmer_size::Int)
    kmers = Vector{LongDNA{4}}()
    for i in 1:kmer_size:(length(x) - kmer_size + 1)
        xi = x[i:(i + kmer_size - 1)]
        if iscertain(xi)
            push!(kmers, xi)
        else # replace all ambigs with nonambigs
            xexp = expand_ambiguous(xi)
            for xi in xexp
                push!(kmers, xi)
            end
        end
    end
    return kmers
end


function format_DF(res::Matrix{Int64}, dist::Int, guides::Vector{LongDNA{4}})
    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, vcat(col_d, :guide))
    select!(res, :guide, Not(:guide)) 
    return res
end


"""
```
filter_overlapping(res::DataFrame, distance::Int)
```

Filter overlapping off-targets. Remember that off-targets have their start relative to the PAM location.

# Arguments
`res` - DataFrame created by one of the off-target finding methods, it contains columns 
    such as `:guide,  :chromosome, :strand, :distance, :start`. 

`distance` - To what distance from the `:start` do we consider the off-target to be overlapping?

# Examples
```julia-repl
$(make_example_doc())
```
"""
function filter_overlapping(res::DataFrame, distance::Int)
    sort!(res, [:guide,  :chromosome, :strand, :distance, :start])

    to_remove = falses(size(res)[1])
    gsc = join(res[1, [:guide, :chromosome, :strand]])
    pos_filter = Set((res.start[1] - distance):(res.start[1] + distance))
    for i in 2:length(to_remove)
        gsci = join(res[i, [:guide, :chromosome, :strand]])
        if gsc == gsci
            if res.start[i] in pos_filter
                to_remove[i] = true
            else
                for j in (res.start[i] - distance):(res.start[i] + distance)
                    push!(pos_filter, j)
                end
            end
        else
            gsc = gsci
            pos_filter = Set((res.start[i] - distance):(res.start[i] + distance))
        end
    end
    res = res[.!to_remove, :]
end


# small helper that is similar to R table function
function counts_by_dist(distances, max_dist)
    map = countmap(distances)
    vec = zeros(Int, max_dist + 1)
    for (dist, count) in map
        vec[dist + 1] = count
    end
    return NamedTuple{Tuple([Symbol("D$i") for i in 0:max_dist])}(Tuple(vec))
end


"""
```
summarize_offtargets(res::DataFrame; distance::Int = maximum(res.distance))
```

Summarize all off-targets into count table from the detail file. This does not automatically filters overlaps.
You can specify distance to filter out some of the higher distances.

# Arguments
`res` - DataFrame created by one of the off-target finding methods, it contains columns 
    such as `:guide,  :chromosome, :strand, :distance, :start`. 

`distance` - What is the maximum distance to assume in the data frame, 
    its possible to specify smaller distance than contained in the `res` DataFrame and autofilter lower distances.

# Examples
```julia-repl
$(make_example_doc())
```
"""
function summarize_offtargets(res::DataFrame; distance::Int = maximum(res.distance))
    df = groupby(res, :guide)
    df = combine(df,  :distance => (x -> counts_by_dist(x, distance)) => AsTable)
    sort!(df, vcat([Symbol("D$i") for i in 0:distance], :guide))
    return df
end