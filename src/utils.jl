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



"
This is a helper function, it can generate distances larger than 1!
"
function comb_of_d1(s::String, alphabet::Vector{Char} = ['A', 'C', 'T', 'G'])
    s_ = collect(s)
    allcomb = Set{String}()
    idx_in_s = combinations(1:length(s), 1)
    alphabet_ = vcat(alphabet, ['-'])
    for i in idx_in_s
        for j in alphabet_
            if s_[i[1]] != j
                if j == '-'
                    scopy_new = copy(s_)
                    deleteat!(scopy_new, i[1])
                    for k in alphabet
                        # gap in the s -> we insert base at the index and truncate to the size
                        scopy_s_ = copy(s_)
                        insert!(scopy_s_, i[1], k)
                        push!(allcomb, join(scopy_s_[1:length(s)]))
                        # gap in the new string -> we delete base at the index and add base at the end
                        scopy_new_ = copy(scopy_new)
                        append!(scopy_new_, k)
                        push!(allcomb, join(scopy_new_))
                        # same as above but add base at the begining
                        scopy_new_ = copy(scopy_new)
                        prepend!(scopy_new_, k)
                        push!(allcomb, join(scopy_new_))
                    end
                else
                    scopy = copy(s_)
                    scopy[i[1]] = j
                    push!(allcomb, join(scopy))
                end
            end
        end
    end
    # for now we need this - suboptimal method
    #allcomb = collect(allcomb)
    #dist = [levenshtein(LongDNA{4}(s), LongDNA{4}(x), 1) == 1 for x in allcomb]
    return allcomb #[dist]
end


# this is only working for distance 1, and it ignores
# bulges on the reference, as these can be also covered by mismatches 
function comb_of_d1(s::LongDNA{4})
    alphabet = [DNA_A, DNA_C, DNA_T, DNA_G] # we can't deal with ambig at this moment here
    allcomb = Set{LongDNA{4}}()
    allcomb1 = Set{LongDNA{4}}()
    for (i, si)  in enumerate(s)
        for ai in alphabet
            if isequal(si, ai) # skip a base and add all combinations
                s_ = copy(s)
                deleteat!(s_, i)
                push!(s_, DNA_A) 
                push!(allcomb, copy(s_))
                s_[end] = DNA_C
                push!(allcomb, copy(s_))
                s_[end] = DNA_T
                push!(allcomb, copy(s_))
                s_[end] = DNA_G
                push!(allcomb, copy(s_))
            else 
                # make all possible 1 mismatch
                s_ = copy(s)
                s_[i] = ai
                push!(allcomb, s_)
            end
            # bulge on the reference! - size +1
            s_ = copy(s)
            insert!(s_, i, ai)
            push!(allcomb1, s_)
        end
    end
    # filter redundant cases that are covered by mismatch
    allcomb1 = collect(allcomb1)
    allcomb1 = filter(x -> !(x[1:end-1] in allcomb), allcomb1)
    return setdiff(allcomb, Set([s])), allcomb1
end


# same as above, but not truncating longer distances!
function comb_of_d1_extended(s::String, alphabet::Vector{Char} = ['A', 'C', 'T', 'G'])
    s_ = collect(s)
    allcomb = Set{String}()
    idx_in_s = combinations(1:length(s), 1)
    alphabet_ = vcat(alphabet, ['-'])
    for i in idx_in_s
        for j in alphabet_
            if s_[i[1]] != j
                if j == '-'
                    scopy_new = copy(s_)
                    deleteat!(scopy_new, i[1])
                    for k in alphabet
                        # gap in the s -> we insert base at the index and truncate to the size
                        scopy_s_ = copy(s_)
                        insert!(scopy_s_, i[1], k)
                        push!(allcomb, join(scopy_s_))
                        # gap in the new string -> we delete base at the index and add base at the end
                        scopy_new_ = copy(scopy_new)
                        append!(scopy_new_, k)
                        push!(allcomb, join(scopy_new_))
                        # same as above but add base at the begining
                        scopy_new_ = copy(scopy_new)
                        prepend!(scopy_new_, k)
                        push!(allcomb, join(scopy_new_))
                    end
                else
                    scopy = copy(s_)
                    scopy[i[1]] = j
                    push!(allcomb, join(scopy))
                    deleteat!(scopy, i[1])
                    push!(allcomb, join(scopy))
                end
            end
        end
    end
    return allcomb
end


# This is returning all length(s) + 1 size possible off-targets
# from the perspective of the reference
function comb_of_d1_extended_ref(s::String, alphabet::Vector{Char} = ['A', 'C', 'T', 'G'])
    s_ = collect(s)
    allcomb = Set{String}()
    idx_in_s = 1:length(s_)
    bulge_guide = collect(permutations(alphabet, 2))
    # bulge_guide = filter(x -> s_[end] != x[1], bulge_guide)
    alphabet_ = vcat(alphabet, ['-'])
    for i in idx_in_s
        for j in alphabet_
            if s_[i] != j
                if j == '-'
                    scopy = copy(s_)
                    deleteat!(scopy, i) # bulge on the guide
                    for k in bulge_guide
                        push!(allcomb, join(vcat(scopy, k)))
                    end
                else
                    scopy = copy(s_)
                    scopy[i] = j # mismatch + extension on the reference
                    for k in alphabet
                        push!(allcomb, join(vcat(scopy, k)))
                    end

                    scopy = copy(s_)
                    insert!(scopy, i, j) # bulge on the reference
                    push!(allcomb, join(scopy))
                end
            end
        end
    end
    return allcomb
end


"
0 - is not within distance d
1 - is within distance d
2 - is within distance d, but assuming bulge and that guide
last position will match on the ref (which we don't know)

d is assumed to be > 0
"
function is_within_d(s::LongDNA{4}, x::LongDNA{4}, d::Int)
    dist = levenshtein(s, x, d)
    if dist == d
        return 1
    # add corner cases where we can have e.g. d = 1
    # GGN ACTGA 
    # GGNGACTGA
    # which gives reference guide of GACTG 
    # which is normally dist 2        ACTGA
    # but is valid as reference also has A there, unfortunatelly
    # we have to assume it might be the case and we count those off-targets
    elseif dist - 1 == d
        if levenshtein(s[1:end-1], x, d) == d
            return 2
        end
    end
    return 0
end


"
Create a list of possible strings of levenshtein distance d
toward the string s. Don't include combinations 
which are smaller and larger than d, but also include 
corner cases.
Assume PAM inside s is on the left!
'-' in alphabet will be treated as indel, don't use it.
"
function comb_of_d(s::String, d::Int = 1, alphabet::Vector{Char} = ['A', 'C', 'T', 'G'])
    if d == 0 
        return ([s], [])
    end
    comb = comb_of_d1(s, alphabet)
    for i in 1:(d-1)
        comb = foldxt(union, Map(x -> comb_of_d1(x, alphabet)), comb)
    end

    comb = ThreadsX.collect(comb)
    dist = ThreadsX.collect(is_within_d(LongDNA{4}(s), LongDNA{4}(x), d) for x in comb)
    return (comb[dist .== 1], comb[dist .== 2])
end


"""
`minkmersize(len::Int = 20, d::Int = 4)`

Pidgeon hole principle: minimum
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

Transform DNA to UInt128, usefull for hashing, and space saving.
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


@inline function BioSequences.LongDNA{4}(x::UInt64, len::Int)
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
    for comb in Iterators.product(amb_dna...)
        for (idx, dna) in zip(amb_idx, comb)
            res[i][idx] = dna
        end
        i += 1
    end
    return res
end


"
This simplistic strategy seems to compress around 50% more compared to
super fast sort. Selection of the starting index does not seem to influence
the compression greatly.
"
function order_by_hamming_and_prefix(guides::Vector{LongDNA{4}}, i::Int64 = 1)
    guides_len = length(guides)
    final_order = zeros(Int64, guides_len)
    is_done = zeros(Bool, guides_len)
    all_done = 1
    final_order[all_done] = i
    
    while all_done < guides_len
        is_done[i] = 1
        g = guides[i]
        g_h = ThreadsX.map(1:guides_len) do x 
            if is_done[x]
                return 0
            else
                return ARTEMIS.hamming(guides[x], g)
            end
        end
    
        g_h_min = Vector{Int64}()
        for w in 1:length(g)
            g_h_min = ThreadsX.findall(x -> x == w, g_h)
            if length(g_h_min) > 0 
                break
            end
        end
    
        prefix_len = ThreadsX.map(x -> ARTEMIS.commonprefix(x, g), guides[g_h_min])
        i = g_h_min[argmax(prefix_len)]
        all_done += 1
        final_order[all_done] = i
    end
    
    return final_order
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
# make a temporary directory
tdir = tempname()
ldb_path = joinpath(tdir, "linearDB")
mkpath(ldb_path)

# use ARTEMIS example genome
genome = joinpath(
    vcat(
        splitpath(dirname(pathof(ARTEMIS)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))

# finally, build a linearDB
build_linearDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0), 
    ldb_path)
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
function summarize_offtargets(res::DataFrame, distance::Int)
    df = groupby(res, :guide)
    df = combine(df,  :distance => (x -> counts_by_dist(x, distance)) => AsTable)
    sort!(df, vcat([Symbol("D$i") for i in 0:distance], :guide))
    return df
end