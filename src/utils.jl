" Don't overflow the typemax. "
safeadd(x::T, y::T) where {T} = ifelse(x + y ≥ x, x + y, typemax(T))


"
Randomize sequence of length `n` from `letters`.
"
function getseq(n = 20, letters = ['A', 'C', 'G', 'T'])
    return LongDNASeq(randstring(letters, n))
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
    #dist = [levenshtein(LongDNASeq(s), LongDNASeq(x), 1) == 1 for x in allcomb]
    return allcomb #[dist]
end


"
0 - is not within distance d
1 - is within distance d
2 - is within distance d, but assuming bulge and that guide
last position will match on the ref (which we don't know)

d is assumed to be > 0
"
function is_within_d(s::LongDNASeq, x::LongDNASeq, d::Int)
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
    dist = ThreadsX.collect(is_within_d(LongDNASeq(s), LongDNASeq(x), d) for x in comb)
    return (comb[dist .== 1], comb[dist .== 2])
end


"
Pidgeon hole principle: minimum
k-mer size that is required for two strings of
size `len` to be aligned within distance of `d`.
"
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
function Base.convert(::Type{UInt128}, x::LongDNASeq)
    gap_idx = findfirst(DNA_Gap, x)
    len = length(x)
    x = BioSequences.encoded_data(x)
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
@inline function Base.convert(::Type{LongDNASeq}, x::UInt128)
    x = bitstring(x)
    x_seq = LongDNASeq("")
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


@inline function BioSequences.LongDNASeq(x::UInt128, len::Int)
    y = LongDNASeq(len)
    y.data = [convert(UInt64, (x << 64) >> 64), convert(UInt64, x >> 64)]
    return y
end


function isambig(x::LongDNASeq)
    return count(isambiguous, x) > 0
end


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


function expand_ambiguous(x::LongDNASeq)
    amb_dna = Vector{Vector{DNA}}()
    amb_idx = Vector{Int64}()
    for (i, dna) in each(isambiguous, x)
        push!(amb_idx, i)
        push!(amb_dna, FROM_AMBIGUOUS[dna])
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
function order_by_hamming_and_prefix(guides::Vector{LongDNASeq}, i::Int64 = 1)
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
                return CRISPRofftargetHunter.hamming(guides[x], g)
            end
        end
    
        g_h_min = Vector{Int64}()
        for w in 1:length(g)
            g_h_min = ThreadsX.findall(x -> x == w, g_h)
            if length(g_h_min) > 0 
                break
            end
        end
    
        prefix_len = ThreadsX.map(x -> CRISPRofftargetHunter.commonprefix(x, g), guides[g_h_min])
        i = g_h_min[argmax(prefix_len)]
        all_done += 1
        final_order[all_done] = i
    end
    
    return final_order
end