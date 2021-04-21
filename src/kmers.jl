"
Return ordered `k`-mers based on a string `seq`.
"
function getkmers(seq::String, k::Int = 4)
    return [seq[i:i+k-1] for i in 1:(length(seq)-k+1)]
end


"
Return a set of `k`-grams based on a string `seq`.
These will be k non-overlapping consecutive slices of `seq`.
"
function getkgrams(seq::String, k::Int = 5)
    len = Int(floor(length(seq)/k))
    kgrams = Set([seq[1:len], seq[len * (k - 1) + 1:end]])
    for i in 2:(k-1)
        push!(kgrams, seq[len * (i - 1) + 1:len*i])
    end
    return kgrams
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
Appends occurence number to the kmer.
"
function kmer_with_order(x)
    order = zeros(Int, length(x))
    x_ = Set(x)
    if length(x_) == length(x)
        return [xi * "_" * string(order[i]) for (i, xi) in enumerate(x)]
    end
    for k in x_
        k_idx = findall(isequal(k), x)
        order[k_idx] = collect(0:(length(k_idx) - 1))
    end
    return [xi * "_" * string(order[i]) for (i, xi) in enumerate(x)]
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
    end
end


"
Generate all possible k-mers.
"
function getkmers(k::Int = 5; alphabet::Vector{Char} = ['A', 'C', 'T', 'G'])
    return join.(collect(multiset_permutations(repeat(alphabet, k), k)))
end



"
Stores all kmers in the data and all the counts of possible transitions
between kmers along the positions of the original string.
"
struct KmerTransitions # TODO replace UInt16?!
    kmers::Vector{String} # all possible kmers
    transitions_by_pos::Array{UInt16, 3} # dimensions: kmers, ACTG, transition position
end


function Base.getindex(kt::KmerTransitions, from_kmer::String, to_kmer::String, pos::Int)
    return kt.transitions_by_pos[findfirst(isequal(from_kmer), kt.kmers), base_to_idx(to_kmer[end]), pos]
end


function Base.getindex(kt::KmerTransitions, s::String)
    kmer_len = length(first(kt.kmers))
    kmers = getkmers(s, kmer_len)
    kt_min = kt[kmers[1], kmers[2], 1]
    for i in 2:(length(kmers) - 1)
        kt_min_i = kt[kmers[i], kmers[i + 1], i]
        if kt_min > kt_min_i
            kt_min = kt_min_i
        end
    end
    return kt_min
end


function add!(kt::KmerTransitions, from_kmer::String, to_kmer::String, pos::Int, count::Int)
    kmer_idx = findfirst(isequal(from_kmer), kt.kmers) # TODO replace with dict for speed?!
    base_idx = base_to_idx(to_kmer[end])
    @inbounds old_val = kt.transitions_by_pos[kmer_idx, base_idx, pos]
    @inbounds kt.transitions_by_pos[kmer_idx, base_idx, pos] = 
        safeadd(old_val, convert(eltype(kt.transitions_by_pos), count))
    return nothing
end


### Kmer based measures - probably useless
"
Jaccard index
"
function jcd(x, y)
    x_ = Set(x)
    y_ = Set(y)
    return length(x_ ∩ y_) / length(x_ ∪ y_)
end


"
Weighted Jaccard index
"
function wjcd(x, y)
    all_kmers = Set(vcat(x, y))
    wjcd_min = 0
    wjcd_max = 0
    for k in all_kmers
        k_in_x = count(isequal(k), x)
        k_in_y = count(isequal(k), y)
        if k_in_x >= k_in_y
            wjcd_min += k_in_y
            wjcd_max += k_in_x
        else
            wjcd_min += k_in_x
            wjcd_max += k_in_y
        end
    end
    return wjcd_min / wjcd_max
end


"
Ordered Weighted Jaccard index.
Locality-sensitive hashing for the edit distance. Marcais et al. 2019
"
function owjcd(x, y)
    # kmer + kmer occurence
    x_ = kmer_with_order(x)
    y_ = kmer_with_order(y)
    return wjcd(x_, y_)
end