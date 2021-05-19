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


function is_within_d(s::LongDNASeq, x::LongDNASeq, d::Int)
    dist = levenshtein(s, x, d)
    if dist == d
        return true

    # add corner cases where we can have e.g. d = 1
    # GGN ACTGA 
    # GGNGACTGA
    # which gives reference guide of GACTG 
    # which is normally dist 2        ACTGA
    # but is valid as reference also has A there, unfortunatelly
    # we have to assume it might be the case and we count those off-targets
    elseif dist - 1 == d
        return levenshtein(s[1:end-1], x, d) == d
    end
    return false
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
    comb = comb_of_d1(s, alphabet)
    for i in 1:(d-1)
        new_comb = Set{String}()
        for el in comb
            new_comb = union(new_comb, comb_of_d1(el, alphabet))
        end
        comb = new_comb
    end

    comb = collect(comb)
    dist = [is_within_d(LongDNASeq(s), LongDNASeq(x), d) for x in comb]
    return comb[dist]
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