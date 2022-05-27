"""
`isinclusive(x::S, y::S) where {S<:BioSymbol}`

Check whether `x` contains all options of `y` or `y` contains all options of `x`.
Can be used as replacement in place of `BioSequences.iscompatible` or `BioSequences.isequal`.

# Examples
```jldoctest
julia> isinclusive(DNA_N, DNA_A)
true

julia> isinclusive(DNA_W, DNA_R) # W is (A/T), R is (A/G)
false

julia> isinclusive(DNA_W, DNA_A)
true
```
"""
@inline function isinclusive(x::S, y::S) where {S<:BioSymbol}
    if x == DNA_Gap || y == DNA_Gap 
        return false
    end
    return all(compatbits(x) & compatbits(y) == compatbits(x)) |
           all(compatbits(x) & compatbits(y) == compatbits(y))
end


@inline function commonprefix(
    guide::T,
    ref::K,
    ismatch::Function = iscompatible) where {T <: BioSequence, K <: BioSequence}
    l = 0
    for (i, g) in enumerate(guide)
        !ismatch(g, ref[i]) && break
        l += 1
    end
    return l
end


"""
`hamming(s1::T, s2::K, ismatch = iscompatible) where {T <: BioSequence, K <: BioSequence}`

Check Hamming distance between `s1` and `s2` using as matching definition `ismatch`.
Should be faster than `pairalign(HammingDistance(), s1, s2)`. Make sure inputs are of the 
same length.

# Examples
```jldoctest
julia> hamming(dna"ACGC", dna"AWRC")
1
```
"""
function hamming(s1::T, s2::K, ismatch = iscompatible) where {T <: BioSequence, K <: BioSequence}
    return count(!ismatch, s1, s2)
end


"""
```
levenshtein(guide::T, ref::K, k::Int = 4,
    ismatch::Function = iscompatible) where {T <: BioSequence, K <: BioSequence}
```

Calculate Levenshtein distance bounded by `k` maximum edit distance with a twist for guide + off-target comparisons.
    
Levenshtein distance is the minimum number of operations (consisting of insertions, deletions,
substitutions of a single character) required to change one string into the other.
`guide` input seqeunce is a gRNA sequence, `ref` input is reference sequence with expansion on the 3'
end of `k` bases. This extension will not count toward the score, if it is not covered with aligned guide.
Return k + 1, if distance is higher than k and terminate early. This function should be 10x faster than
`pairalign(LevenshteinDistance(), guide, ref, distance_only = true)`.

Notice that `guide` and `ref` have to be oriented towards left side e.g. PAM-guide, and PAM-offtarget-extension!
Take a look at the examples below to understand why.

# Examples
```jldoctest
julia> levenshtein(dna"ATGA", dna"AGACCT") # CCT as extension of the ref does not count towards the score in optimal alignment
1

julia> levenshtein(dna"ATGATCG", dna"AGAAATCGATG") # ATG does not count in optimal alignment ATG--ATCG/A-GAAATCG
3
```
"""
function levenshtein(
    guide::T,
    ref::K,
    k::Int = 4,
    ismatch::Function = iscompatible) where {T <: BioSequence, K <: BioSequence}

    len1, len2 = length(guide), length(ref)
    # prefix common to both strings can be ignored
    f = commonprefix(guide, ref, ismatch)
    f == len1 && return 0
    guide, ref = guide[f+1:len1], ref[f+1:len2]
    len1, len2 = len1 - f, len2 - f
    if k > len2
        k = len2
    end

    v = collect(1:len2)
    v_min = len2
    j_start_max = len2 - k
    if !(j_start_max > 0)
        j_start_max = 1
    end
    j_end = k

    @inbounds for (i, ch1) in enumerate(guide)
        # minimum v value for calculated rows - for early termination 
        # when we are outside of k dist
        v_min = len2

        j_start = max(1, i - k)
        if j_start > j_start_max
            j_start = j_start_max
        end
        if j_end < len2
            j_end += 1
        end

        top_left = j_start > 1 ? v[j_start - 1] : i - 1
        left = i # this in first iteration is always i at best
        @inbounds for j = j_start:j_end
            # compute future above
            current = top_left
            if !ismatch(ch1, ref[j])
                # mismatch when all options equal
                if current < v[j]
                    if current < left # mismatch
                        current += 1
                    else
                        current = left + 1
                    end
                else
                    if v[j] < left
                        current = v[j] + 1
                    else
                        current = left + 1
                    end
                end
            end

            if current < v_min
                v_min = current
            end

            top_left = v[j]
            left = current
            v[j] = current          
        end
        # we return v_min as we can "truncate" reference
        # which means we are interested in the best alignment 
        # for the lowest row
        v_min > k && return k + 1
    end
    return v_min # if we want full alignment without ref skip return v[end]
end


# same as above but additionally
# returns how many bases of the ref were NOT used!
# how many bases of ref are left out! can be 0 - can be more than extension
function levenshtein2(
    guide::T,
    ref::K,
    k::Int = 4,
    ismatch::Function = iscompatible) where {T <: BioSequence, K <: BioSequence}

    len1, len2 = length(guide), length(ref)
    # prefix common to both strings can be ignored
    f = commonprefix(guide, ref, ismatch)
    f == len1 && return 0, (len2 - len1)
    guide, ref = guide[f+1:len1], ref[f+1:len2]
    len1, len2 = len1 - f, len2 - f
    if k > len2
        k = len2
    end

    v = collect(1:len2)
    v_min = len2
    j_start_max = len2 - k
    if !(j_start_max > 0)
        j_start_max = 1
    end
    j_end = k

    @inbounds for (i, ch1) in enumerate(guide)
        # minimum v value for calculated rows - for early termination 
        # when we are outside of k dist
        v_min = len2

        j_start = max(1, i - k)
        if j_start > j_start_max
            j_start = j_start_max
        end
        if j_end < len2
            j_end += 1
        end

        top_left = j_start > 1 ? v[j_start - 1] : i - 1
        left = i # this in first iteration is always i at best
        @inbounds for j = j_start:j_end
            # compute future above
            current = top_left
            if !ismatch(ch1, ref[j])
                # mismatch when all options equal
                if current < v[j]
                    if current < left # mismatch
                        current += 1
                    else
                        current = left + 1
                    end
                else
                    if v[j] < left
                        current = v[j] + 1
                    else
                        current = left + 1
                    end
                end
            end

            if current < v_min
                v_min = current
            end

            top_left = v[j]
            left = current
            v[j] = current          
        end
        # we return v_min as we can "truncate" reference
        # which means we are interested in the best alignment 
        # for the lowest row
        v_min > k && return k + 1, (len2 - findlast(isequal(v_min), v))
    end
    return v_min, (len2 - findlast(isequal(v_min), v))
end


# left topleft top
@enum AlnPath ref_ nogap g_

"""
```
struct Aln
    guide::String
    ref::String
    dist::Int
end
```

Simple data structure to hold information on the optimal alignment. Therefore `guide` and
`ref` may contain `DNA_Gap` as these are aligned sequences. Function `align` returns this object.
"""
struct Aln
    guide::String
    ref::String
    dist::Int
end


"""
```
align(guide::T, ref::K, k::Int = 4,
    ismatch::Function = iscompatible) where {T <: BioSequence, K <: BioSequence}
```

Calculate Levenshtein distance bounded by `k` maximum edit distance with a twist for guide + off-target comparisons, but
return also the alignment as an `Align` object.

Levenshtein distance is the minimum number of operations (consisting of insertions, deletions,
substitutions of a single character) required to change one string into the other.
`guide` input seqeunce is a gRNA sequence, `ref` input is reference sequence with expansion on the 3'
end of `k` bases. This extension will not count toward the score, if it is not covered with aligned guide.
Return k + 1, if distance is higher than k and terminate early.

Notice that `guide` and `ref` have to be oriented towards left side e.g. PAM-guide, and PAM-offtarget-extension!
Take a look at the examples below to understand why.

# Examples
```jldoctest
julia> align(dna"ATGA", dna"AGACCT") # CCT as extension of the ref does not count towards the score in optimal alignment
Aln("ATGA", "A-GA", 1)

julia> align(dna"ATGATCG", dna"AGAAATCGATG") # ATG does not count
Aln("ATG--ATCG", "A-GAAATCG", 3)
```
"""
function align(
    guide::T,
    ref::K,
    k::Int = 4, # k has to be > 0
    ismatch::Function = iscompatible) where {T <: BioSequence, K <: BioSequence}

    len1, len2 = length(guide), length(ref)
    # prefix common to both strings can be ignored
    p = commonprefix(guide, ref, ismatch)
    p == len1 && return Aln(guide, ref[1:len1], 0)
    p_guide, p_ref = guide[p + 1:len1], ref[p + 1:len2]
    len1, len2 = len1 - p, len2 - p
    if k > len1
        k = len1
    end

    # we initialize all array with the top values
    v = Array(transpose(repeat(0:len1, 1, len2 + 1)))
    # when we are at the top right corner we have to go left
    # going left means here gap in the reference
    align = fill(ref_, (len2 + 1, len1 + 1))
    j_start_max = max(1, len1 - k)
    j_end = k

    # we will keep track of the smallest distance in the last column
    v_min_row, v_min_col = len1, len1
    v_min_col_idx = 1
    @inbounds for (i, ch1) in enumerate(p_ref)
        v_min_row = len1 # minimum score for each row
        j_start = max(1, i - k)
        if j_start > j_start_max
            j_start = j_start_max
        end
        if j_end < len1 
            j_end += 1
        end

        # sets top_left and left at the begining of the iteration
        #@info "$j_start $j_end"
        v[i + 1, j_start] = i
        # going top is gap in the guide
        align[i + 1, j_start] = g_

        @inbounds for j = j_start:j_end
            if !ismatch(ch1, p_guide[j])
                # mismatch when all options equal
                # top_left < top
                if v[i, j] < v[i, j + 1]
                    # top_left < left
                    if v[i, j] < v[i + 1, j] # mismatch
                        v[i + 1, j + 1] = v[i, j] + 1
                        align[i + 1, j + 1] = nogap
                    else # gap in ref
                        v[i + 1, j + 1] = v[i + 1, j] + 1
                        align[i + 1, j + 1] = ref_
                    end
                else
                    # top < left
                    if v[i, j + 1] < v[i + 1, j] # gap in guide
                        v[i + 1, j + 1] = v[i, j + 1] + 1
                        align[i + 1, j + 1] = g_
                    else # gap in ref
                        v[i + 1, j + 1] = v[i + 1, j] + 1
                        align[i + 1, j + 1] = ref_
                    end
                end
            else # match
                v[i + 1, j + 1] = v[i, j]
                align[i + 1, j + 1] = nogap
            end

            if v[i + 1, j + 1] < v_min_row
                v_min_row = v[i + 1, j + 1]
            end

            # we calculate for last column
            if j == len1 && v[i + 1, end] < v_min_col
                v_min_col = v[i + 1, end]
                v_min_col_idx = i + 1
            end
        end
        # v_min_row keeps smallest value for this row
        # we don't return anything yet as 
        # the smallest is in the last column not row!
        v_min_row > k && break
    end
    v_min_col > k && return Aln("", "", k + 1)

    j = len1 + 1 # j along guide, i along reference
    aln_g, aln_ref = "", ""
    while v_min_col_idx > 1 || j > 1
        if align[v_min_col_idx, j] == nogap
            aln_g *= string(p_guide[j - 1])
            aln_ref *= string(p_ref[v_min_col_idx - 1])
            v_min_col_idx -= 1
            j -= 1
        elseif align[v_min_col_idx, j] == g_
            aln_g *= "-"
            aln_ref *= string(p_ref[v_min_col_idx - 1])
            v_min_col_idx -= 1
        else
            aln_g *= string(p_guide[j - 1])
            aln_ref *= "-"
            j -= 1
        end
    end
    return Aln(string(guide[1:p]) * reverse(aln_g), 
               string(ref[1:p]) * reverse(aln_ref), 
               v_min_col)
end


## PREFIX & SUFFIX PARTIAL ALIGNMENT
struct PrefixAlignment{T <: BioSequence, K <: BioSequence}
    v::Array{Int64, 2}
    align::Array{AlnPath, 2}
    prefix::T
    prefix_len::Int
    guide::K
    guide_len::Int
    suffix_len::Int
    v_min_col::Int
    v_min_col_idx::Int
    k::Int
    isfinal::Bool # true when no more alignment is needed
end

function prefix_align(
    guide::T,
    prefix::K,
    suffix_len::Int,
    k::Int = 4, # asumption that k has to be less than prefix length
    ismatch::Function = iscompatible) where {T <: BioSequence, K <: BioSequence}

    guide_len, prefix_len = length(guide), length(prefix)
    ref_len = prefix_len + suffix_len 
    if k > prefix_len
        k = prefix_len
    end

    # we initialize all array with the top values
    v = Array(transpose(repeat(0:guide_len, 1, ref_len + 1)))
    # when we are at the top right corner we have to go left
    # going left means here gap in the reference
    align = fill(ref_, (ref_len + 1, guide_len + 1))
    j_start_max = max(1, guide_len - k)
    j_end = k

    # we will keep track of the smallest distance in the last column
    v_min_row, v_min_col = guide_len, guide_len
    v_min_col_idx = 1
    @inbounds for (i, ch1) in enumerate(prefix)
        v_min_row = guide_len # minimum score for each row
        j_start = max(1, i - k)
        if j_start > j_start_max
            j_start = j_start_max
        end
        if j_end < guide_len 
            j_end += 1
        end

        # sets top_left and left at the begining of the iteration
        v[i + 1, j_start] = i
        # going top is gap in the guide
        align[i + 1, j_start] = g_

        @inbounds for j = j_start:j_end
            if !ismatch(ch1, guide[j])
                # mismatch when all options equal
                # top_left < top
                if v[i, j] < v[i, j + 1]
                    # top_left < left
                    if v[i, j] < v[i + 1, j] # mismatch
                        v[i + 1, j + 1] = v[i, j] + 1
                        align[i + 1, j + 1] = nogap
                    else # gap in ref
                        v[i + 1, j + 1] = v[i + 1, j] + 1
                        align[i + 1, j + 1] = ref_
                    end
                else
                    # top < left
                    if v[i, j + 1] < v[i + 1, j] # gap in guide
                        v[i + 1, j + 1] = v[i, j + 1] + 1
                        align[i + 1, j + 1] = g_
                    else # gap in ref
                        v[i + 1, j + 1] = v[i + 1, j] + 1
                        align[i + 1, j + 1] = ref_
                    end
                end
            else # match
                v[i + 1, j + 1] = v[i, j]
                align[i + 1, j + 1] = nogap
            end

            if v[i + 1, j + 1] < v_min_row
                v_min_row = v[i + 1, j + 1]
            end

            # we calculate for last column
            if j == guide_len && v[i + 1, end] < v_min_col
                v_min_col = v[i + 1, end]
                v_min_col_idx = i + 1
            end
        end
        v_min_row > k && break
    end
    return PrefixAlignment(v, align,
                           prefix, prefix_len, 
                           guide, guide_len, 
                           suffix_len, 
                           v_min_col, v_min_col_idx, k, v_min_row > k)
end


function suffix_align(
    suffix::K,
    pA::PrefixAlignment,
    ismatch::Function = iscompatible) where {K <: BioSequence}

    # we initialize all array with the top values
    v = copy(pA.v)
    align = copy(pA.align)
    j_start_max = max(1, pA.guide_len - pA.k)
    j_end = min(pA.k + pA.prefix_len, pA.guide_len)

    # we will keep track of the smallest distance in the last column
    v_min_row = pA.guide_len
    v_min_col = pA.v_min_col
    v_min_col_idx = pA.v_min_col_idx
    if !pA.isfinal
        @inbounds for (l, ch1) in enumerate(suffix)
            i = l + pA.prefix_len
            v_min_row = pA.guide_len # minimum score for each row
            j_start = max(1, i - pA.k)
            if j_start > j_start_max
                j_start = j_start_max
            end
            if j_end < pA.guide_len 
                j_end += 1
            end

            # sets top_left and left at the begining of the iteration
            v[i + 1, j_start] = i
            # going top is gap in the guide
            align[i + 1, j_start] = g_
            @inbounds for j = j_start:j_end
                if !ismatch(ch1, pA.guide[j])
                    # mismatch when all options equal
                    # top_left < top
                    if v[i, j] < v[i, j + 1]
                        # top_left < left
                        if v[i, j] < v[i + 1, j] # mismatch
                            v[i + 1, j + 1] = v[i, j] + 1
                            align[i + 1, j + 1] = nogap
                        else # gap in ref
                            v[i + 1, j + 1] = v[i + 1, j] + 1
                            align[i + 1, j + 1] = ref_
                        end
                    else
                        # top < left
                        if v[i, j + 1] < v[i + 1, j] # gap in guide
                            v[i + 1, j + 1] = v[i, j + 1] + 1
                            align[i + 1, j + 1] = g_
                        else # gap in ref
                            v[i + 1, j + 1] = v[i + 1, j] + 1
                            align[i + 1, j + 1] = ref_
                        end
                    end
                else # match
                    v[i + 1, j + 1] = v[i, j]
                    align[i + 1, j + 1] = nogap
                end

                if v[i + 1, j + 1] < v_min_row
                    v_min_row = v[i + 1, j + 1]
                end

                # we calculate for last column
                if j == pA.guide_len && v[i + 1, end] < v_min_col
                    v_min_col = v[i + 1, end]
                    v_min_col_idx = i + 1
                end
            end
            # v_min_row keeps smallest value for this row
            # we don't return anything yet as 
            # the smallest is in the last column not row!
            v_min_row > pA.k && break
        end
    end
    if v_min_col > pA.k
        return Aln("", "", pA.k + 1)
    end

    j = pA.guide_len + 1 # j along guide, i along reference
    aln_g, aln_ref = "", ""
    ref = string(pA.prefix) * string(suffix)
    while v_min_col_idx > 1 || j > 1
        if align[v_min_col_idx, j] == nogap
            aln_g *= string(pA.guide[j - 1])
            aln_ref *= string(ref[v_min_col_idx - 1])
            v_min_col_idx -= 1
            j -= 1
        elseif align[v_min_col_idx, j] == g_
            aln_g *= "-"
            aln_ref *= string(ref[v_min_col_idx - 1])
            v_min_col_idx -= 1
        else
            aln_g *= string(pA.guide[j - 1])
            aln_ref *= "-"
            j -= 1
        end
    end
    return Aln(reverse(aln_g), reverse(aln_ref), v_min_col)
end


"
Mutating pA constantly, restarts alignment from `start_from`.
`start_from` should be relative to the suffix.
"
function suffix_align!(
    suffix::K,
    pA::PrefixAlignment,
    start_from::Int = 1,
    ismatch::Function = iscompatible) where {K <: BioSequence}

    # we initialize all array with the top values
    v = pA.v
    align = pA.align
    
    j_start_max = max(1, pA.guide_len - pA.k)
    ref_start_from = pA.prefix_len + start_from
    j_end = min(pA.k + ref_start_from - 1, pA.guide_len)

    # we will keep track of the smallest distance in the last column
    v_min_row = pA.guide_len
    v_min_col, v_min_col_idx = findmin(v[1:ref_start_from, end])
    if !pA.isfinal
        @inbounds for l in start_from:pA.suffix_len
            i = l + pA.prefix_len
            v_min_row = pA.guide_len # minimum score for each row
            j_start = max(1, i - pA.k)
            if j_start > j_start_max
                j_start = j_start_max
            end
            if j_end < pA.guide_len 
                j_end += 1
            end

            # sets top_left and left at the begining of the iteration
            v[i + 1, j_start] = i
            # going top is gap in the guide
            align[i + 1, j_start] = g_
            @inbounds for j = j_start:j_end
                if !ismatch(suffix[l], pA.guide[j])
                    # mismatch when all options equal
                    # top_left < top
                    if v[i, j] < v[i, j + 1]
                        # top_left < left
                        if v[i, j] < v[i + 1, j] # mismatch
                            v[i + 1, j + 1] = v[i, j] + 1
                            align[i + 1, j + 1] = nogap
                        else # gap in ref
                            v[i + 1, j + 1] = v[i + 1, j] + 1
                            align[i + 1, j + 1] = ref_
                        end
                    else
                        # top < left
                        if v[i, j + 1] < v[i + 1, j] # gap in guide
                            v[i + 1, j + 1] = v[i, j + 1] + 1
                            align[i + 1, j + 1] = g_
                        else # gap in ref
                            v[i + 1, j + 1] = v[i + 1, j] + 1
                            align[i + 1, j + 1] = ref_
                        end
                    end
                else # match
                    v[i + 1, j + 1] = v[i, j]
                    align[i + 1, j + 1] = nogap
                end

                if v[i + 1, j + 1] < v_min_row
                    v_min_row = v[i + 1, j + 1]
                end

                # we calculate for last column
                if j == pA.guide_len && v[i + 1, end] < v_min_col
                    v_min_col = v[i + 1, end]
                    v_min_col_idx = i + 1
                end
            end
            # v_min_row keeps smallest value for this row
            # we don't return anything yet as 
            # the smallest is in the last column not row!
            v_min_row > pA.k && break
        end
    end
    if v_min_col > pA.k
        return Aln("", "", pA.k + 1)
    end

    j = pA.guide_len + 1 # j along guide, i along reference
    aln_g, aln_ref = "", ""
    ref = string(pA.prefix) * string(suffix)
    while v_min_col_idx > 1 || j > 1
        if align[v_min_col_idx, j] == nogap
            aln_g *= string(pA.guide[j - 1])
            aln_ref *= string(ref[v_min_col_idx - 1])
            v_min_col_idx -= 1
            j -= 1
        elseif align[v_min_col_idx, j] == g_
            aln_g *= "-"
            aln_ref *= string(ref[v_min_col_idx - 1])
            v_min_col_idx -= 1
        else
            aln_g *= string(pA.guide[j - 1])
            aln_ref *= "-"
            j -= 1
        end
    end
    return Aln(reverse(aln_g), reverse(aln_ref), v_min_col)
end


"
Created for testing purposes, performs prefix alignment followed
by suffix alignment.
"
function pa_sa(guide::LongDNA{4}, ref::LongDNA{4}, d::Int, prefix_len::Int)
    if prefix_len > length(guide)
        prefix_len = length(guide)
    end
    prefix = ref[1:prefix_len]
    suffix = ref[prefix_len+1:end]
    pa = prefix_align(guide, prefix, length(suffix), d)
    return suffix_align(suffix, pa)
end


## BITPARALEL

"
Compute levenshtein distance.

Input: guide (p), ref (t)
Returns: true if it can find an alignment where score <= k

Four times faster than dynamic programing version, but
does not work for our guides as it will ignore left extra bases of the ref
Example of two valid alignments where this implementation would fail:
   AACTG
 AAAAC GTGAAA
 -- these As should be counted toward the score and are not
    (score is 1, position is 6)

    AAC TG
  AAAACGTGAAA
  -- also skipped As (score is 1, position is 8)

 No idea how to resolve this problem.

Based on Mayers bit-parallel alghoritm. Some references:
# Some examples and bandend version when n > 64
# http://www.mi.fu-berlin.de/wiki/pub/ABI/AdvancedAlgorithms11_Searching/script-03-ShiftOrUkkonenBitVecMyers.pdf
# https://research.ijcaonline.org/volume72/number14/pxc3889214.pdf
# http://www.mi.fu-berlin.de/wiki/pub/ABI/RnaSeqP4/myers-bitvector-verification.pdf
# https://dl.acm.org/doi/10.1145/316542.316550 or https://www.win.tue.nl/~jfg/educ/bit.mat.pdf
# transposition -> file:///home/ai/Downloads/10.1.1.19.7158.pdf
"
function levenshtein_bp(guide::String, ref::String, k::Int = 4)
    m = length(guide)
    n = length(ref)
    @assert m <= 64
    bit_alphabet = zeros(UInt64, 4)

    @inbounds for i in 1:m
        bit_alphabet[base_to_idx(guide[i])] |= UInt64(1) << (i - 1)
    end

    vp = (UInt64(1) << m) - UInt64(1)
    vn = UInt64(0)
    score = m

    @inbounds for pos in 1:n
        x = bit_alphabet[base_to_idx(ref[pos])] | vn
        d0 = ((vp + (x & vp)) ⊻ vp) | x

        hn = vp & d0
        hp = vn | ~(vp | d0)

        if hp & (UInt64(1) << (m - 1)) != 0
            score += 1
        elseif hn & (UInt64(1) << (m - 1)) != 0
            score -= 1
        end

        x = hp << 1
        vn = x & d0
        vp = (hn << 1) | ~(x | d0)

        if score <= k
            return true
        end
    end
    return false
end
