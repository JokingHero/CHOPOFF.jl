"
x contains all options of y
or
y contains all options of x

e.g.
N      , anything -> True
W (A/T), R (A/G)  -> False
W      , A        -> True
"
@inline function isinclusive(x::S, y::S) where {S<:BioSymbol}
    return all(compatbits(x) & compatbits(y) == compatbits(x)) |
           all(compatbits(x) & compatbits(y) == compatbits(y))
end


@inline function commonprefix(
    guide::T,
    ref::K,
    ismatch::Function = isinclusive) where {T <: BioSequence, K <: BioSequence}
    l = 0
    for (i, g) in enumerate(guide)
        !ismatch(g, ref[i]) && break
        l += 1
    end
    return l
end


"
Hamming distance

Using BioSequences is faster than other implementations. Using LongDNASeq is
faster than DNAmer and allows for all IUPAC codes. This is faster than
pairalign(HammingDistance(), s, t)
and faster than
pairalign(EditDistance(), s3, t3, CostModel(match=0, mismatch=1, insertion=1, deletion=1))

For all human genome we have 304794990 guides NGG, median hamming distance speed
is 24e-9 seconds -> 7.3 second per query on full linear scan for all
off-targets.
"
function hamming(s1::T, s2::K, ismatch = isinclusive) where {T <: BioSequence, K <: BioSequence}
    return count(!ismatch, s1, s2)
end


"""
Levenshtein distance with a twist for guide + off-target comparisons.
Keeps the triangle inequality principle with bounded by k maximum edit
distance.

guide ATGA
ref   A-GACCT
Score: 1 (CCT as extension of the genomic reference
does not count towards the score)

guide ATGA--TCG
ref   A-GAAATCGATG
Score: 3 (ATG does not count)

Levenshtein distance is the minimum number of operations
(consisting of insertions, deletions, substitutions of a single character)
required to change one string into the other.
Left input seqeunce is guide sequence, right input is reference sequence
with expansion on the 3' end of max_dist bases. This extension will not
count toward the score if it is not covered with aligned guide.

Return k + 1 if distance higher than k and terminate early.
Medium benchmark on random 20bp guide vs 24bp ref:

Random.seed!(42)
function getSeq(N = 20, letters = ['A', 'C', 'G', 'T'])
	return LongDNASeq(randstring(letters, N))
end

setSeq =(g=getSeq();r=getSeq(24))
m1 = median(@benchmark pairalign(LevenshteinDistance(), g, r, distance_only = true) setup=setSeq)
m2 = median(@benchmark levenshtein(g, r) setup=setSeq)

BenchmarkTools.TrialEstimate:
  time:             225.705 ns
  gctime:           0.000 ns (0.00%)
  memory:           336 bytes
  allocs:           3

print(judge(m1, m2))
TrialJudgement(+1126.73% => regression)
"""
function levenshtein(
    guide::T,
    ref::K,
    k::Int = 4,
    ismatch::Function = isinclusive) where {T <: BioSequence, K <: BioSequence}

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
    offset_k = k - (len2 - len1)
    j_start = 1
    j_start_max = len2 - k
    j_end = k

    @inbounds for (i, ch1) in enumerate(guide)
        # minimum v value for calculated rows - for early termination 
        # when we are outside of k dist
        v_min = len2
        v_min_idx = len2

        if (i - 1 > offset_k && j_start < j_start_max)
            j_start += 1
        end
        if j_end < len2
            j_end += 1
        end

        top_left = i - 1
        left = i

        @info "$j_start $j_end"

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
        #v_min > k && return k + 1
    end
    # ??
    #v_min > k && return k + 1
    @info "$v"
    return v_min # if we want full alignment without ref skip return v[end]
end


# left topleft top
@enum AlnPath ref_ nogap g_

struct Aln
    # can contain "-"
    guide::String
    ref::String
    k::Int
    dist::Int
end

"
Same as above, BUT 
v is along the guide, not the reference,
large loop is along reference
small loop along guide.

Also return alignment.
"
function levenshtein2(
    guide::T,
    ref::K,
    k::Int = 4,
    ismatch::Function = isinclusive) where {T <: BioSequence, K <: BioSequence}

    len1, len2 = length(guide), length(ref)
    # prefix common to both strings can be ignored
    p = commonprefix(guide, ref, ismatch)
    p == len1 && return Aln(guide, ref[1:len1], k, 0)
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

    v_min = len1
    for (i, ch1) in enumerate(p_ref) # add @inbounds
        v_min = len1 # minimum score for each row
        j_start = max(1, i - k)
        if j_start > j_start_max
            j_start = j_start_max
        end
        j_end = i + k # k + 1 on the right
        if j_end > len1 
            j_end = len1
        end

        # sets top_left and left at the begining of the iteration
        #@info "$j_start $j_end"
        v[i + 1, j_start] = i
        # going top is gap in the guide
        align[i + 1, j_start] = g_

        for j =  j_start:j_end # @in_bounds
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

            if v[i + 1, j + 1] < v_min
                v_min = v[i + 1, j + 1]
            end
        end
        # v_min keeps smallest value for this row
        v_min > k && break
    end
    
    # find smallest value in the last column, 
    # we can truncate reference to minimize edit dist
    # exclude _ref gaps in this column!
    min_dist = v[end, end]
    i = len2 + 1
    for (idx, dist) in enumerate(v[:, end])
        if dist < min_dist && align[idx, end] != ref_
            min_dist = dist
            i = idx
        end
    end
    min_dist > k && return Aln("", "", k, k + 1)

    j = len1 + 1 # j along guide, i along reference
    aln_g, aln_ref = "", ""
    while i > 1 || j > 1
        if align[i, j] == nogap
            aln_g *= string(p_guide[j - 1])
            aln_ref *= string(p_ref[i - 1])
            i -= 1
            j -= 1
        elseif align[i, j] == g_
            aln_g *= "-"
            aln_ref *= string(p_ref[i - 1])
            i -= 1
        else
            aln_g *= string(p_guide[j - 1])
            aln_ref *= "-"
            j -= 1
        end
    end
    return Aln(string(guide[1:p]) * reverse(aln_g), 
               string(ref[1:p]) * reverse(aln_ref), 
               k, min_dist)
end


## PREFIX & SUFFIX PARTIAL ALIGNMENT
struct PrefixAlignment
    v::Vector{Int}
    prefixlen::Int
    dist::Int
    k::Int
    isfinal::Bool # true when no more alignment is needed
end

function prefix_levenshtein(
    guide::T,
    prefix::K,
    suffix_len::Int,
    k::Int = 4,
    ismatch::Function = isinclusive) where {T <: BioSequence, K <: BioSequence}

    len1, len2 = length(guide), length(prefix) + suffix_len
    # prefix common to both strings can be ignored
    f = commonprefix(guide, prefix, ismatch)
    f == len1 && return 0
    guide, prefix = guide[f+1:len1], prefix[f+1:len2]
    len1, len2 = len1 - f, len2 - f
    if k > len2
        k = len2
    end

    v = collect(1:len2)
    v_min = len2
    offset_k = k - (len2 - len1) # for speed instead of max(1, i - k)
    j_start = 1
    j_end = k

    @inbounds for (i, ch1) in enumerate(guide)
        # minimum v value for calculated rows - for early termination 
        # when we are outside of k dist
        v_min = len2
        v_min_idx = len2

        if (i - 1 > offset_k)
            j_start += 1
        end
        if j_end < len2
            j_end += 1
        end

        top_left = i - 1
        left = i

        @inbounds for j = j_start:j_end
            # compute future above
            current = top_left
            if !ismatch(ch1, prefix[j])
                # mismatch when all options equal
                if current < v[j]
                    if current < left # mismatch
                        current += 1
                    else # gap in guide
                        current = left + 1
                    end
                else
                    if v[j] < left # gap in prefix
                        current = v[j] + 1
                    else # mismatch
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
        v_min > k && return PrefixAlignment(v, len2, k + 1, k, true)
    end

    v_min > k && return PrefixAlignment(v, len2, k + 1, k, true)
    return PrefixAlignment(v, len2, v[end], k, false)
end


function suffix_levenshtein(
    guide::T,
    suffix::K,
    pA::PrefixAlignment,
    k::Int = 4,
    ismatch::Function = isinclusive) where {T <: BioSequence, K <: BioSequence}

    len1, len2 = length(guide), length(suffix)
    if k >= len2
        k = len2
    end
    # large loop on the reference restricted by deletions
    # small loop on the guide restricted by k
    v = copy(pA.v)
    current = pA.dist
    v_min_idx = 0
    #@inbounds for i in (pA.prefixlen + 1):(len2 + pA.prefixlen)
    #    left = current = v_min = i - 1
    #    @inbounds for j = max(1, i - k - 1):min(len1, i + k)
    #        above, current, left = current, left, v[j]
    
    
    len1, len2 = length(guide), length(suffix)
    if k > len2
        k = len2
    end

    v = collect(1:len2)
    v_min = len2
    offset_k = k - (len2 - len1) # for speed instead of max(1, i - k)
    j_start = 1
    j_end = k

    @inbounds for (i, ch1) in enumerate(guide)
        # minimum v value for calculated rows - for early termination 
        # when we are outside of k dist
        v_min = len2
        v_min_idx = len2

        if (i - 1 > offset_k)
            j_start += 1
        end
        if j_end < len2
            j_end += 1
        end

        top_left = i - 1
        left = i

        @inbounds for j = j_start:j_end
            # compute future above
            current = top_left
            if !ismatch(ch1, ref[j]) !ismatch(suffix[i - pA.prefixlen], ref[j])
                # mismatch when all options equal
                if current < v[j]
                    if current < left # mismatch
                        current += 1
                    else # gap in guide
                        current = left + 1
                    end
                else
                    if v[j] < left # gap in ref
                        current = v[j] + 1
                    else # mismatch
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
        v_min > k && return k + 1
    end
    # we return v_min as we can "truncate" reference
    # which means we are interested in the best alignment 
    # for the lowest row
    v_min > k && return k + 1
    return v_min # if we want full alignment without ref skip return v[end]
end


## BITPARALEL

function get_idx(letter::Char)
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
        bit_alphabet[get_idx(guide[i])] |= UInt64(1) << (i - 1)
    end

    vp = (UInt64(1) << m) - UInt64(1)
    vn = UInt64(0)
    score = m

    @inbounds for pos in 1:n
        x = bit_alphabet[get_idx(ref[pos])] | vn
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
