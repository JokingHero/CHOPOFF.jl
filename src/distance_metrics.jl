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

    if occursin("-", string(guide))
        error("Guide sequence ", string(guide), "contains -.")
    elseif occursin("-", string(ref))
        error("Reference sequence ", string(ref), "contains -.")
    end

    len1, len2 = length(guide), length(ref)
    # prefix common to both strings can be ignored
    f = commonprefix(guide, ref, ismatch)
    f == len1 && return 0
    guide, ref = guide[f+1:len1], ref[f+1:len2]
    len1 = length(guide)
    if k >= len1
        k = len1
    end

    # large loop on the reference restricted by deletions
    # small loop on the guide restricted by k
    v = collect(1:len1)
    current = 0
    v_min_idx = 0
    @inbounds for (i, ch1) in enumerate(ref)
        left = current = v_min = i - 1
        @inbounds for j = max(1, i - k):min(len1, i + k)
            above, current, left = current, left, v[j]
            if !ismatch(ch1, guide[j])
                # mismatch when all options equal
                if current < above
                    if current < left # mismatch
                        current = current + 1
                    else # gap in guide
                        current = left + 1
                    end
                else
                    if above < left # gap in ref
                        current = above + 1
                    else # mismatch
                        current = left + 1
                    end
                end
            end

            if v_min < left
                v_min_idx = j - 1
            else
                v_min = left
                v_min_idx = j
            end
            v[j] = current
        end
        # we aligned all of guide - only ref is left
        # return smallest distance in this row
        v_min_idx == len1 && return v_min
        v_min > k && return k + 1
    end
    current > k && return k + 1
    return current
end


## PREFIX & SUFFIX PARTIAL ALIGNMENT
struct PrefixAlignment
    v::Vector{Int}
    prefixlen::Int
    dist::Int
    isfinal::Bool # true when no more alignment is needed
end

function prefix_levenshtein(
    guide::T,
    prefix::K,
    k::Int = 4,
    ismatch::Function = isinclusive) where {T <: BioSequence, K <: BioSequence}
    len1, len2 = length(guide), length(prefix)
    if k >= len1
        k = len1
    end

    # large loop on the reference restricted by deletions
    # small loop on the guide restricted by k
    v = collect(1:len1)
    current = 0
    v_min_idx = 0
    @inbounds for (i, ch1) in enumerate(prefix)
        println("Done $i")
        left = current = v_min = i - 1
        @inbounds for j = max(1, i - k):min(len1, i + k)
            above, current, left = current, left, v[j]
            if !ismatch(ch1, guide[j])
                # mismatch when all options equal
                if current < above
                    if current < left # mismatch
                        current = current + 1
                    else # gap in guide
                        current = left + 1
                    end
                else
                    if above < left # gap in prefix
                        current = above + 1
                    else # mismatch
                        current = left + 1
                    end
                end
            end

            if v_min < left
                v_min_idx = j - 1
            else
                v_min = left
                v_min_idx = j
            end
            v[j] = current
        end
        # we aligned all of guide - only prefix is left
        # return smallest distance in this row
        v_min_idx == len1 && return PrefixAlignment(v, len2, v_min, true)
        v_min > k && return PrefixAlignment(v, len2, k + 1, true)
    end
    return PrefixAlignment(v, len2, current, false)
end


function suffix_levenshtein(
    guide::T,
    suffix::K,
    pA::PrefixAlignment,
    k::Int = 4,
    ismatch::Function = isinclusive) where {T <: BioSequence, K <: BioSequence}

    len1, len2 = length(guide), length(suffix)
    if k >= len1
        k = len1
    end
    # large loop on the reference restricted by deletions
    # small loop on the guide restricted by k
    v = copy(pA.v)
    current = pA.dist
    v_min_idx = 0
    @inbounds for i in (pA.prefixlen + 1):(len2 + pA.prefixlen)
        left = current = v_min = i - 1
        @inbounds for j = max(1, i - k):min(len1, i + k)
            above, current, left = current, left, v[j]
            if !ismatch(suffix[i - pA.prefixlen], guide[j])
                # mismatch when all options equal
                if current < above
                    if current < left # mismatch
                        current = current + 1
                    else # gap in guide
                        current = left + 1
                    end
                else
                    if above < left # gap in suffix
                        current = above + 1
                    else # mismatch
                        current = left + 1
                    end
                end
            end

            if v_min < left
                v_min_idx = j - 1
            else
                v_min = left
                v_min_idx = j
            end
            v[j] = current
        end
        # we aligned all of guide - only suffix is left
        # return smallest distance in this row
        v_min_idx == len1 && return v_min
        v_min > k && return k + 1
    end
    current > k && return k + 1
    return current
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
