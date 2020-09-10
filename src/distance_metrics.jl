#
# using BioSymbols
# using Random
# using BioSequences
# using BioAlignments
# using BenchmarkTools

"
x contains all options of y
or
y contains all options of x

e.g.
N      , anything -> True
W (A/T), R (A/G)  -> False
W      , A        -> True
"
@inline function isinclusive(x::S, y::S) where S <: BioSymbol
	return all(compatbits(x) & compatbits(y) == compatbits(x)) |
	 all(compatbits(x) & compatbits(y) == compatbits(y))
end


@inline function commonprefix(guide::NucleotideSeq, ref::NucleotideSeq,
	ismatch::Function=isinclusive)
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
function hamming(s1::NucleotideSeq, s2::NucleotideSeq, ismatch = isinclusive)
	return count(!ismatch, s1, s2)
end


"""
Levenshtein distance with a twist for guide + off-target comparisons
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
(consisting of insertions, deletions, substitutions of a single character) required to change one string into the other.
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
function levenshtein(guide::NucleotideSeq, ref::NucleotideSeq,
	k::Integer = 4, ismatch::Function=isinclusive)

	len1, len2 = length(guide), length(ref)
	# prefix common to both strings can be ignored
    f = commonprefix(guide, ref, ismatch)
    f == len1 && return 0
	guide, ref = guide[f+1:end], ref[f+1:end]
	len1 = length(guide)
	if k >= len1
		k = len1
	end

	# large loop on the reference restricted by deletions
	# small loop on the guide restricted by k
    v = collect(1:len1)
	v_min_idx = 0
	current = 0
    for (i, ch1) in enumerate(ref)
        left = current = v_min = i - 1
        for j in max(1, i - k):min(len1, i + k)
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

# some tests
# Random.seed!(42)
# function getSeq(N = 20, letters = ['A', 'C', 'G', 'T'])
# 	return LongDNASeq(randstring(letters, N))
# end
#
# s, t = getSeq(), getSeq(24)
# println(s)
# println(t)
#
# aln = pairalign(LevenshteinDistance(), s, t)
# println(aln)
# levenshtein(s, t, 4)
