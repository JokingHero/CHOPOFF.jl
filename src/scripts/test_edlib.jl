
using BioSequences
using BioAlignments
using CRISPRofftargetHunter
using Edlib

using BenchmarkTools

# N guides
# make distance 4 sequence
seq = getSeq()
ref_seq = getSeq(1) * seq[1:16] * getSeq(3) * seq[17:end]

function withk(g::String, r::String)
    g = getkgrams(g, 5)
    r = getkgrams(r, 5)
    return !isdisjoint(g, r)
end

withk(string(seq), string(ref_seq))

m1 = median(@benchmark pairalign(LevenshteinDistance(), g, r, distance_only = true) setup=(g=seq; r=ref_seq))
m2 = median(@benchmark levenshtein(g, r) setup=(g=seq; r=ref_seq))
m3 = median(@benchmark Edlib.edit_distance(g, r, max_distance = 4, mode = :prefix) setup=(g=string(seq); r=string(ref_seq)))
m4 = median(@benchmark levenshtein_bp(g, r) setup=(g=string(seq); r=string(ref_seq)))
m5 = median(@benchmark isdisjoint(g, r) setup=(g=getkmers(string(seq)); r=getkmers(string(ref_seq))))
m6 = median(@benchmark withk(g, r) setup=(g=string(seq); r=string(ref_seq)))

show(m1)
show(m2)
show(m3)
show(m4) # 4 times faster! wow
show(m5)
show(m6) # much slower

TP = 139655
TN = 933042475
FP = 2114660545
FN = 0

# what will be the numbers for this with qgram?!

# Is it better to ignore pidgeon hole speedup?
(m6.time * (FP + TN) + FP * m4.time) * 1e-9 # ns
(m4.time * (FP + TN)) * 1e-9 # ns

!isdisjoint(getkmers(string(seq)), getkmers(string(ref_seq)))
!isdisjoint(getkmers(string(getSeq())), getkmers(string(getSeq())))

# what % of kmer overlaps is actually inside the distance?


@benchmark levenshtein(setSeq.g, setSeq.r, 5)
@benchmark Edlib.edit_distance(setSeqStr.g, setSeqStr.r, max_distance = 5, mode = :prefix)

guides = map(x -> getSeq(20), 1:100000)
ref = getSeq(24)

function test(guides, ref)
    return map(g -> CRISPRofftargetHunter.levenshtein(g, ref, 4), guides)
end

test(guides[1:2], ref)

@benchmark test($guides, $ref)
#0.87s
#0.22s - 4 times as fast?!

using FASTX
f = open(FASTA.Writer, "../../../Soft/edlib/test_g.fa")
for (i, x) in enumerate(guides)
    write(f, FASTA.Record(string(i), x))
end
close(f)
