
using BioSequences
using BioAlignments
using CRISPRofftargetHunter
using Edlib

using BenchmarkTools

# N guides
# make distance 4 sequence
seq = getSeq()
ref_seq = getSeq(1) * seq[1:16] * getSeq(3) * seq[17:end]

m1 = median(@benchmark pairalign(LevenshteinDistance(), g, r, distance_only = true) setup=(g=seq; r=ref_seq))
m2 = median(@benchmark levenshtein(g, r) setup=(g=seq; r=ref_seq))
m3 = median(@benchmark Edlib.edit_distance(g, r, max_distance = 4, mode = :prefix) setup=(g=string(seq); r=string(ref_seq)))
m4 = median(@benchmark levenshtein_bp(g, r) setup=(g=string(seq); r=string(ref_seq)))

show(m1)
show(m2)
show(m3)
show(m4) # 4 times faster! wow

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
