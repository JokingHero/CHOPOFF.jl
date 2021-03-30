using CRISPRofftargetHunter
using BioSequences
using BenchmarkTools

path = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..","test", "sample_data", "sample_genome.fa")
path = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..","test", "sample_data", "semirandom.fa")
#path = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"

# should not be a part of dbi
db_path = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..","test", "sample_db", "linearDB")
#db_path = "/home/ai/Projects/uib/crispr/chopchop_genomes/linearDB"
#buildlinearDB("hg38v34", path, Motif("Cas9"), 7, db_path)

guides = ["TGGGGGAGTTAGAGTTCTCC",
          "AAAGAATGTCCAATTACTGC",
          "CGCCCAGGCGGCCGCCAGAC",
          "GTATATCCAACAAAAACTTT",
          "TCTGTATCTCCCCATTGGTG",
          "CCCTGGTGACCACAGCGTAG",
          "CCATTTCTGGCCTCAGAGAC",
          "TGGCCTATGCGGAAGTAACC",
          "CGTCTCTGCGCTGGCAGATG",
          "TTGTACTCACAGCATGGGGA"]

# guides on fwd that we should find for sure
# CCCCCGCCAGCTTTCTCAGAAGTC chr1:155 CGG
#  CCCCGCCAGCTTTCTCAGAAGTCC chr1:156 GGG
# TCTCAGAAGTCCGGGTGAAAACCG chr1:165 AGG - check whether it should not be 168!
guides = [dna"CGCCAGCTTTCTCAGAAGTC"]

#sdb = load("/home/ai/.julia/dev/CRISPRofftargetHunter/test/sample_db/linearDB/" * string(reverse(guides[1])[1:7]) * ".bin")
#prefix_aln = Base.map(g -> prefix_levenshtein(g, a[1:7], 4), guides)
#isfinal = Base.map(x -> x.isfinal, prefix_aln)
#suffix_levenshtein(guides[1], sdb.suffix[2], prefix_aln[1], 4)

#1284.190450 seconds (24.00 G allocations: 666.087 GiB, 11.00% gc time)
# CTGAAGA

#g = mer"CTGAAGACTCTTTCGACCGC"
#suf =      mer"CTCTTTCGACCGCCCCC"
#pa = CRISPRofftargetHunter.PrefixAlignment([3, 3, 5, 6, 4, 5, 6, 6, 11, 6, 6, 6, 6, 16, 6, 6, 19, 20, 5, 4], 7, 3, false)
#suffix_levenshtein(g, suf, pa, 4)

#res = searchlinearDB(db_path, 1, guides)


#CRISPRofftargetHunter.levenshtein2(dna"ACTG", dna"GACTG")
#g = dna"CAGAGCCTCTAAGGTGAGC"
#ref = dna"ACAGACAGAGCGCTTTTGGAAAATGCG"
#@assert CRISPRofftargetHunter.levenshtein2(g, ref, 8) == CRISPRofftargetHunter.levenshtein(g, ref, 8)

levenshtein(dna"GCTGAAA", dna"ACTG", 4)
levenshtein2(dna"GCTGAAA", dna"ACTG", 4)
pa_sa(dna"GCTGAAA", dna"ACTG", 4, 4)


# CACT CTCCCTTT 4 3
g = dna"CACT"
ref = dna"CTCCCTTT"
#  12 3 45 67  8
# "CCTACTTGAACGCGT"
# "--T-C--G--CG-GT"
# hm the difference here is from swallowing the first letters!
align(g, ref, 4)
levenshtein(g, ref, 4)
pa_sa(g, ref, 4, 3)

# fast version with one vector instead of matrix
#median(@benchmark levenshtein(g2, r2, 8) setup=(g2=g; r2=ref))
#median(@benchmark levenshtein2_simple(g2, r2, 8) setup=(g2=g; r2=ref))

prefix_align(dna"AAAAAAA", dna"GGGGGGG", 16, 4)