using BioSequences
using CSV
using DataFrames
using FASTX
using CRISPRofftargetHunter
using Random

# load guides
p = FASTA.Reader(open("/home/ai/Projects/uib/crispr/crispr_offtarget_search_benchmark/guides/risearch2_input.fa"))
guides = Base.map(x -> reverse_complement(FASTA.sequence(x))[4:end], p)
close(p)

guides = unique(rand(guides, 500))
all_guides = joinpath(
    "/home/ai/Projects/uib/crispr/",
    "CRISPRofftargetHunter/hg38v34_db.csv")

res = iterate_over_offtargets(guides, all_guides, "offtargets_3.csv")

# save results
res = DataFrame(res)
col_d = [Symbol("D$i") for i in 0:4]
rename!(res, col_d)
res.guide = guides
sort!(res, col_d)
CSV.write("guides_iterative_hg38v34_3.csv", res)
