using BioSequences
using CSV
using DataFrames
using CRISPRofftargetHunter

genome_test = "./test/sample_data/samirandom.fa"

#f = FASTA.Reader(open(genome_test))
#Base.map(x -> FASTA.hasidentifier(x), f)
#close(f)
# run hg38v34 on all top guides from the top genes on 18.02/16:35

genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
motif = Motif("Cas9")
motif_ref = Motif("Cas9ref", "NNNNNNNNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXXXXXXXNGG")
max_dist = 6

# real rank table
guides = Vector{String}()
gatherofftargets!(genome_test, motif_ref, max_dist, guides)
guides = [LongDNASeq(g)[(max_dist + 1):end] for g in guides]
# 121 625 guides

# refg is faster!
# code another paralel version of findofftargets - paralel over guides - not chromosomes
p = CSV.read("test/sample_data/guides_linear_semirandom.csv", DataFrame)
brute = findofftargets_p_refg(genome, motif_ref, max_dist, guides[1:2])
any([occursin("-", g) for g in brute[!, "guide"]])
CSV.write("test/sample_data/guides_linear_hg38v34.csv", brute)

db = findofftargets(genome, motif, 3) # lets try this setting for hg38
db_path = "experiments_old/db_cas9_dist3_noExtension_unique_kmer.bin"
saveDB(db, db_path)
# lets save sketchDB


#p = CSV.read("test/sample_data/guides_linear_semirandom.csv", DataFrame)
#normal = [!occursin("-", x) for x in p[!, "guide"]]
#p = p[normal, :]

# we asume our p table is without wierd stuff...
iters = 100 # how many random iterations
k = 1000 # how many guides

subset_guides = p[rand(1:nrow(p), k), :]
sort!(subset_guides, [Symbol("D$i") for i in 0:max_dist])




# create example genome with 1k guides that are in complex relations!

# quick query based on pidgeon hole
a = "ACTGACTGACTGACTGACTG"
k = 4
s1 = Set()



function is_within_pidgeon_hole(guide::Set{String}, ref::Set{String})
    isdisjoint()
end

a = Set(each(DNAKmer{3}, dna"ATCCTANAGNTACT", 3))
