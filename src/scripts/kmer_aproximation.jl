using BioSequences
using CSV
using DataFrames
using CRISPRofftargetHunter

genome_test = "./test/sample_data/samirandom.fa"
genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
motif = Motif("Cas9")
motif_ref = Motif("Cas9ref", "NNNNNNNNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXXXXXXXNGG")
max_dist = 6

# real rank table
guides = Vector{String}()
findofftargets!(genome_test, motif_ref, max_dist, guides)
guides = [LongDNASeq(g)[(max_dist + 1):end] for g in guides]

#ref = open(genome, "r")
#reader = FASTA.Reader(ref)
#reader = collect(enumerate(reader))
#close(ref)

p = CSV.read("test/sample_data/guides_linear_semirandom.csv", DataFrame)
wierd = Vector{LongDNASeq}() # 1422 wtf
for g in guides
    if !(string(g) in p[!, "guide"])
        push!(wierd, g)
    end
end

# 11:06
brute = findofftargets(genome_test, motif_ref, max_dist, wierd)
#CSV.write("test/sample_data/guides_linear_semirandom.csv", brute)




# create example genome with 1k guides that are in complex relations!
db = findofftargets(genome, motif_ref, )

n = 1000 # number of guides for randomization
max_dist = 4
total_seq = dna""

guides = Vector()
for i in 1:n
    println(i)
    rand_g = getSeq(20 + max_dist) # random guide
    offt = rand(0:20) # random number of offtargets
    offtargets = 0
    while offtargets < offt
        rand_off = getSeq(20 + max_dist)
        dist = levenshtein(rand_g[(max_dist + 1):end], rand_off, max_dist)
        if dist <= max_dist
            offtargets += 1
            push!(guides, rand_off)
        end
    end
    push!(guides, rand_g)
end

using Random
shuffle!(guides)

function wrapdna(guide, pfx_sfx = 20:30)
    return getSeq(rand(pfx_sfx)) * guide * getSeq(1) * dna"GG" * getSeq(rand(pfx_sfx))
end

guides = map(wrapdna, guides)
guides = reduce(*, guides)

middle = Int(ceil(length(guides) / 2))
chrom1 = repeat("N", 100) * string(guides[1:middle]) * repeat("N", 100)
chrom2 = repeat("N", 100) * string(guides[middle:end]) * repeat("N", 100)

using FASTX
f = open(FASTA.Writer, "test/sample_data/samirandom.fa")
write(f, FASTA.Record("semirandom1", LongDNASeq(chrom1)))
write(f, FASTA.Record("semirandom2", LongDNASeq(chrom2)))
close(f)

# test kmer strategy with

# run kmer strategy over night on hg38v34
