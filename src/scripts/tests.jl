using BioSequences
using CRISPRofftargetHunter
using CSV

### TEST
genome = "./test/sample_data/sample_genome.fa"
motif = Motif("Cas9")
max_dist = 3

#genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
db = findofftargets(genome, motif, 1)
fillrate(db)

# my example guides that are inside
# ATGGAATTTTTGATTAAGTT AGG # chr 1
# TATATAGCTATTGAACATAA AGG # chr 2
# ACTGACTGACTGACTGACTG GGG # random string
guides = [dna"ATGGAATTTTTGATTAAGTT",
          dna"GGTAAAATTGGATTATTCGG",
          dna"ACTGACTGACTGACTGACTG"]

dt = estimate(db, guides)


#CSV.write("FileName.csv",  dt, writeheader=true)

db_path = "experiments_old/cas9_example_db.bin"
saveDB(db, db_path)
b = loadDB(db_path)
dt2 = estimate(b, guides)
dt2 == dt

# make 2bit work
# write some tests
# try building database on the human! - lets see how big is it

# 1. now lets make a file of all guides from our example genome as
output = Vector{String}()
findofftargets!(genome, motif, max_dist, output)
uguides = unique(output)
#open("test/sample_data/guides.txt", "w") do io
#      write(io, join(uguides, "\n"))
#end
# lets rank them by our method
guides = [LongDNASeq(g) for g in uguides]
dt = estimate(db, guides)
CSV.write("test/sample_data/guides_not_norm.csv", dt)

# 2. example input + some guides that are not there as negative input
neg_randoms = Vector{LongDNASeq}()
deletion_perm = deletion_permutations(20, max_dist, 1)
deletion_perm = [sort(vcat(x...)) for x in deletion_perm]
deletion_perm = [setdiff(collect(1:20), x) for x in deletion_perm]
uguides_deletes = Set(unique(vcat(output...)))
while length(neg_randoms) < 100
   rand = String(getSeq())
   # compute deletes for this
   deletes = Set([rand[x] for x in deletion_perm])
   if isdisjoint(deletes, uguides_deletes)
      push!(neg_randoms, LongDNASeq(rand))
   end
end

# check negative randoms!
dt_neg = estimate(db, neg_randoms)
CSV.write("test/sample_data/guides_neg.csv", dt_neg)
open("test/sample_data/guides_neg.txt", "w") do io
      write(io, join(neg_randoms, "\n")) # add NNN for crispritz
end

# 3. use CRISPRItz to build database and then query it for up to D4
