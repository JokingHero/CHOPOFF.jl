using BioSequences
using CRISPRofftargetHunter
using CSV

### TEST
genome = "./test/sample_data/sample_genome.fa"
motif = Motif("Cas9")
max_dist = 3

# 1. Try using biore"" eachmatch instead!

function findofftargets!2(
    genome::String,
    motif::Motif,
    max_dist::Int,
    output::Vector{T}) where T<:Union{CountMinSketch, HyperLogLog}

    ext = extension(genome)
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        throw(
        "Wrong extension of the genome.",
        "We can parse only .fa/.fna/.fasta or .2bit references.")
    end

    ref = open(genome, "r")
    reader = is_fa ? FASTA.Reader(ref) : TwoBit.Reader(ref)

    for record in reader
        chrom = is_fa ? FASTA.sequence(record) : TwoBit.sequence(record)

        if FASTA.hasidentifier(record) # TODO 2bit?!
            println("Working on: ", FASTA.identifier(record))
        end

        pushdeletes!(chrom, motif.fwd, motif.pam_loci_fwd, max_dist, false, output)
        pushdeletes!(chrom, motif.rve, motif.pam_loci_rve, max_dist, true, output)
    end

    close(ref)
    return
end

BioSequences.RE.Regex{DNA}("ACTG", :pcre)

#genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
db = findofftargets(genome, motif, max_dist)
# 1543


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

# now lets make a file of all guides from our example genome as
# example input + some guides that are not there as negative input

# use CRISPRItz to build database and then query it for up to D4
