#=
using ARTEMIS
using BioSequences
using CSV
using DataFrames
using FASTX
using StatsBase
using StaticArrays
using BenchmarkTools

cd("test")

genome = joinpath(dirname(pathof(ARTEMIS)), "..", 
    "test", "sample_data", "genome", "semirandom.fa")
#genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNA{4}.(guides_s)

motif = Motif("Cas9")
tdir = tempname()
mkpath(tdir)

fmidir = tempname()
mkpath(fmidir)

motifpospath = build_motifDB("testCas9", genome, motif, tdir; store_kmers = true)
mdb_res = search_motifDB(tdir, guides, 3)
fmidbpath = build_fmiDB(genome, fmidir)
#res = search_fmiDB(fmidbpath, tdir, guides, 4)

=#
#=
dbi = ARTEMIS.DBInfo(genome, "tests", motif)

ref = open(dbi.gi.filepath, "r")
reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
chrom_name = dbi.gi.chrom[1]
record = reader[chrom_name]
chrom = dbi.gi.is_fa ? FASTA.sequence(record) : TwoBit.sequence(record)
close(ref)
# there is a lot of N in the genome - treat them as mismatches

composition(chrom)
# 16 r 32 -> 215mb
# 16 r 512 -> 187mb
# 16 r 1024 -> 186mb - why is it twice?!
index = ARTEMIS.FMIndex(chrom, 16; r = 32)
pam_pos = ARTEMIS.locateall(dna"AGG", index) # 3.7 Milion! - slowish
pos = ARTEMIS.locateall(dna"ACTG", index) # 1M - quite FAST
pos = ARTEMIS.locateall(dna"ACTGT", index) # 265k - very FAST
=#

# for distance of 4 we search for 20/5 = 4
# for distance of 3 we search for 20/4 = 5

# we overlap indexes with previously stored PAM indexes - in this fashion we can quickly
# reject what is wrong quickly
