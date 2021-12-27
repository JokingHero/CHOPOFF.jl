
#=
using CRISPRofftargetHunter
using BioSequences
using CSV
using DataFrames
using FASTX
using StatsBase
using StaticArrays
using BenchmarkTools

cd("test")

genome = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..", 
    "test", "sample_data", "genome", "semirandom.fa")
#genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)

tdir = "/tmp/jl_AKyJgh"
fmidbpath = "/tmp/jl_Uj4tcC"

#guides = [dna"AGAGCGCCTGTGGTTGCCGG"] # "GGCCGTTGGTGTCCGCGAGACTCG" 42622 - on minus
# isplus 85215, 85216
# pos 108117 
# chrom semirandom3 0x03 0x07
#res = search_fmiDB(fmidbpath, tdir, guides, 4)
res = search_fmiDB_raw(fmidbpath, genome, Motif("Cas9"), guides)
mdb_res = search_motifDB(tdir, guides, 4)

=#