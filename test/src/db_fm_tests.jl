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

motif = Motif("Cas9")
motif = CRISPRofftargetHunter.setdist(motif, 2)

CRISPRofftargetHunter.as_partial_alignments("AAAAAAAAAAAAAAAAAAAA", motif, 10)
#fmidir = tempname()
#mkpath(fmidir)
#fmidbpath = build_fmiDB(genome, fmidir)

tdir = "/tmp/jl_AKyJgh"
fmidbpath = "/tmp/jl_1PWg8q"

#guides = [dna"AGAGCGCCTGTGGTTGCCGG"] # "GGCCGTTGGTGTCCGCGAGACTCG" 42622 - on minus
# isplus 85215, 85216
# pos 108117 
# chrom semirandom3 0x03 0x07
#res = search_fmiDB(fmidbpath, tdir, guides, 4)
res = search_fmiDB_raw(fmidbpath, genome, motif, guides)
#mdb_res = search_motifDB(tdir, guides, 4)
=#