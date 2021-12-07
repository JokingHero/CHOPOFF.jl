#=
using CRISPRofftargetHunter
using BioSequences
using CSV
using DataFrames
using FASTX
using StatsBase
using StaticArrays

cd("test")
genome = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..", 
    "test", "sample_data", "genome", "semirandom.fa")
genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)
nhdb="/home/ai/tests_hg38v34/db/noHashDB"
#search_noHashDB(nhdb, guides)

sdb = CRISPRofftargetHunter.load(joinpath(nhdb, "noHashDB.bin"))
cidx = CRISPRofftargetHunter.ColumnIdx(sdb.guides, sdb.counts, 21)

CRISPRofftargetHunter.findbits(cidx, LongDNASeq(sdb.guides[1], 21))


=#

#=
genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
motif = Motif("Cas9")

if motif.distance != 1 || motif.ambig_max != 0
    motif = CRISPRofftargetHunter.setdist(motif, 1)
    motif = CRISPRofftargetHunter.setambig(motif, 0)
end
dbi = CRISPRofftargetHunter.DBInfo(genome, "test", motif)

guides = Vector{UInt64}()
ambig = CRISPRofftargetHunter.gatherofftargets!(guides, dbi)
=#



#=
cd("test")

## SET WD when debugging

genome = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..", 
    "test", "sample_data", "genome", "semirandom.fa")
genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)
tdir = tempname()
mkpath(tdir)

nhdb_path = joinpath(tdir, "noHashDB")
mkpath(nhdb_path)
build_noHashDB(
    "samirandom", genome, 
    Motif("Cas9"; distance = 1, ambig_max = 0), 
    nhdb_path)

guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)

search_noHashDB(nhdb_path, guides)
=#

#=
bdb_path = joinpath(tdir, "binDB")
mkpath(bdb_path)
build_binDB(
    "samirandom", genome, 
    Motif("Cas9", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 0, true, 0), 
    bdb_path)
bdb_res = search_binDB(bdb_path, guides, 1)

nhdb_res = search_hashDB(nhdb_path, guides, true)

using BioSequences
using CRISPRofftargetHunter
s = reverse(dna"ACCTAATTTTGGGGGGTCGG")

ext = CRISPRofftargetHunter.comb_of_d1_extended(string(s)) # all possible alignments with d1
#ext = collect(ext)
#ext = CRISPRofftargetHunter.comb_of_d(string(s), 1)
#ext = Set(vcat(ext[1], ext[2]))

extd1 = CRISPRofftargetHunter.comb_of_d1_extended_ref(string(s))
extd1_20 = Set(String.(SubString.(collect(extd1), 1, 20)))

i = setdiff(ext, extd1)
i = setdiff(i, extd1_20)
filter(x -> length(x) == 21, collect(i))

# bulge on the reference
# GAATGCGCCTATGGATGCGG - guide g20
# GAATGCGCCTATGGATGCGTG - ref 

s = dna"GAATGCGCCTATGGATGCGG"
ext, bord1 = CRISPRofftargetHunter.comb_of_d1(s)
norm, bord = CRISPRofftargetHunter.comb_of_d(string(s), 1)
norm = Set(vcat(LongDNASeq.(norm), LongDNASeq.(bord)))
setdiff(ext, norm)

norm = setdiff(norm, ext)
length.(norm) # all length 20!!!
norm
=#