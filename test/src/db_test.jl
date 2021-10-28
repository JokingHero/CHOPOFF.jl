#=
using CRISPRofftargetHunter
using BioSequences
using CSV
using DataFrames

cd("test")

## SET WD when debugging

genome = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..", 
    "test", "sample_data", "genome", "semirandom.fa")
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)
tdir = tempname()
mkpath(tdir)

nhdb_path = joinpath(tdir, "hashDB")
mkpath(nhdb_path)
build_hashDB(
    "samirandom", genome, 
    Motif("Cas9", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 1, true, 0), 
    nhdb_path)

nhdb_path = "/tmp/jl_r4Wgh4/hashDB"
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)

#=
bdb_path = joinpath(tdir, "binDB")
mkpath(bdb_path)
build_binDB(
    "samirandom", genome, 
    Motif("Cas9", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 0, true, 0), 
    bdb_path)
bdb_res = search_binDB(bdb_path, guides, 1)
=#
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