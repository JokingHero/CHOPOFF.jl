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

# make and run default noHashDB
nhdb_path = joinpath(tdir, "noHashDB")
mkpath(nhdb_path)
build_noHashDB(
    "samirandom", genome, 
    Motif("Cas9", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 1, true, 0), 
    nhdb_path)

nhdb_path = "/tmp/jl_ce1oAx/noHashDB"
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)

nhdb_res = search_noHashDB(nhdb_path, guides)
                          
# should not be in there ACTCAATCATGTTTCCCGTC
# should be only once: AGAAGTCCGGGTGAAAACCG
#                     CAGAAGTCCGGGTGAAAACCG - what we store
=#