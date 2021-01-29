__precompile__(true)

module CRISPRofftargetHunter

using Random
using Combinatorics
using BioSymbols
using BioSequences
using Serialization
using FASTX
using TwoBit
using Probably
using DataFrames

include("utils.jl")
include("distance_metrics.jl")
include("motif.jl")
include("bitoperations.jl")
include("persistence.jl")
include("find_offtargets.jl")

export getSeq, file_read, file_write, file_add, bucket_path, deleterange
export isinclusive, commonprefix, hamming, levenshtein # distance_metrics
export Motif # motif
export SketchDB, saveDB, loadDB # persistence
export deletion_permutations # bitoperations
export findofftargets, findofftargets!, estimate, fillrate # find_offtargets

end
