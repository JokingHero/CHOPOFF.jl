__precompile__(true)

module CRISPRofftargetHunter

using Random
using BioSymbols
using BioSequences
using Serialization

include("utils.jl")
include("distance_metrics.jl")

export getSeq, file_read, file_write, file_add
export isinclusive, commonprefix, hamming, levenshtein
end
