__precompile__(true)

module CRISPRofftargetHunter

using Random
using BioSymbols
using BioSequences

include("utils.jl")
include("distance_metrics.jl")

export isinclusive, commonprefix, hamming, levenshtein
end
