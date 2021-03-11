__precompile__(true)

module CRISPRofftargetHunter

using Statistics
using Random
using Combinatorics
using BioSymbols
using BioSequences
using Serialization
using FASTX
using TwoBit
using Probably
using DataFrames
using ThreadsX

import Base.findall

include("utils.jl")
include("distance_metrics.jl")
include("motif.jl")
include("bitoperations.jl")
include("persistence.jl")
include("find_offtargets.jl")
include("find_offtargets_p.jl")

export getSeq, file_read, file_write, file_add, bucket_path, deleterange, getkmers, minkmersize, getkgrams # utils
export isinclusive, commonprefix, hamming, levenshtein, levenshtein_bp # distance_metrics
export Motif # motif
export SketchDB, saveDB, loadDB # persistence
export deletion_permutations # bitoperations
export gatherofftargets, gatherofftargets!, estimate, fillrate, iterate_over_offtargets # find_offtargets
export findofftargets_p_chrom, findofftargets_p_refg # find_offtargets_p

end
