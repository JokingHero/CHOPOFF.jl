__precompile__(true)

module CRISPRofftargetHunter

using CRC32c
using Dates
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
using Transducers
using ThreadsX

include("utils.jl")
include("persistence.jl")

include("distance_metrics.jl")
include("motif.jl")
include("db_info.jl")

include("find_offtargets.jl")
include("db_helpers.jl")
include("db_sketch.jl")
include("db_linear.jl")
include("db_tree.jl")

export Motif # motif
export gatherofftargets! # find_offtargets
export build_linearDB, search_linearDB # db_linear
export build_sketchDB, search_sketchDB # db_sketch
export build_treeDB, search_treeDB, inspect_treeDB # db_tree

end
