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

import Base.findall

include("utils.jl")
include("distance_metrics.jl")
include("motif.jl")
include("db_info.jl")
include("bitoperations.jl")

include("find_offtargets.jl")
include("find_offtargets_p.jl")

include("db_linear.jl")
include("db_sketch.jl")
#include("db_tree.jl")

include("persistence.jl")

# TODO only export small interaface for building databases and interactions with it
export getSeq, file_read, file_write, file_add, bucket_path, deleterange, comb_of_d # utils
export isinclusive, commonprefix, hamming, levenshtein, align, Aln, PrefixAlignment, pa_sa
export levenshtein_bp, suffix_align, prefix_align # distance_metrics
export Motif # motif
export DBInfo, Loc, decode # genomeinfo
export SketchDB, save, load # persistence
export deletion_permutations # bitoperations
export gatherofftargets, gatherofftargets! # find_offtargets
export findofftargets_p_chrom, findofftargets_p_refg # find_offtargets_p

# only this should be exported for the users - other functions are overly complicated
# and not usefull for them anyway but will bloat tutorial etc.
export buildlinearDB, searchlinearDB # db_linear
export buildsketchDB, searchsketchDB # db_sketch

end
