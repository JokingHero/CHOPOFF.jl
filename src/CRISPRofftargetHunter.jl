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
using ThreadsX

import Base.findall

include("utils.jl")
include("distance_metrics.jl")
include("motif.jl")
include("genomeinfo.jl")
include("bitoperations.jl")
include("persistence.jl")

include("find_offtargets.jl")
include("find_offtargets_p.jl")
include("build_db_sketch.jl")
include("search_db_sketch.jl")
include("build_db.jl")
include("search_db.jl")


export getSeq, file_read, file_write, file_add, bucket_path, deleterange, getkmers, minkmersize, getkgrams # utils
export isinclusive, commonprefix, hamming, levenshtein, levenshtein_bp # distance_metrics
export Motif # motif
export GenomeInfo, Locus, decode # genomeinfo
export SketchDB, saveDB, loadDB # persistence
export deletion_permutations # bitoperations
export gatherofftargets, gatherofftargets!, estimate, fillrate, iterate_over_offtargets # find_offtargets
export findofftargets_p_chrom, findofftargets_p_refg # find_offtargets_p

end
