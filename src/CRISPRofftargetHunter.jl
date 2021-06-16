__precompile__(true)

module CRISPRofftargetHunter

using Base: Float32
using ArgParse: command_actions

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
using PkgVersion
using ArgParse

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

## Standalone binary generation
function parse_commandline()
    s = ArgParseSettings(
        prog = "CRISPRofftargetHunter", 
        description = "Fast and reliable off-target detection for CRISPR guideRNAs.",
        epilog = "This standalone is build on top of the CRISPRofftargetHunter.jl package.",
        version = string(@PkgVersion.Version),
        add_version = true)

    # We can't run tests for selfcontained app yet
    @add_arg_table! s begin
        "build", "B"
            help = "Build a database of the guideRNAs."
            action = :command
        "find", "F"
            help = "Find in a database all off-targets."
            action = :command       
    end  

    @add_arg_table! s["build"] begin
        "treeDB"
            action = :command
            help = "treeDB builds Vantage Point tree for efficient search of off-targets."
        "linearDB"
            action = :command
            help = "linearDB utilizes prefixes to decrease search time, but no other optimizations."
        "sketchDB"
            action = :command
            help = "sketchDB is extremally fast" * 
                ", but it can only estimate number of off-targets" * 
                ", however it can never under-estimate, it can only over-estimate."
        "--name"
            help = "How will you shortly name this database?"
            arg_type = String
            required = true
        "--genome"
            help = "Path to the genome, either .fa or .2bit"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Output folder where database will be stored"
            arg_type = String
            required = true
        "--fwd_motif"
            help = """
                Motif that indicates where is PAM inside `fwdpam`.
                For example for Cas9 it is 20*N + XXX:
                NNNNNNNNNNNNNNNNNNNNXXX"""
        "--fwd_pam"
            help = """
                Motif in 5'-3' that will be matched on the reference (without the X).
                For example for Cas9 it is 20*N + NGG:
                XXXXXXXXXXXXXXXXXXXXNGG"""
        "--not_forward"
            help = "If used will not match to the forward reference strand."
            action = :store_true
        "--not_reverse"
            help = "If used will not match to the reverse reference strand."
            action = :store_true
        "--distance"
            help = """How many extra nucleotides are needed for a search? This
            will indicate within what distance we can search for off-targets."""
            arg_type = Int
            default = 4
        "--extend3" 
            help = """Defines how off-targets will be aligned to the guides and where
                extra nucleotides will be added for alignment within distance. Whether
                to extend in the 5' and 3' direction. Cas9 is extend3 = false, default, without this option."""
            action = :store_true
        "--motif"
            help = """Will try to get the Motif from Motif databases avaialble are e.g. Cas9 or Cpf1"""
            arg_type = String
    end

    @add_arg_table! s["build"]["treeDB"] begin
        "--prefix_length"
            help = "Defines length of the prefix. " * 
                "For each possible prefix there will be " *
                "one vantage point tree."
            arg_type = Int
            default = 7
    end

    @add_arg_table! s["build"]["linearDB"] begin
        "--prefix_length"
            help = "Defines length of the prefix. " * 
                "For each possible prefix there will be " *
                "one linearDB instance."
            arg_type = Int
            default = 7
    end

    @add_arg_table! s["build"]["sketchDB"] begin
        "--probability_of_error"
            help = "Defines chance for over-estimating off-target counts. " * 
                "Decreasing this value will increase the size of the sketchDB."
            arg_type = Float64
            default = 0.001
        "--max_count"
            help = "Maximum count value for given off-target sequence. " * 
                "Increasing this value will allow to keep more precise counts of often repeated guide sequences" * 
                ", but increase the size of the database significantly."
            arg_type = Int
            default = 255
    end

    @add_arg_table! s["find"] begin
        "--distance"
            help = "Maximum edit distance to analyze. Must be less or equal to distance used when building db."
            arg_type = Int
            default = 3
        "--detail"
            help = "Path to the file where additional detailed results should be written. Not aplicable for sketchDB."
            arg_type = String
        "database"
            help = "Path to the folder where the database is stored. Same as used when building."
            arg_type = String
        "type"
            help = "Type of the database: treeDB, sketchDB or linearDB."
            arg_type = String
            range_tester = (x -> x == "sketchDB" || x == "treeDB" || x == "linearDB")
        "guides"
            help = "File path to the guides, each row in the file contains a guide WITHOUT PAM."
            arg_type = String
    end

    return parse_args(s)
end


function main()
    args = parse_commandline()
    
    for (arg, val) in args
        println(" $arg  =>  $val")
    end
end


function julia_main()::Cint
    try
        main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end

# julia --project=. --startup-file=no --trace-compile=app_precompile.jl ./src/CRISPRofftargetHunter.jl test
# julia --project=.
# using PackageCompiler
# create_app(".", "./compiled_app"; precompile_statements_file = "./app_precompile.jl", force = true)