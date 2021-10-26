__precompile__(true)

module CRISPRofftargetHunter

using Dates: isequal
using Base: Float32
using ArgParse: command_actions
using BioSymbols: isambiguous

using FastFilter
using StaticArrays

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
using DataFrames
using Transducers
using ThreadsX
using PkgVersion
using ArgParse
using CSV
using ProgressMeter

include("sketches/hyperloglog.jl")
include("sketches/cms.jl")
include("sketches/bloom.jl")
include("ambig_index.jl")

include("utils.jl")
include("persistence.jl")

include("distance_metrics.jl")
include("motif.jl")
include("db_info.jl")

include("find_offtargets.jl")
include("db_helpers.jl")
include("db_dict.jl")
include("db_sketch.jl")
include("db_linear.jl")
include("db_compressed.jl")
include("db_tree.jl")
include("db_bins.jl")
include("db_hash.jl")
include("db_large_nohash.jl")

export Motif # motif
export build_linearDB, search_linearDB # db_linear
export build_compactDB, search_compactDB # db_compressed
export build_dictDB, search_dictDB # db_sketch
export build_treeDB, search_treeDB, inspect_treeDB # db_tree
export build_binDB, search_binDB # db_bins
export build_hashDB, search_hashDB # db_hash
export build_noHashDB, search_noHashDB # db_large_nohash


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
        "search", "S"
            help = "Search database for off-targets."
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
        "dictDB"
            action = :command
            help = "dictDB is a simple dictionary of all unique guides and their counts."
        "binDB"
            action = :command
            help = "binDB is a binned bloom filter, it is very space efficient " * 
                "and carries small error during estimations."
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
            default = ""
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
        "--error_size"
            help = "Transforms into length of the table, it is epsilon indicating what is desired " *
                "error of the estimated count we can afford."
            arg_type = Int
            default = 3
        "--max_count"
            help = "Maximum count value for given off-target sequence. " * 
                "Increasing this value will allow to keep more precise counts of often repeated guide sequences" * 
                ", but increase the size of the database significantly."
            arg_type = Int
            default = 255
    end

    @add_arg_table! s["build"]["dictDB"] begin
        "--max_count"
            help = "Maximum count value for given off-target sequence. " * 
                "Increasing this value will allow to keep more precise counts of often repeated guide sequences" * 
                ", but increase the size of the database significantly."
            arg_type = Int
            default = 255
    end

    @add_arg_table! s["build"]["binDB"] begin
        "--probability_of_error"
            help = "Defines chance for error off-target counts. " * 
                "Decreasing this value will increase the size of the sketchDB."
            arg_type = Float64
            default = 0.001
        "--max_count"
            help = "Maximum count value for given off-target sequence. " * 
                "Increasing this value will allow to keep more precise counts of often repeated guide sequences" * 
                ", but increase the size of the database significantly."
            arg_type = Int
            default = 10
    end

    @add_arg_table! s["search"] begin
        "--distance"
            help = "Maximum edit distance to analyze. Must be less or equal to distance used when building db."
            arg_type = Int
            default = 1
        "--detail"
            help = "Path to the file where additional detailed results should be written. Not aplicable for sketchDB."
            arg_type = String
            default = ""
        "database"
            help = "Path to the folder where the database is stored. Same as used when building."
            arg_type = String
            required = true
        "type"
            help = "Type of the database: treeDB, sketchDB or linearDB."
            arg_type = String
            range_tester = (x -> x == "sketchDB" || x == "binDB" || x == "dictDB" || x == "treeDB" || x == "linearDB")
            required = true
        "guides"
            help = "File path to the guides, each row in the file contains a guide WITHOUT PAM."
            arg_type = String
            required = true
        "output"
            help = "File path to the file where summarized output should be generated."
            arg_type = String
            required = true
    end

    return parse_args(s)
end


function main()
    args = parse_commandline()
    
    for (arg, val) in args
        println(" $arg  =>  $val")
    end

    if args["%COMMAND%"] == "build"
        args = args["build"]
        if args["motif"] != ""
            motif = Motif(args["motif"])
        else
            motif = Motif(
                args["name"], args["fwd_motif"], 
                args["fwd_pam"], !args["not_forward"], !args["not_reverse"],
                args["distance"], !args["extend3"])
        end

        if args["%COMMAND%"] == "treeDB"
            build_treeDB(args["name"], args["genome"], motif, args["output"], 
                args["treeDB"]["prefix_length"])
        elseif args["%COMMAND%"] == "linearDB"
            build_linearDB(args["name"], args["genome"], motif, args["output"], 
                args["linearDB"]["prefix_length"])
        elseif args["%COMMAND%"] == "sketchDB"
            #build_sketchDB(args["name"], args["genome"], motif, args["output"], 
            #    args["sketchDB"]["probability_of_error"], args["sketchDB"]["error_size"]; 
            #    max_count = args["sketchDB"]["max_count"])
        elseif args["%COMMAND%"] == "dictDB"
            build_dictDB(args["name"], args["genome"], motif, args["output"]; 
                max_count = args["sketchDB"]["max_count"])
        elseif args["%COMMAND%"] == "binDB"
            build_binDB(args["name"], args["genome"], motif, args["output"], 
                args["sketchDB"]["probability_of_error"]; 
                max_count = args["sketchDB"]["max_count"])
        else
            throw("Unsupported database type.")
        end
    elseif args["%COMMAND%"] == "search"
        args = args["search"]
        guides = LongDNASeq.(readlines(args["guides"]))
        if args["type"] == "treeDB"
            res = search_treeDB(args["database"], guides, args["distance"]; detail = args["detail"])
        elseif args["type"] == "linearDB"
            res = search_linearDB(args["database"], guides, args["distance"]; detail = args["detail"])
        elseif args["type"] == "sketchDB"
            #res = search_sketchDB(args["database"],  guides, args["distance"])
        elseif args["type"] == "dictDB"
            res = search_dictDB(args["database"],  guides, args["distance"])
        elseif args["type"] == "binDB"
            res = search_binDB(args["database"],  guides, args["distance"])
        else
            throw("Unsupported database type.")
        end
        CSV.write(args["output"], res)
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