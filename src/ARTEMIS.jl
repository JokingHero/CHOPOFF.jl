__precompile__(true)

module ARTEMIS

using Dates: isequal
using Base: Float32
using ArgParse: command_actions

using FastFilter
using StaticArrays

using CRC32c
using Dates
using Statistics
using StatsBase
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
using Pkg
using ArgParse
using CSV
using ProgressMeter
using VariantCallFormat
using CodecZlib
using PathDistribution
using Dates

include("example_doc.jl")
include("utils.jl")

include("FMidx/FMindexes.jl")
using .FMIndexes

include("ambig_index.jl")
include("persistence.jl")

include("distance_metrics.jl")
include("motif.jl")
include("db_info.jl")
include("motif_path_templates.jl")

include("find_offtargets.jl")
include("db_helpers.jl")
include("db_dict.jl")
include("db_motif.jl")
include("db_linear.jl")
include("db_tree.jl")
include("db_hash.jl")
include("db_vcf.jl")

include("db_fmi_helpers.jl")
include("db_fmi.jl")
include("db_fmi_seed.jl")


export Motif, length_noPAM, length, setambig, setdist # motif
export save, load # persistence
export minkmersize, getseq, expand_ambiguous, all_kmers, as_kmers, 
    as_skipkmers, summarize_offtargets, filter_overlapping # utils
# export AmbigIdx, findbits # ambig_index
export DBInfo # db_info
export gatherofftargets! # find_offtargets
export isinclusive, hamming, levenshtein, Aln, align # distance_metrics

export build_linearDB, search_linearDB # db_linear
export build_dictDB, search_dictDB # db_sketch
export build_treeDB, search_treeDB, inspect_treeDB # db_tree
export build_hashDB, search_hashDB # db_hash
export build_vcfDB, search_vcfDB # db_vcf
export build_motifDB, search_motifDB

export build_PathTemplates

export build_fmiDB, search_fmiDB
export build_pamDB, search_fmiDB_seed

## Standalone binary generation
function parse_commandline(args::Array{String})
    s = ArgParseSettings(
        prog = "ARTEMIS", 
        description = "Fast and reliable off-target detection for CRISPR guideRNAs.",
        epilog = "This standalone is build on top of the ARTEMIS.jl package.",
        version = Pkg.TOML.parse(read(joinpath(dirname(pathof(ARTEMIS)), "..", "Project.toml"), String))["version"],
        add_version = true)

    @add_arg_table! s begin
        "build", "B"
            help = "Build a database of the possible gRNA's/off-target sites."
            action = :command
        "search", "S"
            help = "Search database for off-targets."
            action = :command  
        "estimate", "E"
            help = "Estimate number of off-targets."
            action = :command   
        "filter", "F"
            help = "Filter overlapping off-targets, keep the one with the lowest distance."
            action = :command 
        "summarize", "U"
            help = "Summarize off-targets into table of counts by distance for each gRNA."
            action = :command   
    end

    @add_arg_table! s["build"] begin
        "treeDB"
            action = :command
            help = "treeDB uses Vantage Point tree and triangle inequality to find off-targets."
        "linearDB"
            action = :command
            help = "linearDB utilizes prefixes to decrease search time, but no other optimizations."
        "motifDB"
            action = :command
            help = "motifDB utilizes prefixes together with q-gram filtering."
        "hashDB"
            action = :command
            help = "hashDB is extremely fast, but only estimates off-targets within distance of 1"
        "dictDB"
            action = :command
            help = "dictDB is a simple dictionary of all unique guides and their counts."
        "vcfDB"
            action = :command
            help = "vcfDB is similar specialized database to handle .vcf files and personalized off-target search."
        "fmi"
            action = :command
            help = "Build FM-index of a genome"
        "pamDB"
            action = :command
            help = "Build DB of PAMs of a genome"
        "template"
            action = :command
            help = "Build templates with specific Motif."
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
            help = """Within what distance we can search for off-targets?"""
            arg_type = Int
            default = 3
        "--extend3" 
            help = """Defines how off-targets will be aligned to the guides and where
                extra nucleotides will be added for alignment within distance. Whether
                to extend in the 5' and 3' direction. Default is Cas9 with extend3 = false."""
            action = :store_true
        "--ambig_max"
            help = """How many ambiguous bases are allowed inside the guide?"""
            arg_type = Int
            default = 0
        "--motif"
            help = """Will try to get the Motif template based on the standard name e.g. Cas9 or Cpf1"""
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

    @add_arg_table! s["build"]["motifDB"] begin
        "--prefix_length"
            help = "Defines length of the prefix. " * 
                "For each possible prefix there will be " *
                "one motifDB instance."
            arg_type = Int
            default = 7
        "--skipmer_size"
            help = "Defines length of skipmer. "
            arg_type = Int
            required = false
    end


    @add_arg_table! s["build"]["hashDB"] begin
        "--max_count"
            help = "Maximum count value for given off-target sequence. " * 
                "Increasing this value will allow to keep more precise counts of often repeated guide sequences" * 
                ", but increase the size of the database significantly."
            arg_type = Int
            default = 10
        "--max_iterations"
            help = "How many iterations to try before quitting the building of the DB."
            arg_type = Int
            default = 10
        "--seed"
            help = "Initial seed for database build."
            arg_type = UInt64
            default = UInt64(0x726b2b9d438b9d4d)
        "--precision"
            help = "Whether to use UInt8, UInt16 or UInt32 to store the keys."
            arg_type = String
            range_tester = (
                x -> x == "UInt8" || 
                x == "UInt16" ||
                x == "UInt32")
            default = "UInt16"
    end

    @add_arg_table! s["build"]["dictDB"] begin
        "--max_count"
            help = "Maximum count value for given off-target sequence. " * 
                "Increasing this value will allow to keep more precise counts of often repeated guide sequences" * 
                ", but increase the size of the database significantly."
            arg_type = Int
            default = 10
    end

    @add_arg_table! s["build"]["vcfDB"] begin
        "--vcf"
            help = "Path to the .vcf or .vcf.gz file. "
            arg_type = String
            required = true
    end

    @add_arg_table! s["build"]["pamDB"] begin
        "--fmidir"
            help = "Path to the fmi directory. "
            arg_type = String
            required = true
    end

    @add_arg_table! s["search"] begin
        "treeDB"
            action = :command
            help = "treeDB uses Vantage Point tree and triangle inequality to find off-targets."
        "linearDB"
            action = :command
            help = "linearDB utilizes prefixes to decrease search time, but no other optimizations."
        "motifDB"
            action = :command
            help = "motifDB utilizes prefixes together with q-gram filtering."
        "vcfDB"
            action = :command
            help = "vcfDB is similar specialized database to handle .vcf files and personalized off-target search."
        "fmi"
            action = :command
            help = "Search fmi index using brute force method."
        "fmi_seed"
            action = :command
            help = "Search fmi index using lossless 01*0 seed method."
        "--distance"
            help = "Maximum edit distance to analyze. Must be less or equal to the distance that was used when building db."
            arg_type = Int
            default = 3
        "--database"
            help = "Path to THE FOLDER where the database is stored. Same as used when building. For FM-index based solutions this should be path to the FM-index folder."
            arg_type = String
            required = true
        "--guides"
            help = "File path to the guides, each row in the file contains a guide WITHOUT PAM."
            arg_type = String
            required = true
        "--output"
            help = "Path to the file where detailed results should be written."
            arg_type = String
            required = true
    end

    @add_arg_table! s["search"]["motifDB"] begin
        "--adjust"
            help = "Adjust parameter for motifDB."
            arg_type = Int
            default = 0
            required = false
    end

    @add_arg_table! s["search"]["fmi"] begin
        "--template"
            help = "Path to the table with the template. You can build a template with 'build  template'."
            arg_type = String
            required = true
    end

    @add_arg_table! s["search"]["fmi_seed"] begin
        "--genome"
            help = "Path to the genome."
            arg_type = String
            required = true
        "--pamDB"
            help = "Path to the file with pamDB. - Make with build_pamDB."
            arg_type = String
            required = true
    end

    @add_arg_table! s["estimate"] begin
        "--right"
            help = "Directionality in the filter databases."
            action = :store_true
        "--database"
            help = "Path to the file where the database is stored."
            arg_type = String
            required = true
        "--guides"
            help = "File path to the guides, each row in the file contains a guide WITHOUT PAM."
            arg_type = String
            required = true
        "--output"
            help = "File path to the file where summarized output should be generated."
            arg_type = String
            required = true
    end

    @add_arg_table! s["filter"] begin
        "--distance"
            help = "What is the distance for overlap filtering of off-targets. Reasonable selections are, distance used for "
            arg_type = Int
            required = true
        "--detail_file"
            help = "Path to the file where the database is stored."
            arg_type = String
            required = true
        "--output"
            help = "File path to the file where summarized output should be generated."
            arg_type = String
            required = true
    end


    @add_arg_table! s["summarize"] begin
        "--detail_file"
            help = "Path to the file where the database is stored."
            arg_type = String
            required = true
        "--output"
            help = "File path to the file where summarized output should be generated."
            arg_type = String
            required = true
    end

    return parse_args(args, s)
end


function main(args::Array{String})
    args = parse_commandline(args)
    
    for (arg, val) in args
        println(" $arg  =>  $val")
    end

    if args["%COMMAND%"] == "build"
        args = args["build"]

        if args["%COMMAND%"] != "fmi"
            if args["motif"] != ""
                motif = Motif(args["motif"])
                motif = setdist(motif, args["distance"])
                motif = setambig(motif, args["ambig_max"])
            else
                motif = Motif(
                    args["name"], args["fwd_motif"], 
                    args["fwd_pam"], !args["not_forward"], !args["not_reverse"],
                    args["distance"], !args["extend3"], args["ambig_max"])
            end
        end

        if args["%COMMAND%"] == "treeDB"
            build_treeDB(args["name"], args["genome"], motif, args["output"], 
                args["treeDB"]["prefix_length"])
        elseif args["%COMMAND%"] == "template"
            build_PathTemplates(length_noPAM(motif), motif.distance; 
                storagepath = joinpath(args["output"], args["name"] * ".bin"))
        elseif args["%COMMAND%"] == "linearDB"
            build_linearDB(args["name"], args["genome"], motif, args["output"], 
                args["linearDB"]["prefix_length"])
        elseif args["%COMMAND%"] == "motifDB"
            skipmer = args["motifDB"]["skipmer_size"]
            if skipmer === nothing
                skipmer = Int(floor(length_noPAM(motif) / (motif.distance + 3)))
            end
            build_motifDB(args["name"], args["genome"], motif, args["output"], 
                args["motifDB"]["prefix_length"]; skipmer_size = skipmer)
        elseif args["%COMMAND%"] == "hashDB"
            prec = UInt16
            if args["hashDB"]["precision"] == "UInt8"
                prec = UInt8
            elseif args["hashDB"]["precision"] == "UInt32"
                prec = UInt32
            end
            build_hashDB(args["name"], args["genome"], motif;
                storage_path = args["output"], 
                seed = args["hashDB"]["seed"], 
                max_iterations = args["hashDB"]["max_iterations"],
                max_count = args["hashDB"]["max_count"], precision = prec)
        elseif args["%COMMAND%"] == "dictDB"
            build_dictDB(args["name"], args["genome"], motif; storage_path = args["output"])
        elseif args["%COMMAND%"] == "vcfDB"
            build_vcfDB(args["name"], args["genome"], args["vcfDB"]["vcf"], motif; storage_path = args["output"])
        elseif args["%COMMAND%"] == "fmi"
            build_fmiDB(args["genome"], args["output"])
        elseif args["%COMMAND%"] == "pamDB"
            build_pamDB(args["pamDB"]["fmidir"], motif; storage_path = args["output"])
        else
            throw("Unsupported database type.")
        end
    elseif args["%COMMAND%"] == "search"
        args = args["search"]
        guides = LongDNA{4}.(readlines(args["guides"]))
        if args["%COMMAND%"] == "treeDB"
            res = search_treeDB(args["database"], guides, args["output"]; 
                distance = args["distance"])
        elseif args["%COMMAND%"] == "linearDB"
            res = search_linearDB(args["database"], guides, args["output"]; 
                distance = args["distance"])
        elseif args["%COMMAND%"] == "motifDB"
            search_motifDB(
                args["database"], guides, args["output"]; 
                distance = args["distance"], adjust = args["motifDB"]["adjust"])
        elseif args["%COMMAND%"] == "vcfDB"
            db = load(args["database"])
            res = search_vcfDB(db, guides)
        elseif args["%COMMAND%"] == "fmi"
            template = load(args["fmi"]["template"]) # TODO fix motif...
            search_fmiDB(guides, template, motif, args["database"], args["output"];
                distance = args["distance"])
        elseif args["%COMMAND%"] == "fmi_seed"
            pamDB = load(args["fmi_seed"]["pamDB"])
            search_fmiDB_seed(guides, args["database"], args["fmi_seed"]["genome"], pamDB, args["output"];
                distance = args["distance"])
        else
            throw("Unsupported database type.")
        end
        
    elseif args["%COMMAND%"] == "estimate"
        args = args["estimate"]
        guides = LongDNA{4}.(readlines(args["guides"]))
        db = load(args["database"])
        if isa(db, HashDB)
            @info "Right set as: " * string(args["right"])
            res = search_hashDB(db, guides, args["right"])
        elseif isa(db, DictDB)
            res = search_dictDB(db, guides)
        else
            throw("Unsupported database type.")
        end
        CSV.write(args["output"], res)

    elseif args["%COMMAND%"] == "filter"
        args = args["filter"]
        res = DataFrame(CSV.File(args["detail_file"]))
        res = filter_overlapping(res, args["distance"])
        CSV.write(args["output"], res)

    elseif args["%COMMAND%"] == "summarize"
        args = args["summarize"]
        res = DataFrame(CSV.File(args["detail_file"]))
        res = summarize_offtargets(res, maximum(res.distance))
        CSV.write(args["output"], res)
    end
end


function julia_main()::Cint
    try
        main(ARGS)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end


if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

end