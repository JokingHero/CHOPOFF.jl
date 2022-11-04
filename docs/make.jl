push!(LOAD_PATH, "../src/")
using Documenter, ARTEMIS
DocMeta.setdocmeta!(ARTEMIS, :DocTestSetup, :(using ARTEMIS, BioSequences); recursive=true)

makedocs(
    clean = true,
    strict = :doctest,
    doctest = true,
    sitename = "ARTEMIS.jl",
    authors = "Kornel Labun",
    format = Documenter.HTML(
        sidebar_sitename = false,
        footer = "Made with ♥ by Kornel Labun.",
        assets = ["assets/theme.css"],
        disable_git = true
    ),
    modules  = [ARTEMIS],
    pages = [
        "General" => "index.md",
        "API" => Any[
            "Abstract gRNA" => "abstract_gRNA.md",
            "Find potential off-targets" => "find_potential_ot.md",
            "Align gRNA and off-target" => "align_gRNA.md",
            "Alignment-free filters for gRNAs" => "alignment_free.md",
            "Find all off-targets" => "find_ot.md",
            "Utils" => "utils.md"],
        ])

deploydocs(
    repo = "github.com/JokingHero/ARTEMIS.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = true,
)
