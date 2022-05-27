push!(LOAD_PATH, "../src/")
using Documenter, CRISPRofftargetHunter
DocMeta.setdocmeta!(CRISPRofftargetHunter, :DocTestSetup, :(using CRISPRofftargetHunter, BioSequences); recursive=true)

makedocs(
    clean = true,
    strict = :doctest,
    doctest = true,
    sitename = "CRISPRofftargetHunter.jl",
    authors = "Kornel Labun",
    modules  = [CRISPRofftargetHunter],
    pages = [
        "General" => "index.md",
        "API" => "api.md",
        ])

#=
deploydocs(
    repo = "github.com/JuliaDocs/Documenter.jl.git",
    target = "build",
    push_preview = true,
)
=#
