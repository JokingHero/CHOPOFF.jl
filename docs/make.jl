push!(LOAD_PATH, "../src/")
using Documenter, ARTEMIS
DocMeta.setdocmeta!(ARTEMIS, :DocTestSetup, :(using ARTEMIS, BioSequences); recursive=true)

makedocs(
    clean = true,
    strict = :doctest,
    doctest = true,
    sitename = "ARTEMIS.jl",
    authors = "Kornel Labun",
    modules  = [ARTEMIS],
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
