using Documenter, CRISPRofftargetHunter

makedocs(
    sitename = "CRISPRofftargetHunter.jl",
    authors = "Kornel Labun",
    pages = ["Home" => "index.md"])

#=
deploydocs(
    repo = "github.com/JuliaDocs/Documenter.jl.git",
    target = "build",
    push_preview = true,
)
=#