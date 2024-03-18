using Coverage
using Badges
using Pkg
using CHOPOFF

Pkg.test("CHOPOFF"; coverage = true)

exclude = [
    "src/CHOPOFF.jl", # contains hard to test command line interface
    "src/FMidx/saca.jl", # FMidx stuff is tested, but this package is not base CHOPOFF
    "src/FMidx/WaveletMatrices.jl",
    "src/FMidx/FMindexes.jl",
    "src/db_fmi_seed.jl", # methods using FM-index and are edge experimental
    "src/db_hash.jl", # experimental
    "src/example_doc.jl"] # docs

function summary_by_file(x)
    p, c = get_summary(x)
    return p/c, x.filename
end 

coverage = process_folder()
@info summary_by_file.(coverage)

coverage = coverage[Base.map(p -> p.filename ∉ exclude, coverage)]
p, c = get_summary(merge_coverage_counts(coverage))

b = Badge(label="coverage", message = string(Int(ceil(p/c*100; digits = 0))) * "%")
open("coverage/coverage_fraction.svg", "w") do io
    write(io, Badges.render(b))
end

# remove the coverage files
dir = joinpath(splitpath(dirname(pathof(CHOPOFF)))[1:end-1])
cov_files = String[]
for (r, d, f) in walkdir(dir)
    cov = occursin.(".cov", f)
    append!(cov_files, joinpath.(r, f[cov]))
end
Base.map(rm, cov_files)