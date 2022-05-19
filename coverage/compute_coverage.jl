using Coverage
using Badges
using Pkg
using CRISPRofftargetHunter

Pkg.test("CRISPRofftargetHunter"; coverage = true)

exclude = [
    "src/CRISPRofftargetHunter.jl", 
    "src/FMidx/saca.jl",
    "src/db_fmi.jl",
    "src/db_linear.jl",
    "src/sketches/bloom.jl",
    "src/motif.jl",
    "src/db_info.jl"]

function summary_by_file(x)
    p, c = get_summary(x)
    return p/c, x.filename
end

coverage = process_folder()
coverage = coverage[Base.map(p -> p.filename ∉ exclude, coverage)]
p, c = get_summary(merge_coverage_counts(coverage))


b = Badge(label="coverage", message = string(Int(ceil(p/c*100; digits = 0))) * "%")
coverage_badge =   # svg string
open("coverage/coverage_fraction.svg", "w") do io
    write(io, Badges.render(b))
end

# remove the coverage files
dir = joinpath(splitpath(dirname(pathof(CRISPRofftargetHunter)))[1:end-1])
cov_files = String[]
for (r, d, f) in walkdir(dir)
    cov = occursin.(".cov", f)
    append!(cov_files, joinpath.(r, f[cov]))
end
Base.map(rm, cov_files)