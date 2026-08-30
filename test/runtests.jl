# include all files in test/src directory
using CHOPOFF

dir = joinpath(pkgdir(CHOPOFF), "test", "src")
# Filter on the extension: an editor backup or a stray fixture dropped in
# test/src used to be passed to include() and crash the whole suite.
files = sort(filter(f -> endswith(f, ".jl"), readdir(dir)))
for f in files
    include(joinpath(dir, f))
end
