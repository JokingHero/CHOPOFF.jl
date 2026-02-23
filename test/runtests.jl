# include all files in test/src directory
using CHOPOFF

dir = joinpath(pkgdir(CHOPOFF), "test", "src")
files = readdir(dir)
for f in files
    include(joinpath(dir, f))
end
