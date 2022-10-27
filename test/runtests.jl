# include all files in test/src directory
using ARTEMIS

dir = joinpath(pkgdir(ARTEMIS), "test", "src")
files = readdir(dir)
for f in files
    include(joinpath(dir, f))
end

#include(joinpath(dir, "db_fmi.jl"))
#include(joinpath(dir, "utils.jl"))
#include(joinpath(dir, "wavelet.jl"))