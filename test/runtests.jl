# include all files in test/src directory
using CRISPRofftargetHunter
dir = joinpath(pkgdir(CRISPRofftargetHunter), "test", "src")
files = readdir(dir)
for f in files
    include(joinpath(dir, f))
end