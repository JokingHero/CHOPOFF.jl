# include all files in test/src directory
files = readdir("src/")
for f in files
    include("src/" * f)
end