function make_example_doc(method::String = "linearDB")
    return """
# make a temporary directory
tdir = tempname()
db_path = joinpath(tdir, "$method")
mkpath(db_path)

# use ARTEMIS example genome
artemis_path = splitpath(dirname(pathof(ARTEMIS)))[1:end-1]
genome = joinpath(
    vcat(
        artemis_path, 
        "test", "sample_data", "genome", "semirandom.fa"))

# build a $method
build_$method(
    "samirandom", genome, 
    Motif("Cas9"), 
    db_path)

# load up example gRNAs
using BioSequences
guides_s = Set(readlines(joinpath(vcat(artemis_path, "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)
    
# finally, make results!
res_path = joinpath(tdir, "$method", "results.csv")
search_$method(db_path, guides, res_path; distance = 3)

# load results
using DataFrames, CSV
res = DataFrame(CSV.File(res_path))

# filter results by close proximity
res = filter_overlapping(res, 23)

# summarize results into a table of counts by distance
summary = summarize_offtargets(res, 3)
"""
end