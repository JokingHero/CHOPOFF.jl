function make_example_doc(method::String = "linearDB"; search::String = method)
    return """
using ARTEMIS, BioSequences

# make a temporary directory
tdir = tempname()
db_path = joinpath(tdir, "$method")
mkpath(db_path)

# use ARTEMIS example genome
artemis_path = splitpath(dirname(pathof(ARTEMIS)))[1:end-1]
genome = joinpath(vcat(artemis_path, 
    "test", "sample_data", "genome", "semirandom.fa"))

# build a $method
build_$method(
    "samirandom", genome, 
    Motif("Cas9"), 
    db_path)

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(artemis_path, 
    "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)
    
# finally, make results!
res_path = joinpath(tdir, "$search", "results.csv")
search_$search(db_path, guides, res_path; distance = 3)

# load results
using DataFrames, CSV
res = DataFrame(CSV.File(res_path))

# filter results by close proximity
res = filter_overlapping(res, 23)

# summarize results into a table of counts by distance
summary = summarize_offtargets(res; distance = 3)
"""
end

function make_vcf_example_doc()
    return """
# prepare libs
using ARTEMIS, BioSequences

# use ARTEMIS example genome and vcf file
artemis_path = splitpath(dirname(pathof(ARTEMIS)))[1:end-1]
genome = joinpath(vcat(artemis_path, 
    "test", "sample_data", "genome", "semirandom.fa"))
vcf = joinpath(vcat(artemis_path, 
    "test", "sample_data", "artificial.vcf"))

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(artemis_path, 
    "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)

# example VCF file
vcf_db = build_vcfDB(
    "samirandom", genome, vcf,
    Motif("Cas9"; distance = 1, ambig_max = 0))

# search using vcfDB!
vcf_res = search_vcfDB(vcf_db, guides)
"""
end