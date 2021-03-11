using BioSequences
#using CSV
#using DataFrames
#using FASTX
using CRISPRofftargetHunter
using StatsPlots

all_guides = joinpath(
    "/home/ai/Projects/uib/crispr/",
    "CRISPRofftargetHunter/hg38v34_db.csv")
test_guides = joinpath(
    "/home/ai/Projects/uib/crispr/",
    "CRISPRofftargetHunter/hg38v34_db_test.csv")

function prefix_count(ug::Set{String}, prefix_len = 4)
    kmers_str = IdDict{String, Int}()

    for g in ug
        ref = g[1:prefix_len + 1]
        ref_count = get(kmers_str, ref, 0)
        kmers_str[ref] = ref_count + 1
    end
    return kmers_str
end

function load_uguides(guides_table_path::String)
    ug = Set{String}()

    f = open(guides_table_path)
    i = 0
    w = 0
    for ln in eachline(f)
        if ln == "guide,location"
            continue
        end
        i += 1
        w += 1
        if w == 100000
            println("Iter: ", i)
            w = 0
        end
        ln = split(ln, ",")
        push!(ug, String(ln[1][4:end]))
    end
    close(f)
    return ug
end

test_ug = load_uguides(test_guides)
prefix_count(test_ug)
ug = load_uguides(all_guides)
prefixes = prefix_count(ug, 8)

counts = collect(values(prefixes))
maximum(counts) #821950
minimum(counts) #13337
length(prefixes) #1024
sum(map(x -> x == 1, counts))

histogram(counts)

# what happends when we increase k
# what % of trees we have to check for all given kmers
kmers = collect(keys(prefixes))
4^(4 + 1) # maximum unique kmers for given k

dist = zeros(Int, length(kmers), length(kmers))
for (i, kmer1) in enumerate(kmers)
    for (j, kmer2) in enumerate(kmers)
        print(kmer1 * " " * kmer2)
        dist[i, j] = levenshtein(LongDNASeq(kmer1), LongDNASeq(kmer2))
    end
end

perc = []
for i in 1:length(kmers)
    push!(perc, sum(map(x -> x > 4, dist[i, :])))
end
histogram(perc)
