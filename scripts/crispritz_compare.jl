using CSV
using DataFrames
using BioSequences

cz_file = "./test/sample_data/crispritz_results/guides.output.targets.txt"
cz = DataFrame(CSV.File(cz_file))

guides = Set(readlines("./test/sample_data/crispritz_results/guides.txt"))

ldb_file = "/home/ai/semirandom_detailed.csv"
ldb = DataFrame(CSV.File(ldb_file))
# filter out ldb to only guides
ldb = ldb[[x in guides for x in ldb.guide], :]
max_d = maximum(ldb.distance)

cz = cz[cz.Total .<= 3, :]

# remove 3 N at the begining and - to make guide
function asguide(x::String)
    x = x[1:(length(x) - 3)]
    x = replace.(x, "-" => "")
    @assert length(x) == 20
    return x
end
cz.guide = asguide.(cz.crRNA)
#cz.alignment_reference = [uppercase(x[1:(length(x) - 3)]) for x in cz.DNA]
#cz.alignment_guide = [uppercase(x[1:(length(x) - 3)]) for x in cz.crRNA]

function countspaces(x)
    return length(findall("-", x))
end

function ldb_start(pos, rrna_len, czdna_spac, strand)
    start = Vector{Int}()
    for i in 1:length(pos)
        if strand[i] == "+"
            i_start = pos[i] + rrna_len[i] - czdna_spac[i]
        else
            i_start = pos[i] + 1
        end
        push!(start, i_start)
    end
    return start
end
cz.start = ldb_start(cz.Position, length.(cz.crRNA), countspaces.(cz.DNA), cz.Direction)

function rows_not_in(len::Int, rows_in::Vector{Int})
    all_rows = collect(1:len)
    rows_in_ = sort(unique(rows_in))
    if len == length(rows_in_)
        return []
    end
    deleteat!(all_rows, rows_in)
    return all_rows
end

# list all alignments that are in linearDB and not in crizpritz output
for g in guides
    #g = "AGAGCGCCTGTGGTTGCCGG"
    czg = cz[cz.guide .== g, :]
    ldbg = ldb[ldb.guide .== g, :]

    isfound_in_czg = Vector{Int}() # indexes of rows of ldbg that are found in czg
    isfound_in_ldbg = Vector{Int}() # indexes of rows of czg that are found in ldbg
    for j in 1:nrow(ldbg)
        issame = (czg.Total .== ldbg.distance[j]) .&
            (czg.Direction .== ldbg.strand[j]) .&
            (czg.start .== ldbg.start[j]) .&
            (czg.Chromosome .== ldbg.chromosome[j])
        if any(issame)
            j_idx = findall(issame)
            j_idx = j_idx[1]
            push!(isfound_in_czg, j)
            isfound_in_ldbg = vcat(isfound_in_ldbg, findall(
                (czg.start .== czg.start[j_idx]) .&
                (czg.Chromosome .== czg.Chromosome[j_idx]) .&
                (czg.Direction .== czg.Direction[j_idx])))
        end
    end

    ldbg = ldbg[rows_not_in(nrow(ldbg), isfound_in_czg), :]
    czg = czg[rows_not_in(nrow(czg), isfound_in_ldbg), :]

    # crispritz can't handle anything too close to the telomeres NNN
    ldbg = ldbg[ldbg.start .> 125, :]

    if nrow(ldbg) > 0
        @info "Failed ldbg finding $g"
        @info "$ldbg"
        throw("$g")
    end

    if nrow(czg) > 0
        @info "Failed cz finding $g"
        @info "$czg"
        throw("$g")
    end
    #ldbg[j, :]
    #czg[5, [:alignment_guide, :alignment_reference, :Direction, :start, :Chromosome]]

    #CSV.write("/home/ai/czg.csv", czg)
    #CSV.write("/home/ai/ldbg.csv", ldbg)
end
