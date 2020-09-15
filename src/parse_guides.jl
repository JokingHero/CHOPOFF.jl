
using CRISPRofftargetHunter
using BioSymbols
using BioSequences
using Statistics
# using JLD

# hamming, levenshtein

# 1. Calcualte average guide
# 2. Make a run through all guides and try to find off-targets - brute force

all_guides = joinpath(
    "/home/ai/Projects/uib/crispr/",
    "CRISPRofftargetHunter/hg38v34_db_test.csv",
)
k = 4 # max distance
max_dist = [5, 10, 25, 75, 250]

# avg_guide = fill(0, 5, 23)
# global row = 1
# for line in eachline(all_guides)
#     global row += 1
#     println(row)
#     line == "guide,location" && continue
#
#     guide = LongDNASeq(split(line, ",")[1])
#     #guide_dist = fill(0, k + 1)
#     for (i, ch) in enumerate(guide)
#         if ch == DNA_A
#             avg_guide[1, i] += 1
#         elseif ch == DNA_C
#             avg_guide[2, i] += 1
#         elseif ch == DNA_T
#             avg_guide[3, i] += 1
#         elseif ch == DNA_G
#             avg_guide[4, i] += 1
#         else
#             avg_guide[5, i] += 1
#         end
#     end
# end
# save("./avg_guide.jld", "avg_guide", avg_guide)
# saved_data = load("./avg_guide.jld")
# avg_guide = saved_data["avg_guide"]

struct Guide
    seq::LongDNASeq
    loci::Vector{String}
    #fiveprim::LongDNASeq
    #threprim::LongDNASeq
end

function Guide()
    return Guide(CRISPRofftargetHunter.getSeq(23 + 4), ["random"])
end

function Base.isless(x::Guide, y::Guide)
    return x.seq[4:end] > y.seq[4:end]
end

function Base.isequal(x::Guide, y::Guide)
    return isequal(x.seq[4:end], y.seq[4:end])
end

struct Bucket
    guides::Vector{Guide}
    d_to_p::Vector{Int}
end

function Bucket()
    return Bucket(Vector{Guide}(), Vector{Int}())
end

function Base.string(bucket::Bucket)
    return "g: " * string(length(bucket.guides))
end

function Base.push!(bucket::Bucket, guide::Guide, d::Int)
    push!(bucket.guides, guide)
    push!(bucket.d_to_p, d)
    return nothing
end

function Base.unique!(bucket::Bucket)
    ranks = sortperm(bucket.guides)
    # indices to remove
    to_remove = Vector{Int}()
    for i in 2:length(ranks)
        idx_g = findfirst(isequal(i), ranks)
        idx_g2 = findfirst(isequal(i-1), ranks)
        if isequal(bucket.guides[idx_g],
                   bucket.guides[idx_g2])
           # FIXME - how can this happen actually?
           # if !isequal(bucket.d_to_p[idx_g], bucket.d_to_p[idx_g2])
           #     print(bucket.d_to_p[idx_g])
           #     print(bucket.guides[idx_g])
           #     print(bucket.d_to_p[idx_g2])
           #     print(bucket.guides[idx_g2])
           #     error("Here...")
           # end
           push!(to_remove, idx_g2)
           append!(bucket.guides[idx_g].loci, bucket.guides[idx_g2].loci)
        end
    end
    deleteat!(bucket.guides, to_remove)
    deleteat!(bucket.d_to_p, to_remove)
    return nothing
end

struct Node
    guide::Guide
    radius::Int
    #furthest_d::Int
    #closest_d::Int
    inside::Int
    inside_bucket::Bool
    outside::Int
    outside_bucket::Bool
end

function Base.string(node::Node)
    return "r: " * string(node.radius)
end

function Node(guide::Guide)
    return Node(guide, round((length(guide)-7)/2), 0, false, 0, false)
end

function getindex(node::Node, inside::Bool = true)
    if inside
        return node.inside, node.inside_bucket
    else
        return node.outside, node.outside_bucket
    end
end

struct VPTree
    pam_len::Int
    pam_5prim::Bool
    guide_len::Int
    max_dist::Int
    nodes::Vector{Node}
    buckets::Vector{Bucket}
    max_bucket_len::Int
end

function VPTree()
    return VPTree(3, true, 20, 4, Vector{Node}(), Vector{Bucket}(), 500)
end

function updatenode!(
    tree::VPTree,
    node_idx::Int,
    idx::Int,
    inside::Bool,
    hasbucket::Bool = false
)
    if inside
        tree.nodes[node_idx] = Node(
            tree.nodes[node_idx].guide,
            tree.nodes[node_idx].radius,
            idx, hasbucket,
            tree.nodes[node_idx].outside,
            tree.nodes[node_idx].outside_bucket
        )
    else
        tree.nodes[node_idx] = Node(
            tree.nodes[node_idx].guide,
            tree.nodes[node_idx].radius,
            tree.nodes[node_idx].inside,
            tree.nodes[node_idx].inside_bucket,
            idx, hasbucket
        )
    end
    return nothing
end

"
    Will try to find value in `x` that will allow for almost equal
    split of values into buckets.
"
function balance(x::Vector{Int}, max_dist::Int = 4)
    uniq = unique(x)
    counts = [count(y -> y == i, x) for i in uniq]
    balance = argmin(abs.([sum(counts[1:i]) - sum(counts[i:end]) for i in 1:length(counts)]))
    adj_balance = uniq[balance] - max_dist
    return uniq[argmin(abs.(uniq .- adj_balance))]
end

function Base.push!(tree::VPTree, guide::Guide, node_idx::Int = 1)
    guide_inserted = false
    if length(tree.nodes) == 0
        push!(tree.buckets, Bucket())
        push!(tree.buckets, Bucket())
        push!(tree.nodes, Node(guide, round(tree.guide_len/2) - tree.max_dist, 1, true, 2, true))
        return nothing
    end

    while !guide_inserted

        # when guide == offtarget
        if hamming(guide.seq[4:end], tree.nodes[node_idx].guide.seq[4:end]) == 0
            append!(tree.nodes[node_idx].guide.loci, guide.loci)
            guide_inserted = true
        end

        d = levenshtein(guide.seq[4:23], tree.nodes[node_idx].guide.seq[4:end],
                        tree.max_dist + tree.nodes[node_idx].radius)
        isinside = (d - tree.max_dist) > tree.nodes[node_idx].radius
        child_idx, is_bucket = getindex(tree.nodes[node_idx], isinside)

        if is_bucket
            if length(tree.buckets[child_idx].guides) >= tree.max_bucket_len
                # compress guides to unique only
                unique!(tree.buckets[child_idx])
                # not median as we bias with restricted distance metric
                me = balance(tree.buckets[child_idx].d_to_p, tree.max_dist)
                # new split by guide close to median
                split_idx = rand(findall(tree.buckets[child_idx].d_to_p .== me))
                split_guide = tree.buckets[child_idx].guides[split_idx]
                deleteat!(tree.buckets[child_idx].guides, split_idx)
                deleteat!(tree.buckets[child_idx].d_to_p, split_idx)
                # iterate over all guides and compute new d
                new_d_to_p = fill(0, length(tree.buckets[child_idx].d_to_p))
                for (idx_g, g) in enumerate(tree.buckets[child_idx].guides)
                    new_d_to_p[idx_g] = levenshtein(g.seq[4:23], split_guide.seq[4:end], tree.max_dist + me)
                end
                # split into new buckets
                new_inside = (new_d_to_p .- tree.max_dist) .> me
                push!(tree.buckets, Bucket(tree.buckets[child_idx].guides[new_inside], new_d_to_p[new_inside]))
                deleteat!(tree.buckets[child_idx].guides, new_inside)
                deleteat!(tree.buckets[child_idx].d_to_p, new_inside)
                # add initial guide to the bucket
                d = levenshtein(guide.seq[4:23], split_guide.seq[4:end], tree.max_dist + me)
                if ((d - tree.max_dist) > me)
                    push!(tree.buckets[end], guide, d)
                else
                    push!(tree.buckets[child_idx], guide, d)
                end
                # make new node pointing to the buckets
                push!(tree.nodes, Node(split_guide, me, length(tree.buckets), true, child_idx, true))
                # update parent node to point to the above node
                updatenode!(tree, node_idx, length(tree.nodes), isinside)
            else
                push!(tree.buckets[child_idx], guide, d)
            end
            guide_inserted = true
        else
            node_idx = child_idx
        end
    end

    return nothing
end


function shiftdisp(f::String, o::String, xs::Vector{String})
    rep = repeat([o], length(xs) - 1)
    ch = vcat(f, rep)
    return map(*, ch, xs)
end


function drawSubTrees_(tree::VPTree, idx::Vector{Int}, isbucket::Vector{Bool}, isinside::Vector{Bool})
    if length(idx) != 0
        if length(idx) == 1
            return vcat(["│"], shiftdisp("└─ ", "   ",
                        drawTree(tree, idx[1], isbucket[1], isinside[1])))
        else
            return vcat(["│"], shiftdisp("├─ ", "│  ",
                        drawTree(tree, idx[1], isbucket[1], isinside[1])),
                        drawSubTrees_(tree, idx[2:end], isbucket[2:end], isinside[2:end]))
        end
    else
        return Vector{String}()
    end
end


function drawTree(tree::VPTree, idx::Int, isbucket::Bool, inside::Union{Nothing, Bool})
    node_pic = "☀ "
    if !isnothing(inside)
        if isbucket
            node_pic = inside ? "◧ " : "◨ " # left is inside
        else
            node_pic = inside ? "◐ " : "◑ "
        end
    end
    n_str = isbucket ? string(tree.buckets[idx]) : string(tree.nodes[idx])
    n = node_pic * string(idx) * " " * n_str
    if isbucket
        return [n]
    else
        node = tree.nodes[idx]
        subtree = drawSubTrees_(tree, [node.inside, node.outside],
                                [node.inside_bucket, node.outside_bucket],
                                [true, false])
        return vcat([n], subtree)
    end
end

function printVPtree(tree::VPTree, start_node::Int = 1, levels::Int = 10)
    if (start_node > length(tree.nodes))
        print("VP tree is empty.")
        return nothing
    end
    print(join(drawTree(tree, start_node, false, nothing), "\n"))
end

# function Base.show(io::IO, tree::VPTree)
#     printVPtree(tree::VPTree)
# end

vptree = VPTree(3, true, 20, 4, Vector{Node}(), Vector{Bucket}(), 3)
push!(vptree, Guide())
push!(vptree, Guide())
push!(vptree, Guide())
printVPtree(vptree)

global row = 1
global tree = VPTree()
for line in eachline(all_guides)
    global row += 1
    println(row)
    line == "guide,location" && continue

    line = split(line, ",")
    guide = Guide(LongDNASeq(line[1]), [line[2]])
    print(guide)
    push!(tree, guide)
end


# global row = 1
# for line in eachline(all_guides)
#     global row += 1
#     println(row)
#     line == "guide,location" && continue
#
#     guide = LongDNASeq(split(line, ",")[1])
#     #guide_dist = fill(0, k + 1)
#     skip = false
#     haskey(guides_map, guide) && continue
#
#     for line2 in eachline(all_guides)
#         line2 == "guide,location" && continue
#         offtarget = LongDNASeq(split(line2, ",")[1])
#         dist = levenshtein(guide, offtarget, k)
#         (dist == k + 1) && continue
#
#         guide_dist[dist+1] += 1
#
#         if any(guide_dist .> max_dist)
#             skip = true
#             break
#         end
#     end
#
#     if !skip && !haskey(guides_map, guide)
#         guide_dist[guide] = Guide(guide, guide_dist, split(line, ",")[2])
#     end
# end
