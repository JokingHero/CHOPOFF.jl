
using CRISPRofftargetHunter
using BioSymbols
using BioSequences
using Statistics
using Serialization
# using JLD

# hamming, levenshtein

# 1. Calcualte average guide
# 2. Make a run through all guides and try to find off-targets - brute force

all_guides = joinpath("/home/ai/Projects/uib/crispr/", "CRISPRofftargetHunter/hg38v34_db_test.csv")
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
    return Guide(getSeq(23 + 4), ["random"])
end

function Base.length(guide::Guide)
    return length(guide.loci)
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
    return "g:" * string(length(bucket.guides)) * " b:" * string(balance(bucket.d_to_p, 0))
end

function Base.length(bucket::Bucket)
    return length(bucket.guides)
end

function loci_count(bucket::Bucket)
    return sum(length.(bucket.guides))
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
    for i = 2:length(ranks)
        idx_g = findfirst(isequal(i), ranks)
        idx_g2 = findfirst(isequal(i - 1), ranks)
        if isequal(bucket.guides[idx_g], bucket.guides[idx_g2])
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
    sort!(to_remove)
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
    return "r:" * string(node.radius)
end

function loci_count(node::Node)
    return length(node.guide.loci)
end

function Node(guide::Guide)
    return Node(guide, round((length(guide) - 7) / 2), 0, false, 0, false)
end

function getindex(node::Node, inside::Bool = true)
    if inside
        return node.inside, node.inside_bucket
    else
        return node.outside, node.outside_bucket
    end
end

mutable struct VPtree
    pam_len::Int
    pam_5prim::Bool
    guide_len::Int
    max_dist::Int
    nodes::Vector{Node}
    bucket_size::Vector{Int}
    bucket_loci::Vector{Int}
    max_bucket_len::Int
    output_dir::String # persistent buckets are here
end

function VPtree(output_dir)
    return VPtree(3, true, 20, 4, Vector{Node}(), Vector{Int}(), Vector{Int}(), 500, output_dir)
end

function updatenode!(tree::VPtree, node_idx::Int, idx::Int, inside::Bool, hasbucket::Bool = false)
    if inside
        tree.nodes[node_idx] = Node(
            tree.nodes[node_idx].guide,
            tree.nodes[node_idx].radius,
            idx,
            hasbucket,
            tree.nodes[node_idx].outside,
            tree.nodes[node_idx].outside_bucket,
        )
    else
        tree.nodes[node_idx] = Node(
            tree.nodes[node_idx].guide,
            tree.nodes[node_idx].radius,
            tree.nodes[node_idx].inside,
            tree.nodes[node_idx].inside_bucket,
            idx,
            hasbucket,
        )
    end
    return nothing
end

"
    Will try to find value in `x` that will allow for almost equal
    split of values into buckets.
"
function balance(x::Vector{Int}, max_dist::Int = 4)
    if isempty(x)
        return nothing
    end
    uniq = unique(x)
    counts = [count(y -> y == i, x) for i in uniq]
    balance = argmin(abs.([sum(counts[1:i]) - sum(counts[i:end]) for i = 1:length(counts)]))
    adj_balance = uniq[balance] - max_dist
    return uniq[argmin(abs.(uniq .- adj_balance))]
end

function bucket_path(dir::String, idx::Int)
    gp = joinpath(dir, string("bucket_", idx, "_g.bin"))
    dp = joinpath(dir, string("bucket_", idx, "_d.bin"))
    return gp, dp
end

function push_to_bucket!(tree::VPtree, child_idx::Int, guide::Guide, d::Int)
    gp, dp = bucket_path(tree.output_dir, child_idx)
    file_add(gp, guide)
    file_add(dp, d)
    tree.bucket_size[child_idx] += 1
    tree.bucket_loci[child_idx] += length(guide)
    return nothing
end

" Delete previous bucket and push new content."
function push_to_bucket!(tree::VPtree, child_idx::Int, guide::Vector{Guide}, d::Vector{Int})
    gp, dp = bucket_path(tree.output_dir, child_idx)
    file_write(gp, guide)
    file_write(dp, d)
    tree.bucket_size[child_idx] = length(guide)
    tree.bucket_loci[child_idx] = sum(length.(guide))
    return nothing
end

function Base.push!(tree::VPtree, guide::Guide, node_idx::Int = 1)
    guide_inserted = false
    if length(tree.nodes) == 0
        push!(tree.bucket_size, 0)
        push!(tree.bucket_size, 0)
        push!(tree.bucket_loci, 0)
        push!(tree.bucket_loci, 0)
        push!(tree.nodes, Node(guide, round(tree.guide_len / 2) - tree.max_dist, 1, true, 2, true))
        return nothing
    end

    while !guide_inserted

        # when guide == offtarget
        if hamming(guide.seq[4:end], tree.nodes[node_idx].guide.seq[4:end]) == 0
            append!(tree.nodes[node_idx].guide.loci, guide.loci)
            break
        end

        d = levenshtein(
            guide.seq[4:23],
            tree.nodes[node_idx].guide.seq[4:end],
            tree.max_dist + tree.nodes[node_idx].radius,
        )
        isinside = (d - tree.max_dist) > tree.nodes[node_idx].radius
        child_idx, is_bucket = getindex(tree.nodes[node_idx], isinside)
        if is_bucket
            if tree.bucket_size[child_idx] >= tree.max_bucket_len
                # read bucket
                gp, dp = bucket_path(tree.output_dir, child_idx)
                gread = file_read(gp)
                dread = file_read(dp)
                bucket = Bucket(gread, dread)
                # compress guides to unique only
                unique!(bucket)
                # not median as we bias with restricted distance metric
                me = balance(bucket.d_to_p, tree.max_dist)
                # new split by guide close to median
                split_idx = rand(findall(bucket.d_to_p .== me))
                split_guide = bucket.guides[split_idx]
                deleteat!(bucket.guides, split_idx)
                deleteat!(bucket.d_to_p, split_idx)
                # iterate over all guides and compute new d
                new_d_to_p = fill(0, length(bucket.d_to_p))
                for (idx_g, g) in enumerate(bucket.guides)
                    new_d_to_p[idx_g] = levenshtein(g.seq[4:23], split_guide.seq[4:end], tree.max_dist + me)
                end
                # split into new buckets
                new_inside = (new_d_to_p .- tree.max_dist) .> me
                # make new bucket from part of the guides
                push!(tree.bucket_size, 0)
                push!(tree.bucket_loci, 0)
                push_to_bucket!(tree, length(tree.bucket_size),
                                bucket.guides[new_inside], new_d_to_p[new_inside])
                # overwrite old file with part of the guides
                push_to_bucket!(tree, child_idx, bucket.guides[.!new_inside], new_d_to_p[.!new_inside])
                # add initial guide to the bucket
                d = levenshtein(guide.seq[4:23], split_guide.seq[4:end], tree.max_dist + me)
                if ((d - tree.max_dist) > me)
                    push_to_bucket!(tree, length(tree.bucket_size), guide, d)
                else
                    push_to_bucket!(tree, child_idx, guide, d)
                end
                # make new node pointing to the buckets
                push!(tree.nodes, Node(split_guide, me, length(tree.bucket_size), true, child_idx, true))
                # update parent node to point to the above node
                updatenode!(tree, node_idx, length(tree.nodes), isinside)
            else
                push_to_bucket!(tree, child_idx, guide, d)
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

function drawVPSubTrees(tree::VPtree, idx::Vector{Int}, isbucket::Vector{Bool}, isinside::Vector{Bool}, levels::Int)
    if length(idx) != 0
        if length(idx) == 1
            return vcat(["│"], shiftdisp("└─ ", "   ", drawVPtree(tree, idx[1], isbucket[1], isinside[1], levels)))
        else
            return vcat(
                ["│"],
                shiftdisp("├─ ", "│  ", drawVPtree(tree, idx[1], isbucket[1], isinside[1], levels)),
                drawVPSubTrees(tree, idx[2:end], isbucket[2:end], isinside[2:end], levels),
            )
        end
    else
        return Vector{String}()
    end
end


function drawVPtree(tree::VPtree, idx::Int, isbucket::Bool, inside::Union{Nothing,Bool}, levels::Int)
    node_pic = "☀ "
    if !isnothing(inside)
        if isbucket
            node_pic = inside ? "◧ " : "◨ " # left is inside
        else
            node_pic = inside ? "◐ " : "◑ "
        end
    end
    n_str = isbucket ? string(tree.bucket_size[idx]) : string(tree.nodes[idx])
    n = node_pic * string(idx) * " " * n_str
    if isbucket || levels == 0
        return [n]
    else
        levels = levels - 1
        node = tree.nodes[idx]
        subtree = drawVPSubTrees(
            tree,
            [node.inside, node.outside],
            [node.inside_bucket, node.outside_bucket],
            [true, false],
            levels,
        )
        return vcat([n], subtree)
    end
end

function printVPtree(tree::VPtree, start_node::Int = 1, levels::Int = 3, io::IO = stdout)
    println(
        io,
        "Nodes:",
        length(tree.nodes),
        " Guides:",
        sum(tree.bucket_size),
        " Loci:",
        sum(tree.bucket_loci) + sum(loci_count.(tree.nodes)),
    )

    if start_node > length(tree.nodes) || get(io, :compact, false)
        return nothing
    end

    println(io, "\n" * join(drawVPtree(tree, start_node, false, nothing, levels), "\n"))
    return nothing
end

function Base.show(io::IO, mime::MIME"text/plain", tree::VPtree)
    printVPtree(tree, 1, 3, io)
end


function test_on_n(max_in_b::Int=500, n::Int=10000)
    temp = tempname()
    mkpath(temp)
    tree = VPtree(3, true, 20, 4, Vector{Node}(), Vector{Int}(), Vector{Int}(), max_in_b, temp)
    for i in 1:n
        push!(tree, Guide())
    end

    # confirm that tree has as many loci as specified on input
    for (i, s) in enumerate(tree.bucket_size)
        gp, dp = bucket_path(tree.output_dir, i)
        gp_ = file_read(gp)
        dp_ = file_read(dp)
        @assert length(dp_) == s
        @assert length(gp_) == s
        if (length(gp_)) == 0
            @assert tree.bucket_loci[i] == 0
        else
            @assert tree.bucket_loci[i] == sum(length.(gp_))
        end
    end
    @assert n == (sum(tree.bucket_loci) + sum(loci_count.(tree.nodes)))
    return tree
end

function read_in(all_guides::String, max_in_b::Int)
    temp = tempname()
    mkpath(temp)

    lci_count = 0
    tree = VPtree(3, true, 20, 4, Vector{Node}(), Vector{Int}(), Vector{Int}(), max_in_b, temp)
    for line in eachline(all_guides)
        line == "guide,location" && continue
        line = split(line, ",")
        guide = Guide(LongDNASeq(line[1]), [line[2]])
        push!(tree, guide)
        lci_count += 1

        if lci_count != (sum(tree.bucket_loci) + sum(loci_count.(tree.nodes)))
            throw(ErrorException(string("Loci count error ", lci_count)))
        end
    end

    if lci_count != (sum(tree.bucket_loci) + sum(loci_count.(tree.nodes)))
        println("Loci count error ")
    end
    return tree
end


#all_guides = joinpath("/home/ai/Projects/uib/crispr/", "CRISPRofftargetHunter/hg38v34_db.csv")
t = test_on_n()
t = read_in(all_guides, 3000)
