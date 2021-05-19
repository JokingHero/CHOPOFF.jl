"
Each node contains one guide (suffix) and 
indexes of that guide loci (loci_idx).
This guide also contains a radius and pointers
to other Nodes that are inside and outside the 
radius.
"
struct Node
    suffix::LongSequence{DNAAlphabet{4}}
    loci_idx::LociRange
    radius::UInt8
    inside::UInt32
    outside::UInt32
end


"
Final SuffixVPtreeDB unit that contains all guides from
all chromosomes that start with the `prefix` and their locations.
"
struct SuffixTreeDB
    prefix::LongSequence{DNAAlphabet{4}}
    nodes::Vector{Node}
    loci::Vector{Loc}
end


function to_suffixtree(prefix::LongSequence{DNAAlphabet{4}}, d::Dict, ext::Int)
    guides = collect(keys(d))
    loci = collect(values(d))
    lengths = length.(loci)
    cs = cumsum(lengths)
    loci = reduce(vcat, loci)
    stops = cs
    starts = stops .- lengths .+ 1
    loci_range = LociRange.(starts, stops)
    # construct tree nodes
    # we want consecutive guides to be the ones that are splitting
    # this is the order we want for node array
    # we are not guaranteed that a node will have leafs
    #     1
    #  \     \
    #  2     3
    # \  \  \  \
    # 4  5  6  7
    #  \  \ \\ \
    nodes = Vector{Any}([Pair.(guides, loci_range)])
    g_len = length(guides[1]) - ext
    i = 1

    while !isnothing(i)
        # take leftmost from queue that is not a Node
        leaves = nodes[i]
        #@info "leaves: $leaves"
        #parent = popat!(leaves, rand(1:length(leaves)))
        #@info "parent: $parent"
        d_to_p = Base.map(g -> levenshtein(g.first[1:g_len], parent.first, g_len), leaves)
        #@info "leaves: $d_to_p"
        if length(d_to_p) != 0
            radius = balance(d_to_p)
            #@info "radius: $radius"
            inside = leaves[d_to_p .<= radius]
            outside = leaves[d_to_p .> radius]
            # add to the queue future left and right
            if length(inside) > 0
                push!(nodes, inside)
                inside_idx = UInt32(length(nodes))
            else
                inside_idx = UInt32(0)
            end
            if length(outside) > 0
                push!(nodes, outside)
                outside_idx = UInt32(length(nodes))
            else
                outside_idx = UInt32(0)
            end
            nodes[i] = Node(parent.first, parent.second, UInt8(radius), inside_idx, outside_idx)
        else
            nodes[i] = Node(parent.first, parent.second, UInt8(0), UInt32(0), UInt32(0))
        end
        #sleep(60)
        i = findfirst(x -> typeof(x) != Node, nodes)
    end
    nodes = Vector{Node}(nodes)
    return SuffixTreeDB(prefix, nodes, loci)
end


struct TreeDB
    dbi::DBInfo
    prefixes::Set{LongSequence{DNAAlphabet{4}}}
end


"
Build a Vantage Point tree DB of offtargets for the given `motif`,
DB groups off-targets by their prefixes, each prefix has its own
Vantage Point tree.

Will return a path to the database location, same as `storagedir`.

There is an optimization that if the alignment becomes imposible against
the prefix we don't search the off-targets grouped inside the prefix.
Therefore it is advantageous to select larger prefix than maximum 
search distance, however in that case number of files also grows.
"
function build_treeDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storagedir::String,
    prefix_len::Int = 7)

    if prefix_len <= motif.distance
        throw("prefix_len $prefix_len is <= to " * string(motif.distance))
    end

    dbi = DBInfo(genomepath, name, motif)

    # step 1
    @info "Step 1: Searching chromosomes."
    # For each chromsome paralelized we build database
    ref = open(dbi.filepath, "r")
    reader = dbi.is_fa ? FASTA.Reader(ref, index = dbi.filepath * ".fai") : TwoBit.Reader(ref)
    prefixes = Base.map(x -> do_linear_chrom(x, getchromseq(dbi.is_fa, reader[x]), dbi, prefix_len, storagedir), dbi.chrom)
    close(ref)

    prefixes = Set(vcat(prefixes...))

    # step 2
    @info "Step 2: Constructing per prefix db."
    # Iterate over all prefixes and merge different chromosomes
    for prefix in prefixes
        merged = Dict()
        for chrom in dbi.chrom
            p = joinpath(storagedir, string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                merge!(vcat, merged, pdb.suffix)
                rm(p)
            end
        end
        sdb = to_suffixtree(prefix, merged, motif.distance)
        save(sdb, joinpath(storagedir, string(prefix) * ".bin"))
    end

    linDB = TreeDB(dbi, prefixes)
    save(linDB, joinpath(storagedir, "treeDB.bin"))
    @info "Finished constructing treeDB in " * storagedir
    return storagedir
end


function search_prefixtree(
    prefix::LongSequence{DNAAlphabet{4}},
    dist::Int,
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNASeq},
    storagedir::String)
    
    res = zeros(Int, length(guides), dist + 1)
    if detail != ""
        detail_path = joinpath(detail, "detail_" * string(prefix) * ".csv")
        detail_file = open(detail_path, "w")
    end

    # prefix alignment against all the guides
    len_noPAM = length_noPAM(dbi.motif)
    suffix_len =  len_noPAM + dbi.motif.distance - length(prefix)
    isfinal = Base.map(g -> prefix_align(g, prefix, suffix_len, dist).isfinal, guides)

    if all(isfinal)
        return res
    end

    # if any of the guides requires further alignment 
    # load the SuffixDB and iterate
    sdb = load(joinpath(storagedir, string(prefix) * ".bin"))
    for (i, g) in enumerate(guides)
        if !isfinal[i]
            queue = Vector{UInt32}([1])
            while length(queue) > 0
                node = sdb.nodes[popfirst!(queue)]
                dist_i = levenshtein(g, prefix * node.suffix, length(g)) # TODO adjust k, merge prefix-suffix in node?
                if dist_i <= dist
                    res[i, dist_i + 1] += length(node.loci_idx)
                    if detail != ""
                        aln = align(g, prefix * node.suffix, length(g))
                        offtargets = sdb.loci[node.loci_idx.start:node.loci_idx.stop]
                        if dbi.motif.extends5
                            guide_stranded = reverse(g)
                            aln_guide = reverse(aln.guide)
                            aln_ref = reverse(aln.ref)
                        else
                            guide_stranded = g
                            aln_guide = aln.guide
                            aln_ref = aln.ref
                        end
                        noloc = string(guide_stranded) * "," * aln_guide * "," * 
                                aln_ref * "," * string(aln.dist) * ","
                        for offt in offtargets
                            write(detail_file, noloc * decode(offt, dbi) * "\n")
                        end
                    end
                end

                if ((dist_i - dist - 1) <= node.radius) && (node.inside != 0)
                    push!(queue, node.inside)
                end
    
                if ((dist_i + dist + 1) > node.radius)  && (node.outside != 0)
                    push!(queue, node.outside)
                end
            end
        end
    end

    if detail != ""
        close(detail_file)
    end
    return res
end


"
Will search the previously build database for the off-targets of the `guides`. 
Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM

3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide

`dist` - defines maximum levenshtein distance (insertions, deletions, mismatches) for
which off-targets are considered.
`detail` - path and name for the output file. This search will create intermediate 
files which will have same name as detail, but with a sequence prefix. Final file
will contain all those intermediate files. Leave `detail` empty if you are only 
interested in off-target counts returned by the searchDB.
"
function search_treeDB(storagedir::String, guides::Vector{LongDNASeq},  dist::Int = 4; detail::String = "")
    ldb = load(joinpath(storagedir, "treeDB.bin"))
    prefixes = collect(ldb.prefixes)
    if dist > length(first(prefixes)) || dist > ldb.dbi.motif.distance
        error("For this database maximum distance is " * 
              string(min(ldb.dbi.motif.distance, length(first(prefixes)))))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if ldb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    res = zeros(Int, length(guides_), dist + 1)
    for p in prefixes
        res += search_prefixtree(p, dist, ldb.dbi, dirname(detail), guides_, storagedir)
    end
    
    if detail != ""
        # clean up detail files into one file
        open(detail, "w") do detail_file
            write(detail_file, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")
            for prefix in filter(x -> occursin("detail_", x), readdir(dirname(detail)))
                prefix_file = joinpath(dirname(detail), prefix)
                open(prefix_file, "r") do prefix_file
                    for ln in eachline(prefix_file)
                        write(detail_file, ln * "\n")
                    end
                end
                rm(prefix_file)
            end
        end
    end

    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end


# Visualize the tree
function shiftdisp(f::String, o::String, xs::Vector{String})
    rep = repeat([o], length(xs) - 1)
    ch = vcat(f, rep)
    return map(*, ch, xs)
end


function draw_subtrees(tree::SuffixTreeDB, idx::Vector{UInt32}, isinside::Vector{Bool}, levels::Int)
    if length(idx) != 0
        if length(idx) == 1
            return vcat(["│"], shiftdisp("└─ ", "   ", draw_tree(tree, idx[1], isinside[1], levels)))
        else
            return vcat(
                ["│"],
                shiftdisp("├─ ", "│  ", draw_tree(tree, idx[1], isinside[1], levels)),
                draw_subtrees(tree, idx[2:end], isinside[2:end], levels),
            )
        end
    else
        return Vector{String}()
    end
end


function draw_tree(tree::SuffixTreeDB, idx::UInt32, inside::Union{Nothing,Bool}, levels::Int)
    node_pic = "☀ "
    if !isnothing(inside)
            node_pic = inside ? "◐ " : "◑ "
    end
    n = node_pic * string(idx) * " r:" * string(tree.nodes[idx].radius)
    if levels == 0
        return [n]
    else
        levels = levels - 1
        node = tree.nodes[idx]
        leaves = Vector{UInt32}()
        isinside = Vector{Bool}()
        if node.inside != UInt32(0)
            push!(leaves, node.inside)
            push!(isinside, true)
        end
        if node.outside != UInt32(0)
            push!(leaves, node.outside)
            push!(isinside, false)
        end

        subtree = draw_subtrees(tree, leaves, isinside, levels)
        return vcat([n], subtree)
    end
end


function print_treeDB(tree::SuffixTreeDB, start_node::Int = 1, levels::Int = 3, io::IO = stdout)
    println(
        io,
        "Nodes/Guides:",
        length(tree.nodes),
        " Loci:",
        length(tree.loci))

    if start_node > length(tree.nodes) || get(io, :compact, false)
        return nothing
    end

    println(io, "\n" * join(draw_tree(tree, UInt32(start_node), nothing, levels), "\n"))
    return nothing
end


function inspect_treeDB(storagedir::String; levels::Int = 5, inspect_prefix::String = "")
    ldb = load(joinpath(storagedir, "treeDB.bin"))
    prefixes = collect(ldb.prefixes)
    println("Database Info: \n")
    println(ldb.dbi)
    println("\n prefixes: " * string(length(prefixes)) * "\n")
    if inspect_prefix == ""
        inspect_prefix = rand(prefixes)
    end
    sdb = load(joinpath(storagedir, string(inspect_prefix) * ".bin"))
    println("Tree at prefix: " * string(inspect_prefix) * "\n")
    print_treeDB(sdb, 1, levels)
    return nothing
end