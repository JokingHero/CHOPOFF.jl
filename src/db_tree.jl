"
Each node contains one guide (suffix) and 
indexes of that guide loci (loci_idx).
This guide also contains a radius and pointers
to other Nodes that are inside and outside the 
radius.
"
struct Node
    suffix::LongDNA{4}
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
    prefix::LongDNA{4}
    nodes::Vector{Node}
    loci::Vector{Loc}
end


function to_suffixtree(prefix::LongDNA{4}, 
    guides::Vector{LongDNA{4}}, 
    loci::Vector{Loc}, ext::Int)
    
    (guides, loci_range, loci) = unique_guides(guides, loci)
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
    # stop changing type and length of nodes
    n_guides = length(guides)
    g_len = length(guides[1]) - ext
    order = zeros(UInt32, n_guides) # of guides and loci_range
    order[1] = 1 # we take the first guide as leading guide
    order_idx = 1 # keep track of last filled value
    inside = zeros(UInt32, n_guides)
    outside = zeros(UInt32, n_guides)
    radius = zeros(UInt8, n_guides)
    parent = ones(UInt32, n_guides)
    parent[1] = 0

    for i in 1:n_guides # represents slot of the order, order[i] will be parent in this iteration
        # Take all guides that have i'th guide as parent
        #@info "i is: $i"
        #@info "leaves of: " * string(order[i])
        leaves_idx = ThreadsX.findall(isequal(i), parent)
        d_to_p = ThreadsX.map(l -> levenshtein(guides[l][1:g_len], guides[order[i]], g_len), leaves_idx)
        #@info "$d_to_p"
        if length(d_to_p) != 0
            r = balance(d_to_p)
            radius[i] = r
            inside_idx = leaves_idx[d_to_p .<= r]
            outside_idx = leaves_idx[d_to_p .> r]
            if length(inside_idx) > 0
                order_idx += 1
                inside[i] = order_idx
                order[order_idx] = pop!(inside_idx) # future parent
                parent[inside_idx] .= order_idx # future children
            end
            if length(outside_idx) > 0
                order_idx += 1
                outside[i] = order_idx
                order[order_idx] = pop!(outside_idx) # future parent
                parent[outside_idx] .= order_idx # future children
            end
        end
    end

    guides = guides[order]
    loci_range = loci_range[order]
    nodes = ThreadsX.map(x -> Node(guides[x], loci_range[x], radius[x], inside[x], outside[x]), 1:n_guides)
    return SuffixTreeDB(prefix, nodes, loci)
end


struct TreeDB
    dbi::DBInfo
    prefixes::Set{LongDNA{4}}
end


"""
```
build_treeDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    prefix_len::Int = 7)
```  

Build a Vantage Point tree DB of offtargets for the given `motif`,
DB groups off-targets by their prefixes, each prefix has its own
Vantage Point tree.

Will return a path to the database location, same as `storage_dir`.

There is an optimization that if the alignment becomes imposible against
the prefix we don't search the off-targets grouped inside the prefix.
Therefore it is advantageous to select larger prefix than maximum 
search distance, however in that case number of files also grows.

# Arguments

`name` - Your prefered name for this index for easier identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif deines what kind of gRNA to search for.

`storage_dir`  - Folder path to the where index will be saved with name `linearDB.bin` and many prefix files.

`prefix_len`  - Size of the prefix by which off-targets are indexed. Prefix of 8 or larger will be the fastest,
                however it will also result in large number of files. 

# Examples
```julia-repl
$(make_example_doc("treeDB"))
```
"""
function build_treeDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    prefix_len::Int = 7)

    if prefix_len <= motif.distance
        throw("prefix_len $prefix_len is <= to " * string(motif.distance))
    end

    dbi = DBInfo(genomepath, name, motif)

    # step 1
    @info "Step 1: Searching chromosomes."
    # For each chromsome paralelized we build database
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
    prefixes = Base.map(x -> do_linear_chrom(x, getchromseq(dbi.gi.is_fa, reader[x]), dbi, prefix_len, storage_dir), dbi.gi.chrom)
    close(ref)

    prefixes = Set(vcat(prefixes...))

    # step 2
    @info "Step 2: Constructing per prefix db."
    # Iterate over all prefixes and merge different chromosomes
    i = 0
    for prefix in prefixes
        guides = Vector{LongDNA{4}}()
        loci = Vector{Loc}()
        for chrom in dbi.gi.chrom
            p = joinpath(storage_dir, string(prefix), string(prefix) * "_" * chrom * ".bin")
            if ispath(p)
                pdb = load(p)
                append!(guides, pdb.suffix)
                append!(loci, pdb.loci)
            end
        end
        rm(joinpath(storage_dir, string(prefix)), recursive = true)
        sdb = to_suffixtree(prefix, guides, loci, motif.distance)
        save(sdb, joinpath(storage_dir, string(prefix) * ".bin"))
        i += 1
        #@info "Done prefixes: " * string(round(i / length(prefixes); digits = 2) * 100)
    end

    linDB = TreeDB(dbi, prefixes)
    save(linDB, joinpath(storage_dir, "treeDB.bin"))
    @info "Finished constructing treeDB in " * storage_dir
    return storage_dir
end


function search_prefixtree(
    prefix::LongDNA{4},
    dist::Int,
    dbi::DBInfo,
    detail::String,
    guides::Vector{LongDNA{4}},
    storage_dir::String)

    # prefix alignment against all the guides
    len_noPAM = length_noPAM(dbi.motif)
    suffix_len =  len_noPAM + dbi.motif.distance - length(prefix)
    isfinal = Base.map(g -> prefix_align(g, prefix, suffix_len, dist).isfinal, guides)

    if all(isfinal)
        return 
    end

    detail_path = joinpath(detail, "detail_" * string(prefix) * ".csv")
    detail_file = open(detail_path, "w")

    # if any of the guides requires further alignment 
    # load the SuffixDB and iterate
    sdb = load(joinpath(storage_dir, string(prefix) * ".bin"))
    for (i, g) in enumerate(guides)
        if !isfinal[i]
            queue = Vector{UInt32}([1])
            while length(queue) > 0
                node = sdb.nodes[popfirst!(queue)]
                # we only need to know whether to search inside/outside
                dist_i, ref_left = levenshtein2(g, prefix * node.suffix, node.radius + dist)
                ext = 0
                if dbi.motif.distance > ref_left
                    ext = dbi.motif.distance - ref_left
                end

                if dist_i <= dist
                    aln = align(g, prefix * node.suffix, dist)
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

                
                # ext is not zero only when ext is used for alignment between VP and guide
                if (node.inside != 0) && ((dist_i - dist - ext - 1) <= node.radius)
                    push!(queue, node.inside)
                end
    
                if (node.outside != 0) && ((dist_i + dist + ext + 1) > node.radius)
                    push!(queue, node.outside)
                end
            end
        end
    end

    close(detail_file)
    return 
end


"""
```
search_treeDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3)
```

Search previously build treeDB database for the off-targets of the `guides`. 

Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

```
5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM

3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide
```

# Arguments

`dist` - Defines maximum levenshtein distance (insertions, deletions, mismatches) for 
which off-targets are considered.  

`detail` - Path and name for the output file. This search will create intermediate 
files which will have same name as detail, but with a sequence prefix. Final file
will contain all those intermediate files. Leave `detail` empty if you are only 
interested in off-target counts returned by the treeDB. 


# Examples
```julia-repl
$(make_example_doc("treeDB"))
```
"""
function search_treeDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3)
    ldb = load(joinpath(storage_dir, "treeDB.bin"))
    prefixes = collect(ldb.prefixes)
    if distance > length(first(prefixes)) || distance > ldb.dbi.motif.distance
        error("For this database maximum distance is " * 
              string(min(ldb.dbi.motif.distance, length(first(prefixes)))))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if ldb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    ThreadsX.map(p -> search_prefixtree(p, distance, ldb.dbi, dirname(output_file), guides_, storage_dir), prefixes)
    
    cleanup_detail(output_file)
    return
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


"""
```
inspect_treeDB(
    storage_dir::String; 
    levels::Int = 5, 
    inspect_prefix::String = "")
```

See small part of the full vantage point tree of the treeDB.

TreeDB can be split based on the distance to the radius (r) 
into inside (left <= r) and right (outside > r) nodes.

# Examples
```julia-repl
$(make_example_doc("treeDB"))

# finally, view some part of the database!
inspect_treeDB(db_path; inspect_prefix = "CCGTCGC")
```
"""
function inspect_treeDB(storage_dir::String; levels::Int = 5, inspect_prefix::String = "")
    ldb = load(joinpath(storage_dir, "treeDB.bin"))
    prefixes = collect(ldb.prefixes)
    println("Database Info: \n")
    println(ldb.dbi)
    println("\n prefixes: " * string(length(prefixes)) * "\n")
    if inspect_prefix == ""
        inspect_prefix = rand(prefixes)
    end
    sdb = load(joinpath(storage_dir, string(inspect_prefix) * ".bin"))
    println("Tree at prefix: " * string(inspect_prefix) * "\n")
    print_treeDB(sdb, 1, levels)
    return nothing
end