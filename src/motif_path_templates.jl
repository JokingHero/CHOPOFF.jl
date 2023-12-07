#struct Path
#    seq::LongDNA{4}
#    dist::Int
#end


struct PathTemplates
    paths::Matrix{Int}
    distances::Vector{Int}
    len::Int # length without the PAM
    max_distance::Int
    mismatch_only::Bool
end

function restrictDistance(mpt::PathTemplates, distance::Int)
    if distance == mpt.max_distance
        return mpt
    elseif distance > mpt.max_distance
        throw("This PathTemplates only support distances up to " * string(mpt.max_distance))
    elseif distance < 0
        throw("Distance can't be below 0.")
    end
    paths_expanded = mpt.paths
    paths_expanded = paths_expanded[:, 1:end - (mpt.max_distance - distance)]
    not_dups = map(!, BitVector(nonunique(DataFrame(paths_expanded, :auto)))) # how can there be no duplicated function?!
    not_over_dist = BitVector(mpt.distances .<= distance)
    not = not_dups .& not_over_dist
    return PathTemplates(paths_expanded[not, :], mpt.distances[not], mpt.len, distance, mpt.mismatch_only)
end

function remove_1_before_non_horizontal!(x::Vector{Int}, base::Vector{Int})
    x_not_in_base = (.!in.(x, Ref(base))) << 1
    x_not_in_base[end] = false
    deleteat!(x, x_not_in_base)
end


function remove_gap!(x::Vector{Int}, gap_idx::Int)
    deleteat!(x, x .== gap_idx)
end


"""
```
adj_matrix_of_guide(
    len::Int, d::Int; 
    mismatch_only::Bool = false)
```

Builds up a shortened version of the alignment graph.
Bases are numbered as: all guide bases (1:len) + Ending times
all the distance we want to extend to, afterwards we add numbering for 
Insertions, Gap and Mismatches - 3 for each distance.
Then for finding all possible alignment with distance 0, you would check
path from node 1 to node (guide length + 1). For distance 1, you would check
all paths from node 1 to (guide length + 1) * 2.

# Arguments
`len` - length of the sequence (e.g. guide)

`d` - Maximal distance on which to build the graph.

`mismatch_only`   -  Whether to skip insertions/deletions.

"""
function adj_matrix_of_guide(len::Int, d::Int; mismatch_only::Bool = false) 
    # notParent is a description for mismatch - it is a base that is not the base of parent node
    l_g = (len + 1) * (d + 1) # guide bases + Ending/Nothing/E - last position
    # + 1 below is last mm from last base to last base on next distance
    l_idm = l_g + (len * 3) * d #  all ins, del, mm - last position

    # fill up all connections
    # horizontal connections (between) guide bases
    adj = zeros(Bool, l_idm, l_idm)
    for di in 1:(d + 1)
        for i in 1:len
            adj[(len + 1) * (di - 1) + i, (len + 1) * (di - 1) + i + 1] = 1
        end
    end
    # vertical connections (between) guide bases - N/Gap/notParent - base on another level
    for di in 1:d
        for i in 1:len
            parent = (len + 1) * (di - 1) + i
            parent_d_next = (len + 1) * di + i
            n = l_g + len * (di - 1) * 3 + (i - 1) * 3 + 1

            if !mismatch_only
                # N
                adj[parent, n] = 1
                adj[n, parent_d_next] = 1

                # Gap
                adj[parent, n + 1] = 1
                adj[n + 1, parent_d_next + 1] = 1
            end

            # notParent
            adj[parent, n + 2] = 1
            adj[n + 2, parent_d_next + 1] = 1
        end
    end
    return adj
end


# guide + not guide + N + Gap + remove last index as it is ending node
# 1:20    21:40       41  42
# 1:len   len+1:len*2 len*2 + 1, len*2 + 2, len*2 + 3

#   DNA    TwoBit
# A 0x01   0x00
# C 0x02   0x01
# G 0x04   0x02
# T 0x08   0x03

# to guide + not guide 1 + not guide 2 + not guide 3 + A C G T
# and expand N
function expand_ambiguous_paths(x::Vector{Int}, len::Int, max_dist::Int)
    y = vcat(x, repeat([len*2 + 1], len + max_dist - length(x)))
    n_pos = [len*4 + 1, len*4 + 2, len*4 + 3, len*4 + 4]
    notguide_pos = [len, len*2, len*3] # + actual guide position

    combinations = []
    mask = BitVector()
    for yi in y
        if yi == (len * 2 + 1) # this is N
            push!(combinations, n_pos)
            push!(mask, true)
        elseif yi in len+1:len*2 # this is not guide position
            push!(combinations, notguide_pos .+ yi .- len)
            push!(mask, true)
        else
            push!(mask, false)
        end
    end

    iter = Iterators.product(combinations...)
    res = transpose(repeat(y, 1, length(iter)))
    for (i, comb) in enumerate(iter)
        res[i, mask] = collect(comb)
    end
    return res
end


"""
```
build_PathTemplates(len::Int, d::Int; storagepath::String = "", mismatch_only::Bool = false)
```

Builds up a PathTemplates object. Stores 
shortened version of all possible paths within each distance `d`
mapped on the graph of all possible alignments of sequence of length
`len`. Then one can use `templates_to_sequences_extended` or 
`templates_to_sequences` and map guide sequence to all possible alignments quickly.

# Arguments
`len` - length of the sequence (e.g. guide - without PAM)

`d` - Maximal distance on which to build the graph.

`storagepath` - If not empty "", will save the object under given path.

`mismatch_only` - Whether to skip insertions/deletions.

"""
function build_PathTemplates(len::Int, d::Int; storagepath::String = "", mismatch_only::Bool = false)
    # path is mapped to these numbers, path numbers are
    # (len + end) * (dist  + 1) and
    # (Ins (N) + Gap + MM) * len * dist
    # and they should be mapped to
    # guide + not guide + N + Gap + remove last index as it is ending node
    # 1:20    21:40       41  42
    # 1:len   len+1:len*2 len*2 + 1, len*2 + 2, len*2 + 3

    adj = adj_matrix_of_guide(len, d; mismatch_only = mismatch_only)
    ngp = repeat([len * 2 + 1, len * 2 + 2, len * 2 + 3], len * d)
    # replace noParents (mismatches - not guide) with proper links to noParents
    for di in 1:d
        for i in 1:len
            ngp[((len) * (di - 1) * 3) + i * 3] = len + i 
        end
    end
    adj_map_to_guide = vcat(repeat(vcat(1:len, 0), d + 1), ngp)

    paths = IdDict{Int64, Vector{Vector{Int64}}}()
    gap_idx = len * 2 + 2
    is_seq_idx = collect(1:len)
    for di in 1:(d + 1)
        pd = path_enumeration(1, (len + 1) * di, adj)
        pd = map(x -> adj_map_to_guide[x.path[1:end-1]], pd)
        # this is to remove 1bp before insertion/mm/gap
        map(x -> remove_1_before_non_horizontal!(x, is_seq_idx), pd)
        map(x -> remove_gap!(x, gap_idx), pd)
        paths[di - 1] = pd
    end

    # now we need to expland the paths so that no ambioguity is left
    # because not guide can have 3 types A -> C T G we will build our guide sequence on top of that
    # new mapping is 
    # guide + not guide 1 + not guide 2 + not guide 3
    # all N can be explanded into all possible positions
    paths_expanded = [mapreduce(x -> expand_ambiguous_paths(x, len, d), vcat, paths[i]) for i in 0:d]
    distances = vcat([repeat([i], size(paths_expanded[i + 1], 1)) for i in 0:d]...)
    paths_expanded = vcat(paths_expanded...)
    not_dups = map(!, BitVector(nonunique(DataFrame(paths_expanded, :auto)))) # how can there be no duplicated function?!
    distances = distances[not_dups]
    paths_expanded = paths_expanded[not_dups, :]

    paths = PathTemplates(paths_expanded, distances, len, d, mismatch_only)
    if storagepath != ""
        save(paths, storagepath)
    end
    return paths 
end



"""
```
build_PathTemplates(motif::Motif; storagepath::String = "", mismatch_only::Bool = false)
```

Builds up a PathTemplates object. Stores 
shortened version of all possible paths for given `Motif`. 
Afterwards use `templates_to_sequences_extended` or 
`templates_to_sequences` and map guide sequence to all possible alignments quickly.

# Arguments
`motif` - Motif object.

`storagepath` - If not empty "", will save the object under given path.

`mismatch_only` - Whether to skip insertions/deletions.

"""
function build_PathTemplates(motif::Motif; storagepath::String = "", mismatch_only::Bool = false)
    len = length_noPAM(motif)
    d = motif.distance
    return build_PathTemplates(len, d; storagepath = storagepath, mismatch_only = mismatch_only)
end

# already in the UInt64 form for folding into one UInt64 as in 
# Base.convert(::Type{UInt64}, x::LongDNA{4}), but in bulk and applied to templates
ALPHABET_TWOBIT = Dict(
    DNA_A => UInt64(0x00), 
    DNA_C => UInt64(0x01), 
    DNA_G => UInt64(0x02),
    DNA_T => UInt64(0x03))

# eeded for searching inside the FM-index
ALPHABET_UINT8 = Dict(
    DNA_A => reinterpret(UInt8, DNA_A), 
    DNA_C => reinterpret(UInt8, DNA_C), 
    DNA_G => reinterpret(UInt8, DNA_G),
    DNA_T => reinterpret(UInt8, DNA_T))


"""
```
guide_to_template_format(
    guide::LongDNA{4}; 
    alphabet::Dict{DNA, Unsigned} = ALPHABET_TWOBIT)
```

Helper that allows you to create mapping vector for the Paths.
Then enumerating possible alignments becomes simple subsetting.

This g_ is a twobitnucleotide mapping.

   DNA    TwoBit
 A 0x01   0x00
 C 0x02   0x01
 G 0x04   0x02
 T 0x08   0x03

 to guide + not guide 1 + not guide 2 + not guide 3 + A C G T

g_[Path_vector]
"""
function guide_to_template_format(guide::LongDNA{4}; 
    alphabet::Dict{DNA, Unsigned} = ALPHABET_TWOBIT)
    
    len = length(guide)
    g_ = repeat([0xff], len * 4 + 4)
    for (i, base) in enumerate(guide)
        if base == DNA_A
            g_[i] = alphabet[DNA_A]
            g_[len * 1 + i] = alphabet[DNA_C]
            g_[len * 2 + i] = alphabet[DNA_G]
            g_[len * 3 + i] = alphabet[DNA_T]
        elseif base == DNA_C 
            g_[i] = alphabet[DNA_C]
            g_[len * 1 + i] = alphabet[DNA_A]
            g_[len * 2 + i] = alphabet[DNA_G]
            g_[len * 3 + i] = alphabet[DNA_T]
        elseif base == DNA_G
            g_[i] = alphabet[DNA_G]
            g_[len * 1 + i] = alphabet[DNA_A]
            g_[len * 2 + i] = alphabet[DNA_C]
            g_[len * 3 + i] = alphabet[DNA_T]
        elseif base == DNA_T 
            g_[i] = alphabet[DNA_T]
            g_[len * 1 + i] = alphabet[DNA_A]
            g_[len * 2 + i] = alphabet[DNA_C]
            g_[len * 3 + i] = alphabet[DNA_G]
        end
    end
    g_[len * 4 + 1] = alphabet[DNA_A]
    g_[len * 4 + 2] = alphabet[DNA_C]
    g_[len * 4 + 3] = alphabet[DNA_G]
    g_[len * 4 + 4] = alphabet[DNA_T]
    return g_
end


"""
```
asUInt64(x::AbstractVecOrMat)
```

Helper that allows you to create one UInt64 for each row out of each PathTemplate.

    matp = guide_in_template_format[pathTemplate]
    map(as64, eachrow(matp))
"""
function asUInt64(x::AbstractVecOrMat)
    y = zero(UInt64)
    for c in x
        y = (y << 2) | UInt64(c)
    end
    mask = (one(UInt64) << (2 * length(x))) - 1
    return reinterpret(UInt64, y & mask)
end



#= all of this is about to be removed I think
"""
```
templates_to_sequences_extended(
    guide::LongDNA{4}, 
    template::PathTemplates;
    dist::Int = template.distance)
```

Uses PathTemplates object - `template` to map
all possible alignments for given `guide` within distance `dist`.
This method expands sequence to the maximal alignment length:
length of the guide + length of the distance. All returned sequences
will be of the same length. The advantage of that is that outputs are unique.

# Arguments
`guide` - guide sequence, without PAM.

`template` - PathTemplates object build for your specific guide queries.

`dist` - Maximal distance on which to return the possible alignments.

# Return

Returns a Vector{Set{LongDNA{4}}} where distance 0 is located at index 1,
distance 1 all possible alignments are located at distance 2 and so on...

"""
function templates_to_sequences_extended( # TODO make this much faster?!
    guide::LongDNA{4}, 
    template::PathTemplates;
    dist::Int = template.distance)
    len = template.len + template.distance

    if length(guide) != template.len
        throw("Wrong guide length.")
    end

    g_ = guide_to_template_format(guide)

    ps = Vector{Set{LongDNA{4}}}()
    for di in 0:dist
        push!(ps, Set(ThreadsX.mapreduce(
            x -> expand_ambiguous(
                LongDNA{4}(g_[x]) * repeat(dna"N", len - length(x))), 
            vcat,
            template.paths[di]; init = Vector{LongDNA{4}}())))
    end
    # if a sequence can exist in lower distance it belongs there rather than higher distance
    # dist 0 is at position 1 in ps, d1 at 2
    for di in 1:dist
        ps[di + 1] = setdiff(ps[di + 1], union(ps[1:di]...))
    end
    return ps
end


"""
```
templates_to_sequences(
    guide::LongDNA{4}, 
    template::PathTemplates;
    dist::Int = template.distance)
```

Uses PathTemplates object - `template` to map
all possible alignments for given `guide` within distance `dist`.
This method does not expand sequences to the maximal alignment length
as opposed to `templates_to_sequences_extended`. This means some sequences might 
seem redundant, for example:
```
For guide "AAA" and Motif("test"), distance 2:

sequence distance 
AAA      0 
AA       1
AAAA     1
     ...
```

# Arguments
`guide` - guide sequence, without PAM.

`template` - PathTemplates object build for your specific guide queries.

`dist` - Maximal distance on which to return the possible alignments.

# Return

Returns a Vector{Path} sorted by the distance, from 0 to `dist`.

"""
function templates_to_sequences(
    guide::LongDNA{4}, 
    template::PathTemplates;
    dist::Int = template.distance)

    if length(guide) != template.len
        throw("Wrong guide length.")
    end

    g_ = guide_to_template_format(guide)

    ps = Vector{Path}()
    for di in 0:dist
        seq = ThreadsX.mapreduce(
            x -> expand_ambiguous(LongDNA{4}(g_[x])), 
            vcat,
            template.paths[di]; init = Vector{LongDNA{4}}())
        seq = ThreadsX.collect(Set(seq))
        append!(ps, ThreadsX.map(x -> Path(x, di), seq))
    end

    # this will return simply seqeuences which can be repeats, D0 in front
    return ThreadsX.sort!(ps, by = p -> (p.dist, p.seq))
end


# this one appends the PAM in forward fashion
function templates_to_sequences(
    guide::LongDNA{4}, 
    template::PathTemplates,
    motif::Motif;
    dist::Int = template.distance)

    if length_noPAM(motif) != template.len
        throw("Length of the motif is not the same as the template!")
    end

    if length(guide) != template.len
        throw("Wrong guide length.")
    end

    g_ = guide_to_template_format(guide)

    ps = Vector{Path}()
    for di in 0:dist
        seq = ThreadsX.mapreduce(
            x -> expand_ambiguous(
                appendPAM_forward(LongDNA{4}(g_[x]), motif)), 
            vcat,
            template.paths[di]; init = Vector{LongDNA{4}}())
        seq = ThreadsX.collect(Set(seq))
        append!(ps, ThreadsX.map(x -> Path(x, di), seq))
    end

    # this will return simply sequences which can be repeats, D0 in front
    return ThreadsX.sort!(ps, by = p -> (p.dist, p.seq))
end
=#










