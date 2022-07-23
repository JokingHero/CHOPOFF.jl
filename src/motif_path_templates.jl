struct Path
    seq::LongDNA{4}
    seq_len::Int
    dist::Int
    reducible::Int # index inside the outer vector linking to which 
    # Path also will contain this lower distance Path
    # 0 means no reducible
    # reducible points to the sequence which already contains
end

# AAA 4
# AA  3 - lower dist already contains larger - no such cases
# AAAC  5, 1 - higher distance is more specific - this value has to be unclaimed at reducible!
# we need to still add it to the correct


struct MotifPathTemplates
    paths::IdDict{Int64, Vector{Vector{Int64}}}
    motif::Motif
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
`len` - length of the seqeunce (e.g. guide)

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


function build_motifTemplates(motif::Motif; storagepath::String = "", mismatch_only::Bool = false)
    len = length_noPAM(motif)
    d = motif.distance

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

    paths = MotifPathTemplates(paths, motif)
    if storagepath != ""
        save(paths, storagepath)
    end
    return paths 
end


function templates_to_sequences(
    guide::LongDNA{4}, 
    template::CRISPRofftargetHunter.MotifPathTemplates;
    dist::Int = length(template) - 1,
    reducible::Bool = true)

    if length(guide) != CRISPRofftargetHunter.length_noPAM(template.motif)
        throw("Wrong guide length.")
    end

    g_ = copy(guide)
    for (i, base) in enumerate(guide)
        if base == DNA_A
            not_base = DNA_B
        elseif base == DNA_C 
            not_base = DNA_D
        elseif base == DNA_T
            not_base = DNA_V
        elseif base == DNA_G 
            not_base = DNA_H
        end
        push!(g_, not_base)
    end
    push!(g_, DNA_N)
    push!(g_, DNA_Gap)
    g_ = collect(g_)

    ps = Vector{Path}()
    for di in 0:dist
        seq = ThreadsX.mapreduce(
            x -> CRISPRofftargetHunter.expand_ambiguous(LongDNA{4}(g_[x])), 
            vcat,
            template.paths[di]; init = Vector{LongDNA{4}}())
        seq = ThreadsX.collect(Set(seq))
        append!(ps, ThreadsX.map(x -> Path(x, length(x), di, 0), seq))
    end

    if !reducible # this will return simply seqeuences which can be repeats, D0 in front
        return ThreadsX.sort!(ps, by = p -> (p.dist, p.seq))
    end
    ThreadsX.sort!(ps, by = p -> (p.seq, p.dist))

    to_remove = zeros(Bool, length(ps))
    reducible = zeros(Int, length(ps))
    ps_len = length(ps)
    for (i, p) in enumerate(ps)
        j = i + 1
        before_i = sum(to_remove[1:i])
        while j <= ps_len &&
            p.seq_len <= ps[j].seq_len && 
            isequal(p.seq, ps[j].seq[1:p.seq_len]) # if bigger - stop!
            
            if !to_remove[j]
                reducible[j] = i - before_i
                if ps[j].dist >= p.dist
                    to_remove[j] = true
                end
            end
            j += 1
        end
    end

    deleteat!(ps, to_remove)
    deleteat!(reducible, to_remove)
    ps = ThreadsX.map(x -> Path(x[1].seq, x[1].seq_len, x[1].dist, x[2]), zip(ps, reducible))
    return ps
end











