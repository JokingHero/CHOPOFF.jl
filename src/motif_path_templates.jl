"""
```
PathTemplates(
    paths::Matrix{Int},
    distances::Vector{Int},
    mismatch_only::Bool,
    motif::Motif,
    withPAM::Bool)
```

Contains an expanded version of the alignment graph. Expanded means that all
ambiguous bases are expanded into their basic components, e.g. each sequence containing
NGG PAM will actually be represented 4 times with different base A, C, T, G.

# Arguments
`paths` - contains potential off-target seqeunces in the order of columns: PAM + offtarget sequence + Ns to reach 
 number of bases as maximum distance allows. There are two encodings employed here. First encoding up to the `restricted_len`
 is the UNIT8/2BIT encoding where all bases have to be ACTG, the other encoding is keeping ambigous bases.

`distances` - specifies alignment distances for each potential off-target row

`mismatch_only`   -  Whether insertions/deletions are included.

`motif` - Motif object defining maximal distance, length of the sequence without PAM, length of PAM etc.

`withPAM` - Whether the paths and distances were calcualted with or without PAM.

`restrict_to_len` -  To which length ambiguity should be expanded and after which length ambiguity should be collapsed if possible.
For example: ACTG and ANNN with restriction to length 2, would result in these seqeunces: ACNN AANN AGNN ATNN
This length does not includes PAM - it applies directly to the guide seqeunce. The default is full length of the guide and its maximal distance.
"""
struct PathTemplates
    paths::Matrix{<:Unsigned}
    distances::Vector{UInt8}
    #len::Int # length without the PAM
    #max_distance::Int
    mismatch_only::Bool
    motif::Motif
    withPAM::Bool
    restrict_to_len::Int
end


function restrictDistance(mpt::PathTemplates, distance::Int)
    if distance == mpt.motif.distance
        return mpt
    elseif distance > mpt.motif.distance
        throw("This PathTemplates only support distances up to " * string(mpt.motif.distance))
    elseif distance < 0
        throw("Distance can't be below 0.")
    end

    paths_expanded = mpt.paths
    paths_expanded = paths_expanded[:, 1:(end - (mpt.motif.distance - distance))]
    not_dups = map(!, BitVector(nonunique(DataFrame(paths_expanded, :auto)))) # how can there be no duplicated function?!
    not_over_dist = BitVector(mpt.distances .<= distance)
    not = not_dups .& not_over_dist
    return PathTemplates(mpt.paths[not, :], mpt.distances[not], mpt.mismatch_only, setdist(mpt.motif, distance), mpt.withPAM, mpt.restrict_to_len)
end


function removePAM(mpt::PathTemplates)
    if !mpt.withPAM
        return mpt
    end

    paths_expanded = mpt.paths
    paths_expanded = paths_expanded[:, (length(mpt.motif.pam_loci_fwd) + 1):end]
    not_dups = map(!, BitVector(nonunique(DataFrame(paths_expanded, :auto)))) # how can there be no duplicated function?!
    return PathTemplates(paths_expanded[not_dups, :], mpt.distances[not_dups], mpt.mismatch_only, mpt.motif, false)
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
# and expand ambigous in PAM and extension Ns
function expand_ambiguous_paths(x::Vector{Int}, motif::Motif; 
    restrict_to_len::Int = (length_noPAM(motif) + motif.distance), 
    withPAM::Bool = false)

    len = length_noPAM(motif)
    max_dist = motif.distance
    output_type = smallestutype(unsigned(len*4 + 4))

    #        A          C          G          T
    n_pos = [len*4 + 1, len*4 + 2, len*4 + 3, len*4 + 4]
    notguide_pos = [len, len*2, len*3] # + actual guide position

    y = vcat(x, repeat([len*2 + 1], len + max_dist - length(x))) # add Ns in the ambigous encoding
    restrict_to_len = (restrict_to_len > length(y)) ? length(y) : restrict_to_len
    pam = motif.fwd[motif.pam_loci_fwd]
    if motif.extends5
        reverse!(pam)
    end
    # translate pam to our scheme of guide + not guide 1 + not guide 2 + ng3 + A C G T
    # its not perfect what is below but we need to account for an ambiguity in PAM
    combinations = []
    mask = BitVector()
    pam_translated = Vector{Int}()
    # COV_EXCL_START
    if withPAM
        for p in pam
            if p == DNA_A
                push!(pam_translated, n_pos[1])
                push!(mask, false)
            elseif p == DNA_C
                push!(pam_translated, n_pos[2])
                push!(mask, false)
            elseif p == DNA_G
                push!(pam_translated, n_pos[3])
                push!(mask, false)
            elseif p == DNA_T
                push!(pam_translated, n_pos[4])
                push!(mask, false)
            elseif p == DNA_N
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos)
                push!(mask, true)
            elseif p == DNA_Gap
                continue
            elseif p == DNA_B
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[2:4])
                push!(mask, true)
            elseif p == DNA_D
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[1, 3, 4]]) # A G T 
                push!(mask, true)
            elseif p == DNA_H
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[1, 2, 4]]) # A C T 
                push!(mask, true)
            elseif p == DNA_K
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[3, 4]]) # A G T 
                push!(mask, true)
            elseif p == DNA_M
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[1, 2]]) # A C 
                push!(mask, true)
            elseif p == DNA_R
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[1, 3]]) # A G 
                push!(mask, true)
            elseif p == DNA_S
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[2, 3]]) # C G 
                push!(mask, true)
            elseif p == DNA_V
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[1, 2, 3]]) # A C G 
                push!(mask, true)
            elseif p == DNA_W
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[1, 2, 4]]) # A T
                push!(mask, true)
            elseif p == DNA_Y
                push!(pam_translated, 0) # WILL be relaced by combinations
                push!(combinations, n_pos[[2, 4]]) # C T
                push!(mask, true)
            end
        end
    end
    # COV_EXCL_STOP

    for yi in y[1:restrict_to_len]
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
    mask = vcat(mask, BitVector(repeat([false], length(y) - restrict_to_len)))

    if withPAM
        y = vcat(pam_translated, y)
    end
    
    iter = Iterators.product(combinations...)
    res = permutedims(repeat(y, 1, length(iter)))
    for (i, comb) in enumerate(iter)
        res[i, mask] = collect(comb)
    end
    return convert.(output_type, res)
end


"""
```
build_PathTemplates(
    motif::Motif; 
    storagepath::String = "", 
    mismatch_only::Bool = false, 
    restrict_to_len::Int = length_noPAM(motif),
    withPAM::Bool = false)
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

`restrict_to_len` - To which length ambiguity should be expanded and after which length ambiguity should be collapsed if possible.
    For example: ACTG and ANNN with restriction to length 2, would result in these seqeunces: ACNN AANN AGNN ATNN
    This length does not includes PAM - it applies directly to the guide seqeunce. The default is full length of the guide and its maximal distance.

`withPAM` - Whether to include PAM in the paths. Default is false.

"""
function build_PathTemplates(
    motif::Motif; 
    storagepath::String = "", 
    mismatch_only::Bool = false, 
    restrict_to_len::Int = (length_noPAM(motif) + motif.distance),
    withPAM::Bool = false,
    reuse_saved::Bool = true)

    len = length_noPAM(motif)
    d = motif.distance
    length_of_paths = (withPAM ? length(motif) : length_noPAM(motif)) + motif.distance

    # trying to preload paths for most common use case
    if (reuse_saved && (motif.distance <= 4) && (restrict_to_len <= 16) && !withPAM && !mismatch_only)
        m2 = Motif("Cas9")
        if (motif.fwd == m2.fwd && 
            motif.rve == m2.rve &&
            motif.pam_loci_fwd == m2.pam_loci_fwd && 
            motif.pam_loci_rve == m2.pam_loci_rve)

            dir = joinpath(dirname(pathof(CHOPOFF)), "..", "data")
            pfile = joinpath(dir, "Cas9_d4_p16_paths.bin")
            dfile = joinpath(dir, "Cas9_d4_p16_distances.bin")
            if (isfile(pfile) && isfile(dfile))
                @info "Reusing precomputed alignments."
                paths = CHOPOFF.load(pfile)
                distances = CHOPOFF.load(joinpath(dir, "Cas9_d4_p16_distances.bin"))
                paths = paths[:, 1:restrict_to_len]
                not_dups = map(!, BitVector(nonunique(DataFrame(paths, :auto))))
                not_over_dist = BitVector(distances .<= d)
                not = not_dups .& not_over_dist
                paths = paths[not, :]
                distances = distances[not]
                paths = convert.(smallestutype(maximum(paths)), paths)
                paths = PathTemplates(paths, distances, mismatch_only, motif, withPAM, restrict_to_len)
                if storagepath != ""
                    save(paths, storagepath)
                end
                return paths 
            end
        end
    end

    # path is mapped to these numbers, path numbers are
    # (len + end) * (dist  + 1) and
    # (Ins (N) + Gap + MM) * len * dist
    # and they should be mapped to
    # guide + not guide + N + Gap + remove last index as it is ending node
    # 1:20    21:40       41  42
    # 1:len   len+1:len*2 len*2 + 1, len*2 + 2, len*2 + 3
    adj = CHOPOFF.adj_matrix_of_guide(len, d; mismatch_only = mismatch_only)
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
        map(x -> CHOPOFF.remove_1_before_non_horizontal!(x, is_seq_idx), pd)
        map(x -> CHOPOFF.remove_gap!(x, gap_idx), pd)
        paths[di - 1] = pd
    end

    # now we need to expland the paths so that no ambioguity is left
    # because not guide can have 3 types A -> C T G we will build our guide sequence on top of that
    # new mapping is 
    # guide + not guide 1 + not guide 2 + not guide 3
    # all N can be explanded into all possible positions
    paths_expanded = Base.mapreduce(x -> expand_ambiguous_paths(x, motif; restrict_to_len = restrict_to_len, withPAM = withPAM), vcat, paths[0]; 
        init = Matrix{UInt8}(undef, 0, length_of_paths))
    paths_lengths = [size(paths_expanded)[1]]
    for i in 1:d
        paths_i = ThreadsX.mapreduce(x -> expand_ambiguous_paths(x, motif; restrict_to_len = restrict_to_len, withPAM = withPAM), vcat, paths[i]; 
            init = Matrix{UInt8}(undef, 0, length_of_paths))
        push!(paths_lengths, size(paths_i)[1])
        paths_expanded = vcat(paths_expanded, paths_i)
        paths_i = nothing
    end
    paths = nothing
    not_dups = map(!, BitVector(nonunique(DataFrame(paths_expanded, :auto)))) # how can there be no duplicated function?!
    paths_expanded = paths_expanded[not_dups, :] # filters out around 2%

    distances = vcat([repeat([convert(UInt8, i)], paths_lengths[i + 1]) for i in 0:d]...)
    distances = distances[not_dups]
    
    paths_expanded = convert.(smallestutype(maximum(paths_expanded)), paths_expanded)
    distances = UInt8.(distances)
    paths = PathTemplates(paths_expanded, distances, mismatch_only, motif, withPAM, restrict_to_len)
    if storagepath != ""
        save(paths, storagepath)
    end
    return paths 
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
It assumes guides are in PAM - Ns config here. For antisense we 
have to complement the PAM and guide seqeunce here.

This g_ is a twobitnucleotide mapping.

   DNA    TwoBit
 A 0x01   0x00
 C 0x02   0x01
 G 0x04   0x02
 T 0x08   0x03

 to guide + not guide 1 + not guide 2 + not guide 3 + A C G T

g_[Path_vector]
"""
function guide_to_template_format(guide::LongDNA{4}, is_antisense = false; 
    alphabet::Union{Dict{DNA, UInt8}, Dict{DNA, UInt64}} = ALPHABET_TWOBIT)

    if is_antisense
        guide = complement(guide)
    end
    
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
    g_[len * 4 + 1] = is_antisense ? alphabet[DNA_T] : alphabet[DNA_A]
    g_[len * 4 + 2] = is_antisense ? alphabet[DNA_G] : alphabet[DNA_C]
    g_[len * 4 + 3] = is_antisense ? alphabet[DNA_C] : alphabet[DNA_G]
    g_[len * 4 + 4] = is_antisense ? alphabet[DNA_A] : alphabet[DNA_T]
    return g_
end


function guide_to_template_format_ambig(guide::LongDNA{4}, is_antisense = false)

    if is_antisense
        guide = complement(guide)
    end
    
    g_ = copy(guide)
    for (i, base) in enumerate(guide)
        if base == DNA_A
            not_base = DNA_B
        elseif base == DNA_C 
            not_base = DNA_D
        elseif base == DNA_G
            not_base = DNA_H
        elseif base == DNA_T 
            not_base = DNA_V
        end
        push!(g_, not_base)
    end
    push!(g_, DNA_N)
    push!(g_, DNA_Gap)
    return collect(g_)
end


"""
```
asUInt(x::AbstractVecOrMat)
```

Helper that allows you to create one UInt for DNA strings smaller.
There is no checking for the size of the vector.

    matp = guide_in_template_format[pathTemplate]
    map(x -> asUInt(UInt32, x), eachrow(sa))
"""
function asUInt(t::Type{<:Unsigned}, x::AbstractVecOrMat)
    y = zero(t)
    for c in x
        y = (y << 2) | convert(t, c)
    end
    mask = (one(t) << (2 * length(x))) - t(1)
    return reinterpret(t, y & mask)
end


"""
```
duplicated(x::Vector{UInt64})
```

Helper that allows to find which values occur more than once in the Vector. 
Returns BitVector of duplciated positions.
"""
function duplicated(x::Vector{<:Number})
    s = Set(Vector{eltype(x)}())
    b = BitVector(zeros(length(x)))
    for (i, xi) in enumerate(x)
        if xi in s
            b[i] = true
        else
            push!(s, xi)
        end
    end
    return b
end