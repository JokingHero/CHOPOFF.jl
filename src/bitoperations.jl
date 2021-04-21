# using Combinatorics
# using BenchmarkTools
# using BioSymbols
# using BioSequences

function l(x::Int, shift::Int)
    return (x - 1) * shift + 1
end

function r(x::Int, shift::Int)
    return x * shift
end

"
Delete `start` to `stop` bits and shift left
side of that range to the left. Right side of
deleted seqeunce has to stay the same.

`start`: inclusive position of the deletion start counting from right
`stop`: inclusive position of the deletion end counting from right
"
function deleterange(x::UInt64, start, stop)
    left = (x >> stop) << (start - 1)
    right = (x << (64 - start + 1)) >> (64 - start + 1)
    return left | right
end

function deleterange(x::UInt128, start, stop)
    left = (x >> stop) << (start - 1)
    right = (x << (128 - start + 1)) >> (128 - start + 1)
    return left | right
end

"
Will enlist an array of deletion permutations tuple,
where always first value is start of the deletion and second is
end of the deletion, inclusive.

`len`: Length of the delete vector, it assumes each position in
       `len` contains `shift` long encoding, so maximum length
       in the output is `len` * `shift`
`k`: Maximum positions that can be deleted at a time inside `len`.
`shift`: Encoding length that each position of `len` is carrying.
"
function deletion_permutations(len::Int, k::Int = 4, shift::Int = 4)
    del = Vector{Array{Int64,1}}()
    for i in 1:k
        append!(del, collect(combinations(1:len, i)))
    end

    del_ = Vector{Vector{UnitRange{Int}}}()
    for v in del
        v_ = Vector{UnitRange{Int}}()
        if length(v) == 1
            push!(del_, [UnitRange(Int(l(v[1], shift)), Int(r(v[1], shift)))])
            continue
        end
        tuple_start = 0
        skip = false
        for d in 2:length(v)
            if !skip
                tuple_start = v[d - 1]
            end
            if (v[d - 1] + 1) == v[d]
                skip = true
            else
                push!(v_, UnitRange(Int(l(tuple_start, shift)), Int(r(v[d - 1], shift))))
                tuple_start = v[d]
            end
        end
        push!(v_, UnitRange(Int(l(tuple_start, shift)), Int(r(v[end], shift))))
        sort!(v_, rev=true)
        push!(del_, v_)
    end

    return del_
end


"
Translate data encoded in 64 format into 128 format of DNA encodings.
Can work only on DNA sequences, transforming them into their bit encodings.
"
function translate64_to_128(x::UInt64)
    x = Mer{DNAAlphabet{2}, 20}(x)
    x = LongDNASeq(x)
    x = BioSequences.encoded_data(x)
    return parse(UInt128, bitstring(x[2]) * bitstring(x[1]); base = 2)
end

"
Treat input UInt128 as encoded DNA sequence and translate it
into UInt64 encoding.
"
function translate128_to_64(x::UInt128)
    x = bitstring(x)
    x_seq = LongDNASeq("")
    for i in 1:4:length(x)
        xi = reinterpret(DNA, parse(UInt8, x[i:i + 3]; base = 2))
        if xi != DNA_Gap
            push!(x_seq, xi)
        end
    end
    x_seq = reverse(x_seq)
    x_seq = Mer{DNAAlphabet{2}, 20}(x_seq)
    return BioSequences.encoded_data(x_seq)
end

# TEST translate work
# x = Mer{DNAAlphabet{2}, 20}("ACTGACTGACTGACTGACTG")
# y = BioSequences.encoded_data(x)
# y = translate64_to_128(y)
# y = translate128_to_64(y)
# Mer{DNAAlphabet{2}, 20}(y) == x

function translate_to_dna(x::UInt128)
    x = bitstring(x)
    x_seq = LongDNASeq("")
    for i in 1:4:length(x)
        xi = reinterpret(DNA, parse(UInt8, x[i:i + 3]; base = 2))
        if xi != DNA_Gap
            push!(x_seq, xi)
        end
    end
    x_seq = reverse(x_seq)
    return x_seq
end

# for chunk of a data which is Array of UInt128
function generate_deletes(guides::Vector{UInt128}, positions::Vector{Array{Tuple{Int, Int}, 1}})
    dict = IdDict{UInt128, Vector{UInt128}}()
    for g in guides
        for pos in 1:length(positions)
            y = g
            for  r in 1:length(positions[pos])
                start, stop = positions[pos][r]
                y = deleterange(y, start, stop)
            end
            roots = get(dict, y, [])
            dict[y] = vcat(roots, g)
        end
    end

    for k in keys(dict)
        dict[k] = unique(dict[k])
    end
    return dict
end

function generate_deletes(guide::UInt128, positions::Vector{Array{Tuple{Int, Int}, 1}})
    deletes::Vector{UInt128} = fill(0, length(positions))
    for pos in 1:length(positions)
        y = guide
        for  r in 1:length(positions[pos])
            start, stop = positions[pos][r]
            y = deleterange(y, start, stop)
        end
        deletes[pos] = y
    end
    return unique(deletes)
end


# read dict of keys
# cd("/home/ai/.julia/dev/CRISPRofftargetHunter")
# cd("/export/valenfs/projects/Kornel/offtarget_aligner")
# @load "./hg38v34_hmap_g_count.jld2" hmap_g_count
#
# keys_hmap = collect(keys(hmap_g_count)) #[1:10000]
# keys_hmap = [translate64_to_128(x) for x in keys_hmap]
# const del_ = deletion_permutations(20, 4)
#
# deletes_vec = fill(Vector{UInt128}(), length(keys_hmap))
# Threads.@threads for (i, v) in collect(enumerate(keys_hmap))
#     deletes_vec[i] = generate_deletes(v, del_)
# end
# println("Save deletions vector")
# @save "./hg38v34_hmap_deletes_vec.jld2" deletes_vec
#
# # transform deletions vector into Dict
# dict = IdDict{UInt128, Vector{UInt128}}()
# for i in 1:length(keys_hmap)
#     println(i)
#     for del in deletes_vec[i]
#         if !haskey(dict, del)
#             find_first = fill(false, length(keys_hmap))
#             def_any = x -> x == del
#             Threads.@threads for (di, dv) in collect(enumerate(deletes_vec))
#                 find_first[di] = any(def_any, dv)
#             end
#             dict[del] = keys_hmap[find_first]
#         end
#     end
# end
