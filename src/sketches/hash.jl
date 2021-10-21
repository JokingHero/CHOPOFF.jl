using BioSequences
using CRISPRofftargetHunter
using ThreadsX
using Statistics

"
Implementation of hashing function that cares about uniform distribution of the data.

FNV doc: http://www.isthe.com/chongo/tech/comp/fnv/#FNV-1a
"
function fnv_hash_1a(key::T) where T <: BioSequence
    h = 0x811c9dc5 #0xcbf29ce484222325
    for k in key
        @inbounds h = (h ⊻ reinterpret(UInt8, k)) * 0x01000193 # 0x100000001b3
    end
   return h
end


function fnv_hash_1a(key::UInt32)
    key = reinterpret(UInt8, [key])
    h = 0x811c9dc5 # 0xcbf29ce484222325
    for k in key
        @inbounds h = (h ⊻ reinterpret(UInt8, k)) * 0x01000193 #0x100000001b3
    end
   return h
end

"
Lazy mode mapping.
Yield arbitrary range 1-n.
Bias against larger values.
"
function modulo_range(h::K, n::K)  where K <: Unsigned
    return h % n
end


"
Retry method.
"
function retry_range(h::UInt32, n::UInt32)
    retry_level = (typemax(h) / n) * n
    while h >= retry_level
        h = fnv_hash_1a(h)
    end
    return h % n
end

# 

function fnv_hash(key::T, seed::UInt32, n::UInt32) where T <: BioSequence
    h = fnv_hash_1a(key) ⊻ fnv_hash_1a(seed)
    return retry_range(h, n)
end


"
Murmur 64 hash implmentation. Simplified for 
LongDNASeq, here every element of key is UInt8, but 
really encodes only 4 bits.
"
function murmur_hash_64a(key::T, seed::UInt64) where T <: BioSequence
    m = 0xc6a4a7935bd1e995
    r = 47
    h = seed ⊻ m
    for k in key
        k = reinterpret(UInt8, k)
        k *= m
        k ⊻= k >> r 
        k *= m
        h ⊻= k
        h *= m 
    end
 
    h ⊻= h >> r
    h *= m
    h ⊻= h >> r
    return h
end


function murmur_hash(key::T, seed::UInt64, n::UInt64) where T <: BioSequence
    return modulo_range(murmur_hash_64a(key, seed), n)
end


# read in all guides
db_path = joinpath("/home/ai/", "tests_hg38v34")
guides = CRISPRofftargetHunter.load(joinpath(db_path, "uniq_guides_as_MER20.bin"))
open("guides.txt", "w") do f
    # Make sure we write 64bit integer in little-endian byte order
    for i in guides
        write(f, string(BioSequences.encoded_data(i)) * "\n")
    end
end


#guides2 = ThreadsX.map(LongDNASeq, guides)
n = length(guides)
table_height = ceil(UInt32, 1.23 * n)

s1 = rand(UInt32)
s2 = rand(UInt32)
s3 = rand(UInt32)

# simulate 3 selects
h1 = ThreadsX.map(x -> fnv_hash(x, s1, table_height), guides)
h12 = ThreadsX.map(x -> fnv_hash(x, s2, table_height), guides)
h13 = ThreadsX.map(x -> fnv_hash(x, s3, table_height), guides)


# simulate first hash for selection where, 3 followup hashes for 
# selection in the neighbourhood w
function window_hash(direction::Bool, origin::UInt32, shift::UInt32, table_height::UInt32)
    if direction 
        if (origin + shift) < table_height 
            return (origin + shift)
        end
        return (origin - shift)
    else 
        if ((origin - shift) > 1) 
            return (origin - shift)
        end
        return (origin + shift)
    end
end

w = UInt32(101) # use some prime number
hw = ThreadsX.map(x -> fnv_hash(x, rand(UInt32), table_height), guides) # window
d1 = rand(Bool, length(guides)) # directionality
h1 = ThreadsX.map(x -> fnv_hash(x, rand(UInt32), w), guides)
d12 = rand(Bool, length(guides))
h12 = ThreadsX.map(x -> fnv_hash(x, rand(UInt32), w), guides)
d13 = rand(Bool, length(guides))
h13 = ThreadsX.map(x -> fnv_hash(x, rand(UInt32), w), guides)

hw1 = ThreadsX.map(x -> window_hash(d1[x], hw[x], h1[x], table_height), 1:length(guides))
hw2 = ThreadsX.map(x -> window_hash(d12[x], hw[x], h12[x], table_height), 1:length(guides))
hw3 = ThreadsX.map(x -> window_hash(d13[x], hw[x], h13[x], table_height), 1:length(guides))


h2 = ThreadsX.map(x -> murmur_hash(x, 1, table_height), guides)

index(len, h) = reinterpret(Int, Core.Intrinsics.urem_int(h, reinterpret(UInt64, len))) + 1
h3 = ThreadsX.map(x -> index(table_height, hash(BioSequences.encoded_data(x), UInt64(1))), guides)
h4 = ThreadsX.map(x -> index(table_height, hash(LongDNASeq(x), UInt64(1))), guides) # slowest by far

function judge(h, table_height)
    sort!(h)
    hash1, c1 = CRISPRofftargetHunter.ranges(h)

    sh1 = Set(hash1)
    sh_all = Set(1:table_height)

    # first measure is how many are empty slots
    empty_slots = setdiff(sh_all, sh1)
    c1 = length.(c1)
    return "ES: " * string(length(empty_slots) / table_height) * 
        "  Median: " * string(median(c1)) *
        "  Mean:  " * string(mean(c1)) * 
        "  Max:  " * string(maximum(c1))
end

judge(hw1, table_height)
# "ES: 0.4435771475607137  Median: 1.0  Mean:  1.461133411341778  Max:  10"
h = vcat(hw1, hw2, hw3)
hw1 = 0
hw2 = 0
hw3 = 0
judge(h, table_height)

judge(h1, table_height)
# "ES: 0.44  Median: 1.0  Mean:  1.46  Max:  10"
judge(h12, table_height)
judge(h13, table_height)

h = vcat(h1, h12, h13)
h1 = 0
h12 = 0
h13 = 0
judge(h, table_height)
# ES: 0.08724608466149524  Median: 2.0  Mean:  2.6721598629928374  Max:  16

# "ES: 0.44  Median: 1.0  Mean:  1.46  Max:  10"
judge(h2, table_height)
# "ES: 0.44  Median: 1.0  Mean:  1.46  Max:  10"
judge(h3, table_height)
# "ES: 0.44  Median: 1.0  Mean:  1.46  Max:  10
judge(h4, table_height)



