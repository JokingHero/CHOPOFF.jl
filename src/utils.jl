"
Randomize sequence of length `n` from `letters`.
"
function getSeq(n = 20, letters = ['A', 'C', 'G', 'T'])
    return LongDNASeq(randstring(letters, n))
end


"
Return a set of `k`-mers based on a string `seq`.
"
function getkmers(seq::String, k::Int = 4)
    return Set([seq[i:i+k-1] for i in 1:(length(seq)-k+1)])
end

"
Return a set of `k`-grams based on a string `seq`.
These will be k non-overlapping consecutive slices of `seq`.
"
function getkgrams(seq::String, k::Int = 5)
    len = Int(floor(length(seq)/k))
    kgrams = Set([seq[1:len], seq[len * (k - 1) + 1:end]])
    for i in 2:(k-1)
        push!(kgrams, seq[len * (i - 1) + 1:len*i])
    end
    return kgrams
end


"
Return pidgeon hole principle minimum required
k-mer size that is required for two strings of
size `len` to be aligned within distance of `d`.
"
function minkmersize(len::Int = 20, d::Int = 4)
    return Int(floor(len / d + 1))
end


"
Get file extension from the string path `s`.
"
function extension(s::String)
    extension = match(r"\.[A-Za-z0-9]+$", s)
    if extension !== nothing
        return extension.match
    else
        return ""
    end
end


"
Returns smallest possible Unsigned type that can contain
given `max_value`.
"
function smallestutype(max_value::Unsigned)
    if typemax(UInt8) >= max_value
        return UInt8
    elseif typemax(UInt16) >= max_value
        return UInt16
    elseif typemax(UInt32) >= max_value
        return UInt32
    elseif typemax(UInt64) >= max_value
        return UInt64
    elseif typemax(UInt128) >= max_value
        return UInt128
    else
        throw("Too big unsigned value to fit in our types.")
    end
end

####### OLD functions
"
Write vector to file in binary format. This is being
serialized which means it can be read back only by the
same version of julia. Will remove file if it exists
before writing.
"
function file_write(write_path::String, vec::Vector)
    # make sure to delete content of write path first
    rm(write_path, force = true)
    io = open(write_path, "w")
    s = Serializer(io)
    for i in vec
        serialize(s, i)
    end
    close(io)
    return nothing
end

"
Read serialized vector from binary file. It can be read
back only by the same version of julia as when used during saving.
"
function file_read(read_path::String)
    x = Vector()
    io = open(read_path, "r")
    s = Serializer(io)
    while !eof(io)
        push!(x, deserialize(s))
    end
    close(io)
    return x
end

"
Append value to file in binary format. This is being
serialized which means it can be read back only by the
same version of julia.
"
function file_add(write_path::String, value)
    io = open(write_path, "a")
    s = Serializer(io)
    serialize(s, value)
    close(io)
end

"
    Will try to find value in `x` that will allow for almost equal
    split of values into buckets.
"
function balance(x::Vector{Int})
    if isempty(x)
        return nothing
    end
    uniq = unique(x)
    sort!(uniq)
    counts = [count(y -> y == i, x) for i in uniq]
    balance = argmin(abs.([sum(counts[1:i]) - sum(counts[i:end]) for i = 1:length(counts)]))
    return uniq[argmin(abs.(uniq .- uniq[balance]))]
end

"
Provides with the bucket path for guides and distances to parent.
"
function bucket_path(dir::String, idx::Int)
    gp = joinpath(dir, string("bucket_", idx, "_g.bin"))
    dp = joinpath(dir, string("bucket_", idx, "_d.bin"))
    return gp, dp
end
