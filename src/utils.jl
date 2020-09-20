"
Randomize sequence of length `n` from `letters`.
"
function getSeq(n = 20, letters = ['A', 'C', 'G', 'T'])
    return LongDNASeq(randstring(letters, n))
end

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
