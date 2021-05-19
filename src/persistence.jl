"
Uses julia serializer to save the data to binary format.
Read more at https://docs.julialang.org/en/v1/stdlib/Serialization/
Notice that:

1. This function will overwrite `destination`! 
2. This serialization is dependent on julia build! This means
   files can't be reloaded across different julia builds.
"
function save(
    object::Any,
    destination::String)
    rm(destination, force = true)
    io = open(destination, "w")
    s = Serializer(io)
    serialize(s, object)
    close(io)
    return nothing
end


"
Load file saved with `save` function. This can not load 
properly files saved in other julia builds.
"
function load(destination::String)
    io = open(destination, "r")
    s = Serializer(io)
    object = deserialize(s)
    close(io)
    return object
end
