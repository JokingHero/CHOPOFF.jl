# overwrites previous!
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

function load(destination::String)
    io = open(destination, "r")
    s = Serializer(io)
    object = deserialize(s)
    close(io)
    return object
end
