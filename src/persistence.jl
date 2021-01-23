struct SketchDB{T<:Unsigned}
    sketches::Vector{CountMinSketch{T}}
    motif::Motif
end

# overwrites previous!
function saveDB(
    db::SketchDB,
    destination::String)
    rm(destination, force = true)
    io = open(destination, "w")
    s = Serializer(io)
    serialize(s, db)
    close(io)
    return nothing
end

function loadDB(destination::String)
    io = open(destination, "r")
    s = Serializer(io)
    db::SketchDB = deserialize(s)
    close(io)
    return db
end
