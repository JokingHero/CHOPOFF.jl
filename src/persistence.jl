struct SketchDB{T<:Unsigned}
    sketch::CountMinSketch{T}
    kmers::Dict{String, Int}
    motif::Motif
    max_dist::Int
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
