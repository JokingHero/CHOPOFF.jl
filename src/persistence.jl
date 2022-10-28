"""
`save(object::Any, destination::String)`

Uses julia serializer to save the data to binary format.
Read more about [serialization](https://docs.julialang.org/en/v1/stdlib/Serialization/).
Notice that:
    
1. This function will overwrite `destination`! 
2. This serialization is dependent on julia build! This means files can fail to work when reloaded across different julia builds.
"""
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


"""
`load(destination::String)`

Load file saved with `save` function. This **may not** load properly files saved in other julia builds.
"""
function load(destination::String)
    io = open(destination, "r")
    s = Serializer(io)
    object = deserialize(s)
    close(io)
    return object
end


"""
`cleanup_detail(detail::String)`

Merge multiple detail files into one final file. 
This assumes there are no other files in the folder that start with "detail_".
After the merge files that start with "detail_" will be deleted!
"""
function cleanup_detail(detail::String)
    open(detail, "w") do detail_file
        write(detail_file, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")
        for ch in filter(x -> occursin("detail_", x), readdir(dirname(detail)))
            ch = joinpath(dirname(detail), ch)
            open(ch, "r") do ch_file
                for ln in eachline(ch_file)
                    write(detail_file, ln * "\n")
                end
            end
            rm(ch)
        end
    end
    return 
end