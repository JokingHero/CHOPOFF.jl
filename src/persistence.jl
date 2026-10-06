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
`with_atomic_output(f::Function, output_file::String)`

Calls `f(temporary)` with a fresh file next to `output_file`, then renames it
over `output_file`. Readers see either the old file or the complete new one;
if `f` throws, `output_file` is left unchanged.
"""
function with_atomic_output(f::Function, output_file::String)
    output_dir = dirname(abspath(output_file))
    mkpath(output_dir)
    temporary, temporary_io = mktemp(output_dir; cleanup = false)
    close(temporary_io)
    try
        f(temporary)
        Base.Filesystem.rename(temporary, output_file; force = true)
    finally
        ispath(temporary) && rm(temporary; force = true)
    end
    return nothing
end


const DETAIL_FIRST_LINE =
    "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n"

"""
`with_detail_parts(f::Function, output_file::String; first_line::String = DETAIL_FIRST_LINE)`

Calls `f(parts_dir)` with a new private directory next to `output_file`, where
search workers write their partial detail files. The parts are then
concatenated in file-name order under `first_line` and published atomically as
`output_file`. The parts directory is removed even when `f` throws, and no
other file in the output directory is read or deleted.
"""
function with_detail_parts(
    f::Function, output_file::String; first_line::String = DETAIL_FIRST_LINE)

    output_dir = dirname(abspath(output_file))
    mkpath(output_dir)
    parts_dir = mktempdir(output_dir; prefix = ".chopoff_parts_", cleanup = false)
    try
        f(parts_dir)
        with_atomic_output(output_file) do temporary
            open(temporary, "w") do io
                write(io, first_line)
                for part in readdir(parts_dir; join = true)
                    for ln in eachline(part)
                        write(io, ln, "\n")
                    end
                end
            end
        end
    finally
        rm(parts_dir; recursive = true, force = true)
    end
    return nothing
end
