# Multithreaded VCF decompression. BGZF (bgzip output, required for tabix) is a
# series of gzip members of at most 64 KiB, each header holding its compressed
# size, so blocks inflate independently while the caller parses earlier ones.
# Plain gzip falls back to the single-threaded GzipDecompressorStream.

# Blocks per thread per batch; about 1 MiB of text per thread.
const VARIANT_BGZF_BATCH_BLOCKS = 16
const VARIANT_BGZF_HEADER = 18

# Total block size from a BGZF header (gzip, deflate, FEXTRA, BC subfield first,
# as bgzip and htslib write it), or nothing for other input.
function variant_bgzf_block_size(header)
    length(header) == VARIANT_BGZF_HEADER && header[1] == 0x1f && header[2] == 0x8b &&
        header[3] == 0x08 && header[4] & 0x04 != 0 &&
        header[13] == UInt8('B') && header[14] == UInt8('C') || return nothing
    return (Int(header[17]) | (Int(header[18]) << 8)) + 1
end

# Inflates blocks[part] into their slots of text. CodecZlib exposes no
# public call that writes into an existing buffer, so this uses its ZStream.
function variant_bgzf_inflate!(text, blocks, sizes, ends, part, path)
    z = CodecZlib.ZStream()
    CodecZlib.inflate_init!(z, 31) == CodecZlib.Z_OK || error("zlib initialization failed")
    try
        for i in part
            first = i == 1 ? 0 : ends[i - 1]
            CodecZlib.inflate_reset!(z)
            code = GC.@preserve text blocks begin
                z.next_in, z.avail_in = pointer(blocks[i]), sizes[i]
                z.next_out, z.avail_out = pointer(text, first + 1), ends[i] - first
                CodecZlib.inflate!(z, CodecZlib.Z_FINISH)
            end
            code == CodecZlib.Z_STREAM_END && z.avail_out == 0 && z.avail_in == 0 ||
                error("Corrupt BGZF block in $path")
        end
    finally
        CodecZlib.inflate_end!(z)
    end
end

# Reads up to n blocks into the reused blocks and sizes; returns the count read.
function variant_bgzf_read_blocks!(blocks, sizes, io, n, path)
    k = 0
    while k < n && !eof(io)
        k += 1
        k > length(blocks) && (push!(blocks, UInt8[]); push!(sizes, 0))
        block = blocks[k]
        readbytes!(io, block, VARIANT_BGZF_HEADER) == VARIANT_BGZF_HEADER ||
            error("Truncated BGZF block in $path")
        size = variant_bgzf_block_size(@view block[1:VARIANT_BGZF_HEADER])
        size === nothing && error("Malformed BGZF block in $path")
        skip(io, -VARIANT_BGZF_HEADER)
        readbytes!(io, block, size) == size || error("Truncated BGZF block in $path")
        sizes[k] = size
    end
    return k
end

mutable struct VariantBGZFStream <: IO
    io::IOStream
    chunks::Channel{Vector{UInt8}}
    free::Channel{Vector{UInt8}}    # consumed text buffers for the producer to reuse
    buffer::Vector{UInt8}
    position::Int
end

function open_variant_vcf(path)
    io = open(path, "r")
    endswith(lowercase(path), ".gz") || return io
    header = read(io, VARIANT_BGZF_HEADER)
    seekstart(io)
    variant_bgzf_block_size(header) === nothing && return GzipDecompressorStream(io)
    free = Channel{Vector{UInt8}}(Inf)
    # The bound channel rethrows a producer error on the next take.
    chunks = Channel{Vector{UInt8}}(2; spawn = true) do chunks
        batch = VARIANT_BGZF_BATCH_BLOCKS * Threads.nthreads()
        blocks, sizes = Vector{UInt8}[], Int[]
        while !eof(io)
            n = variant_bgzf_read_blocks!(blocks, sizes, io, batch, path)
            # ISIZE trailers give each block's text offset in the batch.
            ends = cumsum([Int(ltoh(reinterpret(UInt32, blocks[i][(sizes[i] - 3):sizes[i]])[1])) for i in 1:n])
            text = resize!(isready(free) ? take!(free) : UInt8[], ends[end])
            @sync for part in Iterators.partition(1:n, cld(n, Threads.nthreads()))
                Threads.@spawn variant_bgzf_inflate!(text, blocks, sizes, ends, part, path)
            end
            put!(chunks, text)
        end
    end
    return VariantBGZFStream(io, chunks, free, UInt8[], 1)
end

function Base.eof(s::VariantBGZFStream)
    while s.position > length(s.buffer)
        next = iterate(s.chunks)
        next === nothing && return true
        put!(s.free, s.buffer)
        s.buffer, s.position = next[1], 1
    end
    return false
end

function Base.readline(s::VariantBGZFStream; keep::Bool = false)
    line = UInt8[]
    while !eof(s)
        j = findnext(==(0x0a), s.buffer, s.position)
        stop = something(j, length(s.buffer))
        isempty(line) ? (line = s.buffer[s.position:stop]) : append!(line, @view s.buffer[s.position:stop])
        s.position = stop + 1
        j === nothing || break
    end
    if !keep && !isempty(line) && line[end] == 0x0a
        pop!(line)
        !isempty(line) && line[end] == 0x0d && pop!(line)
    end
    return String(line)
end

function Base.close(s::VariantBGZFStream)
    close(s.chunks)
    close(s.free)
    close(s.io)
end
