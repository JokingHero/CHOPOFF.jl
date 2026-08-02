# Chunk-oriented UCSC 2bit metadata and range decoding.

const PREFIX_HASH_SCAN_TWOBIT_SIGNATURE = UInt32(0x1a412743)
const PREFIX_HASH_SCAN_TWOBIT_BASES = (
    UInt8('T'), UInt8('C'), UInt8('A'), UInt8('G'))
const PREFIX_HASH_SCAN_TWOBIT_DECODE = ntuple(256) do idx
    byte = UInt8(idx - 1)
    (
        PREFIX_HASH_SCAN_TWOBIT_BASES[Int((byte >> 6) & 0x03) + 1],
        PREFIX_HASH_SCAN_TWOBIT_BASES[Int((byte >> 4) & 0x03) + 1],
        PREFIX_HASH_SCAN_TWOBIT_BASES[Int((byte >> 2) & 0x03) + 1],
        PREFIX_HASH_SCAN_TWOBIT_BASES[Int(byte & 0x03) + 1],
    )
end

struct PrefixHashScanTwoBitIndex
    names::Vector{String}
    lengths::Vector{Int}
    packed_offsets::Vector{Int64}
    nblock_starts::Vector{Vector{UInt32}}
    nblock_sizes::Vector{Vector{UInt32}}
end

@inline function read_prefix_hash_scan_twobit_u32(io, swap::Bool)
    value = read(io, UInt32)
    return swap ? bswap(value) : value
end

function read_prefix_hash_scan_twobit_u32s(io, count::Int, swap::Bool)
    values = Vector{UInt32}(undef, count)
    read!(io, values)
    swap && map!(bswap, values, values)
    return values
end

function read_prefix_hash_scan_twobit_index(filepath::String)
    open(filepath, "r") do io
        signature = read(io, UInt32)
        swap = if signature == PREFIX_HASH_SCAN_TWOBIT_SIGNATURE
            false
        elseif bswap(signature) == PREFIX_HASH_SCAN_TWOBIT_SIGNATURE
            true
        else
            error("Invalid 2bit file signature.")
        end
        read_prefix_hash_scan_twobit_u32(io, swap) == 0 ||
            error("Unsupported 2bit file version.")
        sequence_count = Int(read_prefix_hash_scan_twobit_u32(io, swap))
        read_prefix_hash_scan_twobit_u32(io, swap) == 0 ||
            error("Invalid 2bit file header.")

        names = Vector{String}(undef, sequence_count)
        record_offsets = Vector{UInt32}(undef, sequence_count)
        for idx in 1:sequence_count
            name_size = Int(read(io, UInt8))
            names[idx] = String(read(io, name_size))
            record_offsets[idx] = read_prefix_hash_scan_twobit_u32(io, swap)
        end

        lengths = Vector{Int}(undef, sequence_count)
        packed_offsets = Vector{Int64}(undef, sequence_count)
        nblock_starts = Vector{Vector{UInt32}}(undef, sequence_count)
        nblock_sizes = Vector{Vector{UInt32}}(undef, sequence_count)
        for idx in 1:sequence_count
            seek(io, record_offsets[idx])
            lengths[idx] = Int(read_prefix_hash_scan_twobit_u32(io, swap))
            nblock_count = Int(read_prefix_hash_scan_twobit_u32(io, swap))
            nblock_starts[idx] = read_prefix_hash_scan_twobit_u32s(
                io, nblock_count, swap)
            nblock_sizes[idx] = read_prefix_hash_scan_twobit_u32s(
                io, nblock_count, swap)
            mask_count = Int(read_prefix_hash_scan_twobit_u32(io, swap))
            skip(io, 8 * mask_count)
            read_prefix_hash_scan_twobit_u32(io, swap) == 0 ||
                error("Invalid 2bit sequence record.")
            packed_offsets[idx] = position(io)
        end
        return PrefixHashScanTwoBitIndex(
            names, lengths, packed_offsets, nblock_starts, nblock_sizes)
    end
end

function read_prefix_hash_scan_twobit_range!(
    buffer::Vector{UInt8},
    io,
    index::PrefixHashScanTwoBitIndex,
    chrom_idx::Int,
    first_base::Int,
    last_base::Int)

    1 <= first_base <= last_base <= index.lengths[chrom_idx] ||
        throw(BoundsError(index.lengths[chrom_idx], first_base:last_base))
    first_byte = div(first_base - 1, 4)
    last_byte = div(last_base - 1, 4)
    byte_count = last_byte - first_byte + 1
    resize!(buffer, byte_count)
    seek(io, index.packed_offsets[chrom_idx] + first_byte)
    read!(io, buffer)

    resize!(buffer, 4 * byte_count)
    @inbounds for byte_idx in byte_count:-1:1
        decoded = PREFIX_HASH_SCAN_TWOBIT_DECODE[Int(buffer[byte_idx]) + 1]
        output_idx = 4 * byte_idx - 3
        buffer[output_idx] = decoded[1]
        buffer[output_idx + 1] = decoded[2]
        buffer[output_idx + 2] = decoded[3]
        buffer[output_idx + 3] = decoded[4]
    end

    first_offset = mod(first_base - 1, 4)
    output_length = last_base - first_base + 1
    copyto!(buffer, 1, buffer, first_offset + 1, output_length)
    resize!(buffer, output_length)

    query_first = first_base - 1
    query_last = last_base
    starts = index.nblock_starts[chrom_idx]
    sizes = index.nblock_sizes[chrom_idx]
    block_idx = searchsortedlast(starts, UInt32(query_first))
    if block_idx == 0
        block_idx = 1
    elseif Int(starts[block_idx]) + Int(sizes[block_idx]) <= query_first
        block_idx += 1
    end
    @inbounds while block_idx <= length(starts)
        block_first = Int(starts[block_idx])
        block_first >= query_last && break
        block_last = block_first + Int(sizes[block_idx])
        if block_last > query_first
            local_first = max(block_first, query_first) - query_first + 1
            local_last = min(block_last, query_last) - query_first
            fill!(@view(buffer[local_first:local_last]), UInt8('N'))
        end
        block_idx += 1
    end
    return buffer
end
