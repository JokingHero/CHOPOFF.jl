# Indexed reference reading and global chunk scheduling.

struct PrefixHashScanFASTARecords{R, C, L}
    reader::R
    chrom::C
    lengths::L
end

function prefix_hash_scan_dbinfo(filepath::String, motif::Motif)
    if !is_fasta(filepath)
        index = read_prefix_hash_scan_twobit_index(filepath)
        maxchromlen = isempty(index.lengths) ? 0 : maximum(index.lengths)
        gi = GenomeInfo(
            now(Dates.UTC),
            filepath,
            UInt32(0),
            index.names,
            smallestutype(unsigned(length(index.names))),
            smallestutype(unsigned(maxchromlen)),
            false,
        )
        return DBInfo("prefixHashScan", now(Dates.UTC), gi, "", motif), index.lengths
    end

    fai_path = filepath * ".fai"
    isfile(fai_path) || error("FASTA index not found: " * fai_path)
    index = FASTA.Index(fai_path)
    chrom = Vector{String}(undef, length(index.names))
    for (name, idx) in index.names
        chrom[idx] = name
    end
    maxchromlen = isempty(index.lengths) ? 0 : maximum(index.lengths)
    gi = GenomeInfo(
        now(Dates.UTC),
        filepath,
        UInt32(0),
        chrom,
        smallestutype(unsigned(length(chrom))),
        smallestutype(unsigned(maxchromlen)),
        true,
    )
    return DBInfo("prefixHashScan", now(Dates.UTC), gi, "", motif), index.lengths
end

struct PrefixHashScanIndexedRecords{R, C}
    reader::R
    chrom::C
end

function Base.iterate(records::PrefixHashScanFASTARecords)
    load_start = time_ns()
    next = iterate(records.reader)
    load_ns = time_ns() - load_start
    next === nothing && return nothing
    record, reader_state = next
    isempty(records.chrom) && error("FASTA contains records absent from GenomeInfo.")
    chrom_name = records.chrom[1]
    FASTA.identifier(record) == chrom_name || error("FASTA record order differs from FAI.")
    FASTX.seqsize(record) == records.lengths[1] ||
        error("FASTA record length differs from FAI.")
    return (chrom_name, record, load_ns), (reader_state, 2)
end

function Base.iterate(records::PrefixHashScanFASTARecords, state)
    reader_state, chrom_idx = state
    load_start = time_ns()
    next = iterate(records.reader, reader_state)
    load_ns = time_ns() - load_start
    if next === nothing
        chrom_idx == length(records.chrom) + 1 ||
            error("FASTA contains fewer records than GenomeInfo.")
        return nothing
    end
    chrom_idx <= length(records.chrom) ||
        error("FASTA contains more records than GenomeInfo.")
    record, next_reader_state = next
    chrom_name = records.chrom[chrom_idx]
    FASTA.identifier(record) == chrom_name || error("FASTA record order differs from FAI.")
    FASTX.seqsize(record) == records.lengths[chrom_idx] ||
        error("FASTA record length differs from FAI.")
    return (chrom_name, record, load_ns), (next_reader_state, chrom_idx + 1)
end

function Base.iterate(records::PrefixHashScanIndexedRecords, chrom_idx::Int = 1)
    chrom_idx > length(records.chrom) && return nothing
    chrom_name = records.chrom[chrom_idx]
    load_start = time_ns()
    record = records.reader[chrom_name]
    load_ns = time_ns() - load_start
    return (chrom_name, record, load_ns), chrom_idx + 1
end

function read_prefix_hash_scan_fasta_range!(
    buffer::Vector{UInt8},
    io,
    index::FASTA.Index,
    chrom_idx::Int,
    first_base::Int,
    last_base::Int)

    line_bases, line_width = FASTA.linebases_width(index, chrom_idx)
    first_line, first_offset = divrem(first_base - 1, line_bases)
    last_line, last_offset = divrem(last_base - 1, line_bases)
    physical_first = index.offsets[chrom_idx] +
        first_line * line_width + first_offset
    physical_last = index.offsets[chrom_idx] +
        last_line * line_width + last_offset
    resize!(buffer, physical_last - physical_first + 1)
    seek(io, physical_first)
    read!(io, buffer)

    logical_length = last_base - first_base + 1
    read_idx = 1
    write_idx = 1
    available = line_bases - first_offset
    newline_width = line_width - line_bases
    while write_idx <= logical_length
        count = min(available, logical_length - write_idx + 1)
        copyto!(buffer, write_idx, buffer, read_idx, count)
        write_idx += count
        read_idx += count
        if write_idx <= logical_length
            read_idx += newline_width
            available = line_bases
        end
    end
    resize!(buffer, logical_length)
    return buffer
end

read_prefix_hash_scan_range!(
    buffer::Vector{UInt8}, io, index::FASTA.Index, args...) =
    read_prefix_hash_scan_fasta_range!(buffer, io, index, args...)

read_prefix_hash_scan_range!(
    buffer::Vector{UInt8}, io, index::PrefixHashScanTwoBitIndex, args...) =
    read_prefix_hash_scan_twobit_range!(buffer, io, index, args...)

function stream_prefix_hash_scan_chunk(
    geometry::PrefixScanGeometry,
    io,
    buffer::Vector{UInt8},
    index,
    work::PrefixHashScanChunkWork,
    chromosome_length::Int,
    query,
    dbi::DBInfo,
    guides_::Vector{LongDNA{4}},
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    scratch_plus_hits::Vector{PrefixHashScanHit},
    scratch_minus_hits::Vector{PrefixHashScanHit},
    lookup_scratch::PrefixHashScanLookupScratch,
    ::Val{M},
    stats::S,
    plus::Vector{PrefixHashScanVerifiedHit} = PrefixHashScanVerifiedHit[],
    minus::Vector{PrefixHashScanVerifiedHit} = PrefixHashScanVerifiedHit[]) where {M, S}

    candidate_last_offset = prefix_scan_candidate_last_offset(geometry)
    read_first = max(1, work.core_first - distance)
    read_last = min(
        chromosome_length, work.core_last + candidate_last_offset + distance)
    read_start = prefix_hash_scan_timer(stats)
    raw = read_prefix_hash_scan_range!(
        buffer, io, index, work.chrom_idx, read_first, read_last)
    if stats !== nothing
        read_ns = time_ns() - read_start
        stats.record_io_ns += read_ns
        stats.chrom_load_ns += read_ns
    end

    local_first = work.core_first - read_first + 1
    local_last = work.core_last - read_first + 1
    global_offset = read_first - 1
    if M === :fused
        scan_verify_prefix_raw_range!(
            geometry,
            plus,
            minus,
            raw,
            query,
            local_first,
            local_last,
            local_first,
            local_last,
            local_first,
            local_last,
            global_offset,
            dbi,
            guides_,
            myers_profiles,
            distance,
            stats,
        )
        if dbi.motif.ambig_max > 0
            empty!(scratch_plus_hits)
            empty!(scratch_minus_hits)
            scan_ambiguous_prefix_hits_range!(
                scratch_plus_hits, scratch_minus_hits, raw, geometry, dbi,
                query, geometry.prefix_bases, local_first, local_last,
                local_first, local_last, local_first, local_last,
                Val(dbi.motif.ambig_max), stats)
            evaluate_prefix_hash_scan_hits!(
                plus, raw, geometry, scratch_plus_hits, global_offset, dbi,
                false, guides_, myers_profiles, distance, stats)
            evaluate_prefix_hash_scan_hits!(
                minus, raw, geometry, scratch_minus_hits, global_offset, dbi,
                true, guides_, myers_profiles, distance, stats)
            sort!(plus; by = hit -> hit.pos, alg = QuickSort)
            sort!(minus; by = hit -> hit.pos, alg = QuickSort)
        end
    else
        if M === :buffered_reuse
            motif_candidates = scan_prefix_hits_raw_range!(
                geometry,
                scratch_plus_hits,
                scratch_minus_hits,
                raw,
                query,
                local_first,
                local_last,
                local_first,
                local_last,
                local_first,
                local_last,
            )
            plus_hits = scratch_plus_hits
            minus_hits = scratch_minus_hits
        elseif M === :bucketed_reuse
            motif_candidates = scan_prefix_hits_raw_range_bucketed!(
                geometry,
                scratch_plus_hits,
                scratch_minus_hits,
                lookup_scratch.plus_candidates,
                lookup_scratch.minus_candidates,
                lookup_scratch.plus_radix,
                lookup_scratch.minus_radix,
                lookup_scratch.radix_counts,
                raw,
                query,
                local_first,
                local_last,
                local_first,
                local_last,
                local_first,
                local_last,
            )
            plus_hits = scratch_plus_hits
            minus_hits = scratch_minus_hits
        else
            plus_hits, minus_hits, motif_candidates =
                scan_prefix_hits_raw_range(
                    geometry,
                    raw,
                    query,
                    local_first,
                    local_last,
                    local_first,
                    local_last,
                    local_first,
                    local_last,
                )
        end
        if dbi.motif.ambig_max > 0
            scan_ambiguous_prefix_hits_range!(
                plus_hits, minus_hits, raw, geometry, dbi, query,
                geometry.prefix_bases, local_first, local_last,
                local_first, local_last, local_first, local_last,
                Val(dbi.motif.ambig_max), stats)
        end
        if stats !== nothing
            stats.motif_candidates += motif_candidates
        end
        evaluate_prefix_hash_scan_hits!(
            plus, raw, geometry, plus_hits, global_offset, dbi, false,
            guides_, myers_profiles, distance, stats)
        evaluate_prefix_hash_scan_hits!(
            minus, raw, geometry, minus_hits, global_offset, dbi, true,
            guides_, myers_profiles, distance, stats)
    end
    return PrefixHashScanChromResult(plus, minus, stats)
end

function stream_prefix_hash_scan_chromosome(
    geometry::PrefixScanGeometry,
    io,
    buffer::Vector{UInt8},
    index,
    chrom_idx::Int,
    chromosome_length::Int,
    query,
    dbi::DBInfo,
    guides_::Vector{LongDNA{4}},
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    chunk_bases::Int,
    scratch_plus_hits::Vector{PrefixHashScanHit},
    scratch_minus_hits::Vector{PrefixHashScanHit},
    lookup_scratch::PrefixHashScanLookupScratch,
    mode::Val{M},
    stats::S) where {M, S}

    plus = PrefixHashScanVerifiedHit[]
    minus = PrefixHashScanVerifiedHit[]
    last_candidate =
        chromosome_length - prefix_scan_candidate_last_offset(geometry)
    for core_first in 1:chunk_bases:last_candidate
        work = PrefixHashScanChunkWork(
            chrom_idx, core_first, min(core_first + chunk_bases - 1, last_candidate))
        stream_prefix_hash_scan_chunk(
            geometry, io, buffer, index, work, chromosome_length, query, dbi, guides_,
            myers_profiles, distance, scratch_plus_hits, scratch_minus_hits,
            lookup_scratch, mode, stats, plus, minus)
    end
    return PrefixHashScanChromResult(plus, minus, stats)
end

function prefix_hash_scan_chunk_work(
    reference_lengths, chunk_bases::Int,
    geometry::PrefixScanGeometry = CAS9_D3_PREFIX_SCAN_GEOMETRY)
    work = PrefixHashScanChunkWork[]
    chrom_chunk_ranges = Vector{UnitRange{Int}}(undef, length(reference_lengths))
    for chrom_idx in eachindex(reference_lengths)
        range_first = length(work) + 1
        last_candidate = reference_lengths[chrom_idx] -
            prefix_scan_candidate_last_offset(geometry)
        for core_first in 1:chunk_bases:last_candidate
            push!(work, PrefixHashScanChunkWork(
                chrom_idx,
                core_first,
                min(core_first + chunk_bases - 1, last_candidate),
            ))
        end
        chrom_chunk_ranges[chrom_idx] = range_first:length(work)
    end
    return work, chrom_chunk_ranges
end

function stream_prefix_hash_scan(
    geometry::PrefixScanGeometry,
    genome_path::String,
    reference_lengths,
    query,
    dbi::DBInfo,
    guides_::Vector{LongDNA{4}},
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    chunk_bases::Int,
    scan_threads::Int,
    mode::Val{M},
    stats::S,
    ::Val{Scheduler} = Val(:chunk)) where {M, S, Scheduler}

    index = is_fasta(genome_path) ?
        FASTA.Index(genome_path * ".fai") :
        read_prefix_hash_scan_twobit_index(genome_path)
    if Scheduler === :chunk
        work, chrom_chunk_ranges = prefix_hash_scan_chunk_work(
            reference_lengths, chunk_bases, geometry)
    elseif Scheduler === :chromosome
        work = collect(eachindex(reference_lengths))
        chrom_chunk_ranges = [idx:idx for idx in eachindex(reference_lengths)]
    else
        error("Unknown prefixHashScan stream scheduler: $Scheduler")
    end
    results = Vector{Union{Nothing, PrefixHashScanChromResult}}(
        nothing, length(work))
    next_work = Threads.Atomic{Int}(1)
    worker_count = min(scan_threads, length(work))
    workers = map(1:worker_count) do _
        Threads.@spawn begin
            io = open(genome_path, "r")
            buffer = UInt8[]
            scratch_plus_hits = PrefixHashScanHit[]
            scratch_minus_hits = PrefixHashScanHit[]
            lookup_scratch = PrefixHashScanLookupScratch()
            try
                while true
                    work_idx = Threads.atomic_add!(next_work, 1)
                    work_idx > length(work) && break
                    worker_stats = prefix_hash_scan_worker_stats(stats)
                    if Scheduler === :chunk
                        item = work[work_idx]
                        results[work_idx] = stream_prefix_hash_scan_chunk(
                            geometry, io, buffer, index, item,
                            reference_lengths[item.chrom_idx], query, dbi,
                            guides_, myers_profiles, distance, scratch_plus_hits,
                            scratch_minus_hits, lookup_scratch, mode, worker_stats)
                    else
                        chrom_idx = work[work_idx]
                        results[work_idx] = stream_prefix_hash_scan_chromosome(
                            geometry, io, buffer, index, chrom_idx,
                            reference_lengths[chrom_idx], query, dbi, guides_,
                            myers_profiles, distance, chunk_bases,
                            scratch_plus_hits, scratch_minus_hits, lookup_scratch,
                            mode, worker_stats)
                    end
                end
            finally
                close(io)
            end
        end
    end
    fetch.(workers)
    return results, chrom_chunk_ranges
end

stream_prefix_hash_scan(genome_path::String, args...) =
    stream_prefix_hash_scan(
        CAS9_D3_PREFIX_SCAN_GEOMETRY, genome_path, args...)
