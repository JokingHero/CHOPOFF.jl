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

# First and last non-N base of a chromosome, 1-based, or `nothing` when the
# whole record is unknown sequence.
#
# The whole-chromosome kernels (cas9_prefix_scan_bounds_raw and friends) get the
# entire record in one buffer and trim its terminal N runs, which bounds how far
# a candidate window may reach into unassembled sequence. The chunked scanner
# never sees a whole chromosome, and chunks are claimed out of order, so it
# cannot derive the same bounds locally: a chunk starting inside an *internal* N
# run must not be trimmed, and "was everything before me N?" is not knowable
# from the chunk alone. So the terminal runs are located once per chromosome and
# handed to the workers.
#
# Reads outward in windows rather than pulling the record in whole; terminal N
# runs are short relative to a chromosome, so this normally costs two reads.
function prefix_hash_scan_known_bounds(
    buffer::Vector{UInt8}, io, index, chrom_idx::Int,
    chromosome_length::Int; probe_bases::Int = 1 << 16)

    is_unknown(base::UInt8) = base == UInt8('N') || base == UInt8('n')

    seq_start = 0
    pos = 1
    while pos <= chromosome_length
        stop = min(pos + probe_bases - 1, chromosome_length)
        raw = read_prefix_hash_scan_range!(buffer, io, index, chrom_idx, pos, stop)
        offset = findfirst(!is_unknown, raw)
        if offset !== nothing
            seq_start = pos + offset - 1
            break
        end
        pos = stop + 1
    end
    seq_start == 0 && return nothing

    seq_stop = 0
    pos = chromosome_length
    while pos >= seq_start
        start = max(pos - probe_bases + 1, seq_start)
        raw = read_prefix_hash_scan_range!(buffer, io, index, chrom_idx, start, pos)
        offset = findlast(!is_unknown, raw)
        if offset !== nothing
            seq_stop = start + offset - 1
            break
        end
        pos = start - 1
    end
    return seq_start, seq_stop
end

function prefix_hash_scan_known_bounds_by_chromosome(
    genome_path::String, index, reference_lengths)

    bounds = Vector{Union{Nothing, Tuple{Int, Int}}}(
        nothing, length(reference_lengths))
    buffer = UInt8[]
    open(genome_path, "r") do io
        for chrom_idx in eachindex(reference_lengths)
            bounds[chrom_idx] = prefix_hash_scan_known_bounds(
                buffer, io, index, chrom_idx, reference_lengths[chrom_idx])
        end
    end
    return bounds
end

# Buffer-local scan bounds for one chunk, matching what the whole-chromosome
# kernels compute: a plus-strand window may begin at most `motif.distance` bases
# before the first known base, and must end by the last known base; the mirror
# applies on the minus strand. Returns `nothing` when nothing in this chunk is
# scannable.
function prefix_hash_scan_chunk_bounds(
    geometry::PrefixScanGeometry, dbi::DBInfo,
    known::Union{Nothing, Tuple{Int, Int}},
    read_first::Int, raw_length::Int, local_first::Int, local_last::Int)

    known === nothing && return nothing
    candidate_last_offset = prefix_scan_candidate_last_offset(geometry)
    seq_start = known[1] - read_first + 1
    seq_stop = known[2] - read_first + 1
    motif_distance = dbi.motif.distance

    plus_first = max(seq_start - motif_distance, local_first)
    plus_last = min(seq_stop - candidate_last_offset, local_last)
    minus_first = max(seq_start, local_first)
    minus_last = min(
        min(seq_stop + motif_distance, raw_length) - candidate_last_offset,
        local_last)

    firsts = Int[]
    lasts = Int[]
    plus_first <= plus_last && (push!(firsts, plus_first); push!(lasts, plus_last))
    minus_first <= minus_last && (push!(firsts, minus_first); push!(lasts, minus_last))
    isempty(firsts) && return nothing
    return (minimum(firsts), maximum(lasts),
        plus_first, plus_last, minus_first, minus_last)
end

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
    simd_backend::Val,
    stats::S,
    plus::Vector{PrefixHashScanVerifiedHit} = PrefixHashScanVerifiedHit[],
    minus::Vector{PrefixHashScanVerifiedHit} = PrefixHashScanVerifiedHit[],
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    known_bounds::Union{Nothing, Tuple{Int, Int}} = nothing,
    ) where {M, S}

    chunk_counts = early_stop_state === nothing ? nothing :
        zeros(Int, length(guides_), distance + 1)

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
    bounds = prefix_hash_scan_chunk_bounds(
        geometry, dbi, known_bounds, read_first, length(raw),
        local_first, local_last)
    bounds === nothing && return PrefixHashScanChromResult(plus, minus, stats)
    candidate_first, candidate_last, plus_first, plus_last,
        minus_first, minus_last = bounds
    if M === :fused
        scan_verify_prefix_raw_range!(
            geometry,
            plus,
            minus,
            raw,
            query,
            candidate_first,
            candidate_last,
            plus_first,
            plus_last,
            minus_first,
            minus_last,
            global_offset,
            dbi,
            guides_,
            myers_profiles,
            distance,
            stats,
            simd_backend,
        )
        if dbi.motif.ambig_max > 0
            empty!(scratch_plus_hits)
            empty!(scratch_minus_hits)
            scan_ambiguous_prefix_hits_range!(
                scratch_plus_hits, scratch_minus_hits, raw, geometry, dbi,
                query, geometry.prefix_bases, candidate_first, candidate_last,
                plus_first, plus_last, minus_first, minus_last,
                Val(dbi.motif.ambig_max), stats, simd_backend)
            evaluate_prefix_hash_scan_hits!(
                plus, raw, geometry, scratch_plus_hits, global_offset, dbi,
                false, guides_, myers_profiles, distance, stats,
                early_stop_state, chunk_counts)
            evaluate_prefix_hash_scan_hits!(
                minus, raw, geometry, scratch_minus_hits, global_offset, dbi,
                true, guides_, myers_profiles, distance, stats,
                early_stop_state, chunk_counts)
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
                candidate_first,
                candidate_last,
                plus_first,
                plus_last,
                minus_first,
                minus_last,
                simd_backend,
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
                candidate_first,
                candidate_last,
                plus_first,
                plus_last,
                minus_first,
                minus_last,
                simd_backend,
            )
            plus_hits = scratch_plus_hits
            minus_hits = scratch_minus_hits
        else
            plus_hits, minus_hits, motif_candidates =
                scan_prefix_hits_raw_range(
                    geometry,
                    raw,
                    query,
                    candidate_first,
                    candidate_last,
                    plus_first,
                    plus_last,
                    minus_first,
                    minus_last,
                    simd_backend,
                )
        end
        if dbi.motif.ambig_max > 0
            scan_ambiguous_prefix_hits_range!(
                plus_hits, minus_hits, raw, geometry, dbi, query,
                geometry.prefix_bases, candidate_first, candidate_last,
                plus_first, plus_last, minus_first, minus_last,
                Val(dbi.motif.ambig_max), stats, simd_backend)
        end
        if stats !== nothing
            stats.motif_candidates += motif_candidates
        end
        evaluate_prefix_hash_scan_hits!(
            plus, raw, geometry, plus_hits, global_offset, dbi, false,
            guides_, myers_profiles, distance, stats,
            early_stop_state, chunk_counts)
        evaluate_prefix_hash_scan_hits!(
            minus, raw, geometry, minus_hits, global_offset, dbi, true,
            guides_, myers_profiles, distance, stats,
            early_stop_state, chunk_counts)
    end
    chunk_counts === nothing ||
        merge_prefix_hash_scan_chunk_counts!(early_stop_state, chunk_counts)
    return PrefixHashScanChromResult(plus, minus, stats)
end

function stream_prefix_hash_scan_count_chunk(
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
    simd_backend::Val,
    stats::S,
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    known_bounds::Union{Nothing, Tuple{Int, Int}} = nothing,
    ) where {M, S}

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
    bounds = prefix_hash_scan_chunk_bounds(
        geometry, dbi, known_bounds, read_first, length(raw),
        local_first, local_last)
    bounds === nothing &&
        return PrefixHashScanCountResult(
            zeros(Int, length(guides_), distance + 1), stats)
    candidate_first, candidate_last, plus_first, plus_last,
        minus_first, minus_last = bounds
    if M === :buffered_reuse
        motif_candidates = scan_prefix_hits_raw_range!(
            geometry, scratch_plus_hits, scratch_minus_hits, raw, query,
            candidate_first, candidate_last, plus_first, plus_last,
            minus_first, minus_last, simd_backend)
    elseif M === :bucketed_reuse
        motif_candidates = scan_prefix_hits_raw_range_bucketed!(
            geometry, scratch_plus_hits, scratch_minus_hits,
            lookup_scratch.plus_candidates, lookup_scratch.minus_candidates,
            lookup_scratch.plus_radix, lookup_scratch.minus_radix,
            lookup_scratch.radix_counts, raw, query,
            candidate_first, candidate_last, plus_first, plus_last,
            minus_first, minus_last, simd_backend)
    else
        error("Unknown prefixHashScan count mode: $M")
    end
    if dbi.motif.ambig_max > 0
        scan_ambiguous_prefix_hits_range!(
            scratch_plus_hits, scratch_minus_hits, raw, geometry, dbi, query,
            geometry.prefix_bases, candidate_first, candidate_last,
            plus_first, plus_last, minus_first, minus_last,
            Val(dbi.motif.ambig_max), stats, simd_backend)
    end
    stats === nothing || (stats.motif_candidates += motif_candidates)

    counts = zeros(Int, length(guides_), distance + 1)
    evaluate_prefix_hash_scan_count_hits!(
        counts, raw, geometry, scratch_plus_hits, false,
        myers_profiles, distance, stats, early_stop_state, counts)
    evaluate_prefix_hash_scan_count_hits!(
        counts, raw, geometry, scratch_minus_hits, true,
        myers_profiles, distance, stats, early_stop_state, counts)
    early_stop_state === nothing ||
        merge_prefix_hash_scan_chunk_counts!(early_stop_state, counts)
    return PrefixHashScanCountResult(counts, stats)
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
    simd_backend::Val,
    stats::S,
    known_bounds::Union{Nothing, Tuple{Int, Int}} = nothing) where {M, S}

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
            lookup_scratch, mode, simd_backend, stats, plus, minus,
            nothing, known_bounds)
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
    ::Val{Scheduler} = Val(:chunk),
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    ; simd_backend::Val = default_prefix_hash_scan_simd_backend()) where {M, S, Scheduler}

    index = is_fasta(genome_path) ?
        FASTA.Index(genome_path * ".fai") :
        read_prefix_hash_scan_twobit_index(genome_path)
    known_bounds = prefix_hash_scan_known_bounds_by_chromosome(
        genome_path, index, reference_lengths)
    if Scheduler === :chunk
        work, chrom_chunk_ranges = prefix_hash_scan_chunk_work(
            reference_lengths, chunk_bases, geometry)
    elseif Scheduler === :chromosome
        work = collect(eachindex(reference_lengths))
        chrom_chunk_ranges = [idx:idx for idx in eachindex(reference_lengths)]
    else
        error("Unknown prefixHashScan stream scheduler: $Scheduler")
    end
    if early_stop_state !== nothing
        early_stop_state.work_items_total = length(work)
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
                    prefix_hash_scan_active_mask(early_stop_state) == 0 && break
                    work_idx = Threads.atomic_add!(next_work, 1)
                    work_idx > length(work) && break
                    early_stop_state === nothing ||
                        Threads.atomic_add!(early_stop_state.work_items_claimed, 1)
                    worker_stats = prefix_hash_scan_worker_stats(stats)
                    if Scheduler === :chunk
                        item = work[work_idx]
                        results[work_idx] = stream_prefix_hash_scan_chunk(
                            geometry, io, buffer, index, item,
                            reference_lengths[item.chrom_idx], query, dbi,
                            guides_, myers_profiles, distance, scratch_plus_hits,
                            scratch_minus_hits, lookup_scratch, mode, simd_backend,
                            worker_stats,
                            PrefixHashScanVerifiedHit[], PrefixHashScanVerifiedHit[],
                            early_stop_state, known_bounds[item.chrom_idx])
                    else
                        chrom_idx = work[work_idx]
                        results[work_idx] = stream_prefix_hash_scan_chromosome(
                            geometry, io, buffer, index, chrom_idx,
                            reference_lengths[chrom_idx], query, dbi, guides_,
                            myers_profiles, distance, chunk_bases,
                            scratch_plus_hits, scratch_minus_hits, lookup_scratch,
                            mode, simd_backend, worker_stats,
                            known_bounds[chrom_idx])
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

function stream_prefix_hash_scan_counts(
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
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    ; simd_backend::Val = default_prefix_hash_scan_simd_backend()) where {M, S}

    index = is_fasta(genome_path) ?
        FASTA.Index(genome_path * ".fai") :
        read_prefix_hash_scan_twobit_index(genome_path)
    known_bounds = prefix_hash_scan_known_bounds_by_chromosome(
        genome_path, index, reference_lengths)
    work, _ = prefix_hash_scan_chunk_work(
        reference_lengths, chunk_bases, geometry)
    if early_stop_state !== nothing
        early_stop_state.work_items_total = length(work)
    end
    results = Vector{Union{Nothing, PrefixHashScanCountResult}}(
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
                    prefix_hash_scan_active_mask(early_stop_state) == 0 && break
                    work_idx = Threads.atomic_add!(next_work, 1)
                    work_idx > length(work) && break
                    early_stop_state === nothing ||
                        Threads.atomic_add!(early_stop_state.work_items_claimed, 1)
                    item = work[work_idx]
                    worker_stats = prefix_hash_scan_worker_stats(stats)
                    results[work_idx] = stream_prefix_hash_scan_count_chunk(
                        geometry, io, buffer, index, item,
                        reference_lengths[item.chrom_idx], query, dbi, guides_,
                        myers_profiles, distance, scratch_plus_hits,
                        scratch_minus_hits, lookup_scratch, mode, simd_backend,
                        worker_stats,
                        early_stop_state, known_bounds[item.chrom_idx])
                end
            finally
                close(io)
            end
        end
    end
    fetch.(workers)
    return results
end

stream_prefix_hash_scan(genome_path::String, args...) =
    stream_prefix_hash_scan(
        CAS9_D3_PREFIX_SCAN_GEOMETRY, genome_path, args...)
