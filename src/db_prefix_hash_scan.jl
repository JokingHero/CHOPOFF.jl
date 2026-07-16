using SIMD

struct PrefixScanGeometry
    guide_bases::Int
    pam_bases::Int
    prefix_bases::Int
    distance::Int
end

const CAS9_D3_PREFIX_SCAN_GEOMETRY = PrefixScanGeometry(20, 3, 16, 3)

@inline prefix_scan_candidate_bases(geometry::PrefixScanGeometry) =
    geometry.guide_bases + geometry.pam_bases

@inline prefix_scan_candidate_last_offset(geometry::PrefixScanGeometry) =
    prefix_scan_candidate_bases(geometry) - 1

mutable struct PrefixHashScanStats
    motif_candidates::Int
    ambiguous_prefixes::Int
    prefix_hits::Int
    guide_pairs::Int
    alignment_calls::Int
    distance_calls::Int
    traceback_calls::Int
    emitted_rows::Int
    path_rows::Int
    query_hashes::Int
    bruteforce_guide_pairs::Int
    metadata_ns::UInt64
    query_build_ns::UInt64
    path_load_ns::UInt64
    query_hash_ns::UInt64
    query_format_ns::UInt64
    query_fold_ns::UInt64
    query_dedup_ns::UInt64
    query_insert_ns::UInt64
    query_lookup_ns::UInt64
    chrom_load_ns::UInt64
    record_io_ns::UInt64
    sequence_convert_ns::UInt64
    findguides_ns::UInt64
    candidate_prefix_ns::UInt64
    candidate_hash_ns::UInt64
    candidate_materialize_ns::UInt64
    align_ns::UInt64
    emit_ns::UInt64
    scan_ns::UInt64
    verify_ns::UInt64
    path_source::Symbol
    query_variant::Symbol
    scan_backend::Symbol
end

PrefixHashScanStats() = PrefixHashScanStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, :none, :none, :none)

@inline prefix_hash_scan_timer(::Nothing) = UInt64(0)
@inline prefix_hash_scan_timer(::PrefixHashScanStats) = time_ns()
@inline prefix_hash_scan_worker_stats(::Nothing) = nothing
@inline prefix_hash_scan_worker_stats(::PrefixHashScanStats) = PrefixHashScanStats()

function reset!(stats::PrefixHashScanStats)
    stats.motif_candidates = 0
    stats.ambiguous_prefixes = 0
    stats.prefix_hits = 0
    stats.guide_pairs = 0
    stats.alignment_calls = 0
    stats.distance_calls = 0
    stats.traceback_calls = 0
    stats.emitted_rows = 0
    stats.path_rows = 0
    stats.query_hashes = 0
    stats.bruteforce_guide_pairs = 0
    stats.metadata_ns = 0
    stats.query_build_ns = 0
    stats.path_load_ns = 0
    stats.query_hash_ns = 0
    stats.query_format_ns = 0
    stats.query_fold_ns = 0
    stats.query_dedup_ns = 0
    stats.query_insert_ns = 0
    stats.query_lookup_ns = 0
    stats.chrom_load_ns = 0
    stats.record_io_ns = 0
    stats.sequence_convert_ns = 0
    stats.findguides_ns = 0
    stats.candidate_prefix_ns = 0
    stats.candidate_hash_ns = 0
    stats.candidate_materialize_ns = 0
    stats.align_ns = 0
    stats.emit_ns = 0
    stats.scan_ns = 0
    stats.verify_ns = 0
    stats.path_source = :none
    stats.query_variant = :none
    stats.scan_backend = :none
    return stats
end

struct PrefixHashScanBitmaskQuery{H <: Unsigned}
    masks::Dict{H, UInt64}
    nguides::Int
end

struct PrefixHashScanDirectory
    offsets::Vector{UInt32}
    suffixes::Vector{UInt16}
    masks::Vector{UInt64}
    hash_len::UInt8
    bucket_bases::UInt8
end

struct PrefixHashScanPrefilteredDirectory{D}
    directory::D
    presence::Vector{UInt64}
    prefix_bits::UInt8
end

struct PrefixHashScanHit
    start::Int
    mask::UInt64
end

struct PrefixHashScanLookupScratch
    plus_candidates::Vector{UInt64}
    minus_candidates::Vector{UInt64}
    plus_radix::Vector{UInt64}
    minus_radix::Vector{UInt64}
    radix_counts::Vector{Int}
end

PrefixHashScanLookupScratch() = PrefixHashScanLookupScratch(
    UInt64[], UInt64[], UInt64[], UInt64[], zeros(Int, 2048))

struct PrefixHashScanVerifiedHit
    guide_idx::Int
    pos::Int
    dist::Int
    is_antisense::Bool
    aln_guide::String
    aln_ref::String
end

struct PrefixHashScanChromResult{S}
    plus::Vector{PrefixHashScanVerifiedHit}
    minus::Vector{PrefixHashScanVerifiedHit}
    stats::S
end

struct PrefixHashScanChunkWork
    chrom_idx::Int
    core_first::Int
    core_last::Int
end

struct PrefixHashScanMyersProfile
    eq_by_iupac::NTuple{16, UInt64}
    length::UInt8
    final_bit::UInt64
end


include("prefix_hash_scan/query.jl")
include("prefix_hash_scan/verification.jl")
include("prefix_hash_scan/cas9.jl")
include("prefix_hash_scan/streaming.jl")

"""
    search_prefixHashScan(guides, genome_path, motif, output_file; kwargs...)

Indexless prototype that reuses prefixHashDB's symbolic edit-distance prefix
paths on the query side, scans motif-compatible genome candidates directly, and
exact-verifies candidate hits with `align`.
"""
function search_prefixHashScan(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif,
    output_file::String;
    distance::Int = 3,
    hash_len::Int = min(length_noPAM(motif) - distance, 16),
    early_stopping::Vector{Int} = fill(1_000_000, distance + 1),
    query_variant::Symbol = prefix_hash_scan_query_variant(),
    scan_backend::Symbol = :auto,
    bucket_bases::Int = 11,
    scan_threads::Int = Threads.nthreads(),
    stream_chunk_bases::Int = 2 * 1024 * 1024,
    prefilter_bits::Int = 26,
    query_build_backend::Symbol = :auto,
    lookup_variant::Symbol = :auto,
    verify_variant::Symbol = :auto,
    verbose::Bool = false,
    stats::Union{Nothing, PrefixHashScanStats} = nothing)

    if length(early_stopping) != (distance + 1)
        error("Specify one early stopping condition for a each distance, starting from distance 0.")
    end
    if any(length_noPAM(motif) .!= length.(guides))
        error("Guide queries are not of the correct length to use with this Motif: " * string(motif))
    end
    if hash_len < 1 || hash_len > 16
        error("hash_len must be in 1:16.")
    end
    scan_threads >= 1 || error("scan_threads must be positive.")
    stream_chunk_bases >= 64 ||
        error("stream_chunk_bases must be at least 64.")
    prefilter_bits in (0, 22, 24, 26) ||
        error("prefilter_bits must be 0, 22, 24, or 26.")
    query_build_backend in (:auto, :serial, :parallel) ||
        error("query_build_backend must be :auto, :serial, or :parallel.")
    lookup_variant in (:auto, :inline, :bucketed) ||
        error("lookup_variant must be :auto, :inline, or :bucketed.")
    verify_variant in (:auto, :align, :distance_first, :myers_raw) ||
        error("verify_variant must be :auto, :align, :distance_first, or :myers_raw.")
    if any(isambig.(guides))
        error("search_prefixHashScan does not support ambiguous query guides.")
    end

    if stats !== nothing
        reset!(stats)
    end

    metadata_start = prefix_hash_scan_timer(stats)
    dbi, reference_lengths = prefix_hash_scan_dbinfo(genome_path, motif)
    if stats !== nothing
        stats.metadata_ns += time_ns() - metadata_start
    end
    hash_type = smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
    resolved_query_variant =
        resolve_prefix_hash_scan_query_variant(query_variant, length(guides))
    scan_backend in (
        :auto, :legacy, :fused_dict, :fused_directory, :fused_fasta_simd,
        :streaming_fasta_simd, :streaming_fasta_simd_fused) ||
        error("scan_backend must be :auto, :legacy, :fused_dict, :fused_directory, :fused_fasta_simd, :streaming_fasta_simd, or :streaming_fasta_simd_fused.")
    geometry = CAS9_D3_PREFIX_SCAN_GEOMETRY
    supports_fused = distance == geometry.distance &&
        resolved_query_variant == :bitmask64 &&
        is_cas9_prefix_hash_candidate(dbi, hash_len)
    supports_fasta_simd = supports_fused && dbi.gi.is_fa &&
        hash_len == geometry.prefix_bases &&
        can_use_prefix_hash_scan_simd()
    resolved_scan_backend = if scan_backend == :auto
        supports_fasta_simd ? :streaming_fasta_simd :
            (supports_fused ? :fused_directory : :legacy)
    else
        scan_backend
    end
    if resolved_scan_backend != :legacy && !supports_fused
        error("Fused scan backends require Cas9 distance 3, hash_len <= 16, and at most 64 guides.")
    end
    if resolved_scan_backend in (
            :fused_fasta_simd, :streaming_fasta_simd,
            :streaming_fasta_simd_fused) &&
            !supports_fasta_simd
        resolved_scan_backend = :fused_directory
    end
    resolved_lookup_variant = if lookup_variant == :auto
        resolved_scan_backend == :streaming_fasta_simd &&
            prefilter_bits != 0 && bucket_bases == 11 &&
            stream_chunk_bases + prefix_scan_candidate_last_offset(
                CAS9_D3_PREFIX_SCAN_GEOMETRY) + distance <= typemax(UInt32) ?
            :bucketed : :inline
    else
        lookup_variant
    end
    if resolved_lookup_variant == :bucketed
        resolved_scan_backend == :streaming_fasta_simd ||
            error("lookup_variant=:bucketed requires the buffered streaming FASTA SIMD backend.")
        prefilter_bits != 0 ||
            error("lookup_variant=:bucketed requires a nonzero prefilter.")
        bucket_bases == 11 ||
            error("lookup_variant=:bucketed currently requires bucket_bases=11.")
        stream_chunk_bases + prefix_scan_candidate_last_offset(
            CAS9_D3_PREFIX_SCAN_GEOMETRY) + distance <= typemax(UInt32) ||
            error("lookup_variant=:bucketed requires chunks smaller than 4 GiB.")
    end
    if verify_variant == :myers_raw &&
            !(resolved_scan_backend in (
                :fused_fasta_simd, :streaming_fasta_simd,
                :streaming_fasta_simd_fused))
        error("verify_variant=:myers_raw requires a fused FASTA SIMD backend.")
    end
    if resolved_scan_backend in (
            :streaming_fasta_simd, :streaming_fasta_simd_fused) &&
            !(verify_variant in (:auto, :myers_raw))
        error("The streaming FASTA SIMD backend requires verify_variant=:auto or :myers_raw.")
    end
    if verbose
        worker_count = min(scan_threads, Threads.nthreads(), length(guides))
        resolved_query_build_backend = resolve_prefix_hash_scan_query_build_backend(
            query_build_backend, worker_count)
        scheduler = resolved_scan_backend in (
            :streaming_fasta_simd, :streaming_fasta_simd_fused) ?
            :chunk : :record
        @info(
            "prefixHashScan execution",
            scan_backend = resolved_scan_backend,
            lookup_variant = resolved_lookup_variant,
            query_build_backend = resolved_query_build_backend,
            scheduler = scheduler,
            scan_threads = scan_threads,
            chunk_bases = stream_chunk_bases,
        )
    end

    query_start = prefix_hash_scan_timer(stats)
    if resolved_scan_backend in (
            :fused_directory, :fused_fasta_simd, :streaming_fasta_simd,
            :streaming_fasta_simd_fused)
        query, guides_ = build_prefix_hash_scan_compact_query(
            guides,
            motif,
            distance,
            hash_len,
            stats;
            bucket_bases = bucket_bases,
            prefilter_bits = resolved_scan_backend in (
                :fused_fasta_simd, :streaming_fasta_simd,
                :streaming_fasta_simd_fused) ?
                prefilter_bits : 0,
            query_build_backend = query_build_backend,
            query_threads = scan_threads,
        )
    elseif resolved_query_variant == :bruteforce
        query = nothing
        guides_ = oriented_prefix_hash_scan_guides(guides, motif)
        if stats !== nothing
            stats.query_variant = :bruteforce
        end
    else
        query, guides_ = build_prefix_hash_scan_map(
            guides,
            motif,
            distance,
            hash_len,
            hash_type,
            stats;
            query_variant = resolved_query_variant,
        )
    end
    if stats !== nothing
        stats.scan_backend = resolved_scan_backend
        stats.query_build_ns += time_ns() - query_start
    end

    es_acc = zeros(Int, length(guides), length(early_stopping))
    is_es = falses(length(guides))
    seen = [Set{Tuple{String, Int, String, Int, String, String, String}}() for _ in guides]
    candidate_guides = Int[]
    use_bruteforce_query = resolved_query_variant == :bruteforce
    use_bitmask_query = query isa PrefixHashScanBitmaskQuery
    use_direct_cas9_hash = is_cas9_prefix_hash_candidate(dbi, hash_len)
    use_fused_scan = resolved_scan_backend in (
        :fused_dict, :fused_directory, :fused_fasta_simd,
        :streaming_fasta_simd, :streaming_fasta_simd_fused)
    use_myers_raw = verify_variant == :myers_raw ||
        (verify_variant == :auto && resolved_scan_backend in (
            :fused_fasta_simd, :streaming_fasta_simd,
            :streaming_fasta_simd_fused))
    use_distance_first = verify_variant == :distance_first ||
        (verify_variant == :auto && use_fused_scan && !use_myers_raw)
    myers_profiles = use_myers_raw ?
        build_prefix_hash_scan_myers_profiles(guides_) : nothing
    mkpath(dirname(output_file))
    open(output_file, "w") do out
        write(out, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")

        if resolved_scan_backend in (
                :streaming_fasta_simd, :streaming_fasta_simd_fused)
            scan_start = prefix_hash_scan_timer(stats)
            chunk_results, chrom_chunk_ranges = stream_prefix_hash_scan(
                dbi.gi.filepath,
                reference_lengths,
                query,
                dbi,
                guides_,
                myers_profiles,
                distance,
                stream_chunk_bases,
                scan_threads,
                Val(resolved_scan_backend == :streaming_fasta_simd_fused ?
                    :fused : (resolved_lookup_variant == :bucketed ?
                        :bucketed_reuse : :buffered_reuse)),
                stats,
            )
            for chrom_idx in eachindex(chrom_chunk_ranges)
                chrom_name = dbi.gi.chrom[chrom_idx]
                chunk_range = chrom_chunk_ranges[chrom_idx]
                if stats !== nothing
                    for chunk_idx in chunk_range
                        merge_prefix_hash_scan_worker_stats!(
                            stats, something(chunk_results[chunk_idx]).stats)
                    end
                end
                for strand in (:plus, :minus)
                    for chunk_idx in chunk_range
                        hits = getfield(something(chunk_results[chunk_idx]), strand)
                        for hit in hits
                            commit_prefix_hash_scan_verified!(
                                out,
                                hit,
                                guides[hit.guide_idx],
                                chrom_name,
                                early_stopping,
                                es_acc,
                                is_es,
                                seen,
                                stats,
                            )
                        end
                    end
                end
            end
            if stats !== nothing
                stats.scan_ns += time_ns() - scan_start
            end
            return nothing
        end

        ref = open(dbi.gi.filepath, "r")
        try
            reader = dbi.gi.is_fa ? FASTA.Reader(ref; copy = false) : TwoBit.Reader(ref)
            records = dbi.gi.is_fa ?
                PrefixHashScanFASTARecords(reader, dbi.gi.chrom, reference_lengths) :
                PrefixHashScanIndexedRecords(reader, dbi.gi.chrom)
            scan_start = time_ns()
            for (chrom_name, record, record_io_ns) in records
                if resolved_scan_backend == :fused_fasta_simd
                    raw_part = FASTX.seq_data_part(record, 1:FASTX.seqsize(record))
                    raw = @view record.data[raw_part]
                    if stats !== nothing
                        stats.record_io_ns += record_io_ns
                        stats.chrom_load_ns += record_io_ns
                    end
                    plus_hits, minus_hits = scan_cas9_prefix_hits_raw(
                        raw, dbi, query, stats; scan_threads = scan_threads)
                    for (is_antisense, hits) in ((false, plus_hits), (true, minus_hits))
                        if stats !== nothing
                            stats.prefix_hits += length(hits)
                            stats.guide_pairs +=
                                sum(hit -> count_ones(hit.mask), hits; init = 0)
                        end
                        for hit in hits
                            verify_prefix_hash_scan_bitmask_candidate!(
                                out,
                                raw,
                                hit.start:(hit.start + prefix_scan_candidate_last_offset(
                                    CAS9_D3_PREFIX_SCAN_GEOMETRY)),
                                dbi,
                                is_antisense,
                                hit.mask,
                                guides,
                                guides_,
                                chrom_name,
                                distance,
                                early_stopping,
                                es_acc,
                                is_es,
                                seen,
                                stats;
                                distance_first = use_distance_first,
                                myers_profiles = myers_profiles,
                            )
                        end
                    end
                    continue
                end

                convert_start = time_ns()
                chrom_seq = dbi.gi.is_fa ?
                    FASTA.sequence(LongDNA{4}, record) :
                    TwoBit.sequence(LongDNA{4}, record)
                sequence_convert_ns = time_ns() - convert_start
                if stats !== nothing
                    stats.record_io_ns += record_io_ns
                    stats.sequence_convert_ns += sequence_convert_ns
                    stats.chrom_load_ns += record_io_ns + sequence_convert_ns
                end

                if use_fused_scan
                    plus_hits, minus_hits = scan_cas9_prefix_hits(
                        chrom_seq, dbi, query, hash_len, stats; scan_threads = scan_threads)
                    for (is_antisense, hits) in ((false, plus_hits), (true, minus_hits))
                        if stats !== nothing
                            stats.prefix_hits += length(hits)
                            stats.guide_pairs += sum(hit -> count_ones(hit.mask), hits; init = 0)
                        end
                        for hit in hits
                            candidate_range = hit.start:(hit.start +
                                prefix_scan_candidate_last_offset(
                                    CAS9_D3_PREFIX_SCAN_GEOMETRY))
                            verify_prefix_hash_scan_bitmask_candidate!(
                                out,
                                chrom_seq,
                                candidate_range,
                                dbi,
                                is_antisense,
                                hit.mask,
                                guides,
                                guides_,
                                chrom_name,
                                distance,
                                early_stopping,
                                es_acc,
                                is_es,
                                seen,
                                stats;
                                distance_first = use_distance_first,
                            )
                        end
                    end
                    continue
                end

                for is_antisense in (false, true)
                    findguides_start = time_ns()
                    positions = findguides(dbi, chrom_seq, is_antisense)
                    if stats !== nothing
                        stats.findguides_ns += time_ns() - findguides_start
                    end
                    isempty(positions) && continue
                    strand = is_antisense ? "-" : "+"

                    for candidate_range in positions
                        if stats !== nothing
                            stats.motif_candidates += 1
                        end

                        candidate_mask = zero(UInt64)
                        if use_bruteforce_query
                            guide_count = count(!, is_es)
                            has_candidate_guides = guide_count != 0
                        else
                            hash_start = time_ns()
                            if use_direct_cas9_hash
                                hashes = candidate_prefix_hashes_direct_cas9(chrom_seq, candidate_range, is_antisense, hash_len, hash_type)
                            else
                                hashes = nothing
                            end
                            if hashes === nothing
                                if stats !== nothing
                                    stats.candidate_hash_ns += time_ns() - hash_start
                                end
                                prefix_start = time_ns()
                                prefix = normalized_candidate_prefix(chrom_seq, candidate_range, dbi, is_antisense, hash_len)
                                if stats !== nothing
                                    stats.candidate_prefix_ns += time_ns() - prefix_start
                                end
                                hash_start = time_ns()
                                hashes = candidate_prefix_hashes(prefix, hash_type, stats)
                            end
                            if stats !== nothing
                                stats.candidate_hash_ns += time_ns() - hash_start
                            end
                            lookup_start = time_ns()
                            if use_bitmask_query
                                candidate_mask = prefix_hash_scan_candidate_mask(query, hashes)
                                has_candidate_guides = candidate_mask != 0
                                guide_count = count_ones(candidate_mask)
                            else
                                has_candidate_guides = append_prefix_hash_scan_guides!(candidate_guides, query, hashes)
                                guide_count = length(candidate_guides)
                            end
                            if stats !== nothing
                                stats.query_lookup_ns += time_ns() - lookup_start
                            end
                        end
                        has_candidate_guides || continue

                        if stats !== nothing
                            stats.prefix_hits += 1
                            stats.guide_pairs += guide_count
                            if use_bruteforce_query
                                stats.bruteforce_guide_pairs += guide_count
                            end
                        end

                        materialize_start = time_ns()
                        ot, pos = materialize_normalized_candidate(chrom_seq, candidate_range, dbi, is_antisense)
                        if stats !== nothing
                            stats.candidate_materialize_ns += time_ns() - materialize_start
                        end
                        verify_start = time_ns()
                        if use_bitmask_query
                            mask = candidate_mask
                            while mask != 0
                                guide_idx = trailing_zeros(mask) + 1
                                mask &= mask - 1
                                is_es[guide_idx] && continue

                                if stats !== nothing
                                    stats.alignment_calls += 1
                                    stats.traceback_calls += 1
                                end
                                align_start = time_ns()
                                aln = align(guides_[guide_idx], ot, distance, iscompatible)
                                if stats !== nothing
                                    stats.align_ns += time_ns() - align_start
                                end
                                aln.dist > distance && continue

                                if motif.extends5
                                    aln_guide = reverse(aln.guide)
                                    aln_ref = reverse(aln.ref)
                                else
                                    aln_guide = aln.guide
                                    aln_ref = aln.ref
                                end

                                key = (
                                    string(guides[guide_idx]),
                                    aln.dist,
                                    chrom_name,
                                    Int(pos),
                                    strand,
                                    aln_guide,
                                    aln_ref,
                                )
                                key in seen[guide_idx] && continue
                                push!(seen[guide_idx], key)

                                emit_start = time_ns()
                                print(out, guides[guide_idx], ",", aln_guide, ",", aln_ref, ",",
                                    aln.dist, ",", chrom_name, ",", pos, ",", strand, "\n")
                                if stats !== nothing
                                    stats.emit_ns += time_ns() - emit_start
                                    stats.emitted_rows += 1
                                end
                                es_acc[guide_idx, aln.dist + 1] += 1
                                if es_acc[guide_idx, aln.dist + 1] >= early_stopping[aln.dist + 1]
                                    is_es[guide_idx] = true
                                end
                            end
                        elseif use_bruteforce_query
                            for guide_idx in eachindex(guides_)
                                is_es[guide_idx] && continue

                                if stats !== nothing
                                    stats.alignment_calls += 1
                                    stats.traceback_calls += 1
                                end
                                align_start = time_ns()
                                aln = align(guides_[guide_idx], ot, distance, iscompatible)
                                if stats !== nothing
                                    stats.align_ns += time_ns() - align_start
                                end
                                aln.dist > distance && continue

                                if motif.extends5
                                    aln_guide = reverse(aln.guide)
                                    aln_ref = reverse(aln.ref)
                                else
                                    aln_guide = aln.guide
                                    aln_ref = aln.ref
                                end

                                key = (
                                    string(guides[guide_idx]),
                                    aln.dist,
                                    chrom_name,
                                    Int(pos),
                                    strand,
                                    aln_guide,
                                    aln_ref,
                                )
                                key in seen[guide_idx] && continue
                                push!(seen[guide_idx], key)

                                emit_start = time_ns()
                                print(out, guides[guide_idx], ",", aln_guide, ",", aln_ref, ",",
                                    aln.dist, ",", chrom_name, ",", pos, ",", strand, "\n")
                                if stats !== nothing
                                    stats.emit_ns += time_ns() - emit_start
                                    stats.emitted_rows += 1
                                end
                                es_acc[guide_idx, aln.dist + 1] += 1
                                if es_acc[guide_idx, aln.dist + 1] >= early_stopping[aln.dist + 1]
                                    is_es[guide_idx] = true
                                end
                            end
                        else
                            for guide_idx in candidate_guides
                                is_es[guide_idx] && continue

                                if stats !== nothing
                                    stats.alignment_calls += 1
                                    stats.traceback_calls += 1
                                end
                                align_start = time_ns()
                                aln = align(guides_[guide_idx], ot, distance, iscompatible)
                                if stats !== nothing
                                    stats.align_ns += time_ns() - align_start
                                end
                                aln.dist > distance && continue

                                if motif.extends5
                                    aln_guide = reverse(aln.guide)
                                    aln_ref = reverse(aln.ref)
                                else
                                    aln_guide = aln.guide
                                    aln_ref = aln.ref
                                end

                                key = (
                                    string(guides[guide_idx]),
                                    aln.dist,
                                    chrom_name,
                                    Int(pos),
                                    strand,
                                    aln_guide,
                                    aln_ref,
                                )
                                key in seen[guide_idx] && continue
                                push!(seen[guide_idx], key)

                                emit_start = time_ns()
                                print(out, guides[guide_idx], ",", aln_guide, ",", aln_ref, ",",
                                    aln.dist, ",", chrom_name, ",", pos, ",", strand, "\n")
                                if stats !== nothing
                                    stats.emit_ns += time_ns() - emit_start
                                    stats.emitted_rows += 1
                                end
                                es_acc[guide_idx, aln.dist + 1] += 1
                                if es_acc[guide_idx, aln.dist + 1] >= early_stopping[aln.dist + 1]
                                    is_es[guide_idx] = true
                                end
                            end
                        end
                        if stats !== nothing
                            stats.verify_ns += time_ns() - verify_start
                        end
                    end
                end
            end
            if stats !== nothing
                stats.scan_ns += time_ns() - scan_start
            end
        finally
            close(ref)
        end
    end
    return
end

"""
    search_prefixHashScan(guides, genome_path, output_file;
        early_stopping=fill(1_000_000, 4), scan_threads=Threads.nthreads(),
        verbose=false)

Search a FASTA reference directly for Cas9 off-targets at edit distance 3. This
supported entrypoint accepts 1-64 unambiguous 20-base guides and requires the
standard FASTA `.fai` index. Candidate guide/PAM windows containing non-ACGT
bases are skipped.
"""
function search_prefixHashScan(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    output_file::String;
    early_stopping::Vector{Int} = fill(1_000_000, 4),
    scan_threads::Int = Threads.nthreads(),
    verbose::Bool = false)

    geometry = CAS9_D3_PREFIX_SCAN_GEOMETRY
    1 <= length(guides) <= 64 ||
        error("search_prefixHashScan supports 1-64 guides per search.")
    all(==(geometry.guide_bases), length.(guides)) ||
        error("search_prefixHashScan requires 20-base Cas9 guides.")
    any(isambig.(guides)) &&
        error("search_prefixHashScan does not support ambiguous query guides.")
    length(early_stopping) == geometry.distance + 1 ||
        error("Specify one early stopping condition for each distance from 0 to 3.")
    is_fasta(genome_path) ||
        error("search_prefixHashScan currently requires a FASTA reference.")

    motif = setdist(Motif("Cas9"), geometry.distance)
    return search_prefixHashScan(
        guides,
        genome_path,
        motif,
        output_file;
        distance = geometry.distance,
        hash_len = geometry.prefix_bases,
        early_stopping = early_stopping,
        query_variant = :bitmask64,
        scan_backend = :auto,
        bucket_bases = 11,
        scan_threads = scan_threads,
        stream_chunk_bases = 2 * 1024 * 1024,
        prefilter_bits = 26,
        query_build_backend = :auto,
        lookup_variant = :auto,
        verify_variant = :auto,
        verbose = verbose,
    )
end
