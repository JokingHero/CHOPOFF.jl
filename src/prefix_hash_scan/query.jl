# Query construction and compact lookup directory.

function is_precomputed_cas9_prefix_paths(motif::Motif, hash_len::Int)
    cas9 = Motif("Cas9")
    return hash_len <= 16 &&
        motif.distance <= 4 &&
        length_noPAM(motif) == length_noPAM(cas9) &&
        motif.extends5 == cas9.extends5 &&
        motif.pam_loci_fwd == cas9.pam_loci_fwd &&
        motif.pam_loci_rve == cas9.pam_loci_rve
end

function dedup_prefix_paths(paths, distances, hash_len::Int, distance::Int)
    paths = paths[:, 1:hash_len]
    not_dups = map(!, BitVector(nonunique(DataFrame(paths, :auto))))
    not_over_dist = BitVector(distances .<= distance)
    keep = not_dups .& not_over_dist
    paths = paths[keep, :]
    return convert.(smallestutype(maximum(paths)), paths)
end

function load_prefix_hash_scan_paths(
    motif::Motif,
    distance::Int,
    hash_len::Int,
    stats::Union{Nothing, PrefixHashScanStats} = nothing)

    load_start = prefix_hash_scan_timer(stats)
    source = :generated
    paths = nothing

    precomputed = load_precomputed_prefix_paths(motif, distance, hash_len; need_distances = false)
    if precomputed !== nothing
        paths, _, _ = precomputed
        source = :precomputed
    end

    if paths === nothing
        mpt = build_PathTemplates(motif; restrict_to_len = hash_len, withPAM = false)
        paths = mpt.paths[mpt.distances .<= distance, 1:hash_len]
        not_dups = map(!, BitVector(nonunique(DataFrame(paths, :auto))))
        paths = paths[not_dups, :]
        paths = convert.(smallestutype(maximum(paths)), paths)
    end

    if stats !== nothing
        stats.path_rows = size(paths, 1)
        stats.path_source = source
        stats.path_load_ns += time_ns() - load_start
    end
    return paths, source
end

function fold_prefix_hash(path_row, guide_formatted, hash_type::Type{<:Unsigned})
    h = zero(hash_type)
    @inbounds for path_idx in path_row
        h = (h << 2) | convert(hash_type, guide_formatted[Int(path_idx)])
    end
    return h
end

function resolve_prefix_hash_scan_query_variant(query_variant::Symbol, nguides::Int = 0)
    query_variant == :auto && return nguides <= 64 ? :bitmask64 : :columnwise
    query_variant == :bruteforce && return :bruteforce
    if query_variant == :bitmask64 && nguides > 64
        error("query_variant=:bitmask64 supports at most 64 guides.")
    end
    query_variant in (:baseline, :columnwise, :bitmask64) && return query_variant
    error("query_variant must be :auto, :baseline, :columnwise, :bitmask64, or :bruteforce.")
end

function prefix_hash_scan_query_variant()
    raw = lowercase(strip(get(ENV, "CHOPOFF_PREFIX_HASH_SCAN_QUERY", "auto")))
    raw == "auto" && return :auto
    raw == "baseline" && return :baseline
    raw == "columnwise" && return :columnwise
    raw == "bitmask64" && return :bitmask64
    raw == "bruteforce" && return :bruteforce
    error("Invalid CHOPOFF_PREFIX_HASH_SCAN_QUERY='$raw'. Allowed values: auto, baseline, columnwise, bitmask64, bruteforce.")
end

function oriented_prefix_hash_scan_guides(guides::Vector{LongDNA{4}}, motif::Motif)
    guides_ = copy(guides)
    if motif.extends5
        guides_ = reverse.(guides_)
    end
    return guides_
end

function unique_sorted_prefix_hashes!(hashes::Vector{T}) where {T <: Unsigned}
    isempty(hashes) && return hashes
    write_idx = 1
    @inbounds for read_idx in 2:length(hashes)
        if hashes[read_idx] != hashes[write_idx]
            write_idx += 1
            hashes[write_idx] = hashes[read_idx]
        end
    end
    resize!(hashes, write_idx)
    return hashes
end

function prefix_hashes_baseline(
    paths,
    guide_formatted,
    hash_type::Type{<:Unsigned})

    hashes = Set{hash_type}()
    for path_row in eachrow(paths)
        push!(hashes, fold_prefix_hash(path_row, guide_formatted, hash_type))
    end
    hashes = collect(hashes)
    sort!(hashes)
    return hashes
end

function fill_prefix_hashes_columnwise!(
    hashes::Vector{T},
    paths,
    guide_formatted) where {T <: Unsigned}

    fill!(hashes, zero(T))
    npaths, hash_len = size(paths)
    @inbounds for col in 1:hash_len
        for row in 1:npaths
            hashes[row] = (hashes[row] << 2) | convert(T, guide_formatted[Int(paths[row, col])])
        end
    end
    return hashes
end

function prefix_hashes_columnwise(
    paths,
    guide_formatted,
    hash_type::Type{<:Unsigned})

    hashes = Vector{hash_type}(undef, size(paths, 1))
    fill_prefix_hashes_columnwise!(hashes, paths, guide_formatted)
    sort!(hashes)
    unique_sorted_prefix_hashes!(hashes)
    return hashes
end

function prefix_hash_scan_guide_hashes(
    paths,
    guide::LongDNA{4},
    hash_type::Type{<:Unsigned};
    query_variant::Symbol = prefix_hash_scan_query_variant())

    variant = resolve_prefix_hash_scan_query_variant(query_variant)
    guide_formatted = guide_to_template_format(guide; alphabet = ALPHABET_TWOBIT)
    if variant == :baseline
        return prefix_hashes_baseline(paths, guide_formatted, hash_type)
    end
    return prefix_hashes_columnwise(paths, guide_formatted, hash_type)
end

function build_prefix_hash_scan_map_from_paths(
    paths,
    guides_::Vector{LongDNA{4}},
    hash_type::Type{<:Unsigned},
    stats::Union{Nothing, PrefixHashScanStats} = nothing;
    query_variant::Symbol = prefix_hash_scan_query_variant())

    variant = resolve_prefix_hash_scan_query_variant(query_variant, length(guides_))
    if stats !== nothing
        stats.query_variant = variant
    end

    if variant == :bitmask64
        query = PrefixHashScanBitmaskQuery(Dict{hash_type, UInt64}(), length(guides_))
    else
        query = Dict{hash_type, Vector{Int}}()
    end
    hash_start = prefix_hash_scan_timer(stats)
    for (guide_idx, guide) in enumerate(guides_)
        format_start = prefix_hash_scan_timer(stats)
        guide_formatted = guide_to_template_format(guide; alphabet = ALPHABET_TWOBIT)
        if stats !== nothing
            stats.query_format_ns += time_ns() - format_start
        end

        if variant == :baseline
            fold_start = prefix_hash_scan_timer(stats)
            hashes = prefix_hashes_baseline(paths, guide_formatted, hash_type)
            if stats !== nothing
                stats.query_fold_ns += time_ns() - fold_start
            end
        else
            hashes = Vector{hash_type}(undef, size(paths, 1))
            fold_start = prefix_hash_scan_timer(stats)
            fill_prefix_hashes_columnwise!(hashes, paths, guide_formatted)
            if stats !== nothing
                stats.query_fold_ns += time_ns() - fold_start
            end

            dedup_start = prefix_hash_scan_timer(stats)
            sort!(hashes)
            unique_sorted_prefix_hashes!(hashes)
            if stats !== nothing
                stats.query_dedup_ns += time_ns() - dedup_start
            end
        end

        if stats !== nothing
            stats.query_hashes += length(hashes)
        end

        insert_start = prefix_hash_scan_timer(stats)
        if variant == :bitmask64
            bit = UInt64(1) << (guide_idx - 1)
            for h in hashes
                query.masks[h] = get(query.masks, h, zero(UInt64)) | bit
            end
        else
            for h in hashes
                push!(get!(query, h, Int[]), guide_idx)
            end
        end
        if stats !== nothing
            stats.query_insert_ns += time_ns() - insert_start
        end
    end
    if stats !== nothing
        stats.query_hash_ns += time_ns() - hash_start
    end
    return query
end

function build_prefix_hash_scan_map(
    guides::Vector{LongDNA{4}},
    motif::Motif,
    distance::Int,
    hash_len::Int,
    hash_type::Type{<:Unsigned},
    stats::Union{Nothing, PrefixHashScanStats} = nothing;
    query_variant::Symbol = prefix_hash_scan_query_variant())

    paths, _ = load_prefix_hash_scan_paths(motif, distance, hash_len, stats)
    guides_ = oriented_prefix_hash_scan_guides(guides, motif)
    query = build_prefix_hash_scan_map_from_paths(paths, guides_, hash_type, stats; query_variant = query_variant)
    return query, guides_
end

function normalized_candidate_prefix(
    chrom_seq::LongDNA{4},
    candidate_range::UnitRange{Int64},
    dbi::DBInfo,
    is_antisense::Bool,
    hash_len::Int)

    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    guide = removepam(chrom_seq[candidate_range], pam_loci)
    if dbi.motif.extends5 && is_antisense
        prefix_source = complement(guide)
    elseif dbi.motif.extends5 && !is_antisense
        prefix_source = reverse(guide)
    elseif !dbi.motif.extends5 && is_antisense
        prefix_source = reverse_complement(guide)
    else
        prefix_source = guide
    end
    return prefix_source[1:hash_len]
end

function candidate_prefix_hashes(
    prefix::LongDNA{4},
    hash_type::Type{<:Unsigned},
    stats::Union{Nothing, PrefixHashScanStats})

    if isambig(prefix)
        if stats !== nothing
            stats.ambiguous_prefixes += 1
        end
        expanded, _ = expand_ambiguous(prefix)
        return unique(convert.(hash_type, expanded))
    end
    return hash_type[convert(hash_type, prefix)]
end

function dna_twobit_base(base)
    base == DNA_A && return UInt8(0)
    base == DNA_C && return UInt8(1)
    base == DNA_G && return UInt8(2)
    base == DNA_T && return UInt8(3)
    return nothing
end

function dna_twobit_complement_base(base)
    base == DNA_A && return UInt8(3)
    base == DNA_C && return UInt8(2)
    base == DNA_G && return UInt8(1)
    base == DNA_T && return UInt8(0)
    return nothing
end

function is_cas9_prefix_hash_candidate(dbi::DBInfo, hash_len::Int)
    cas9 = Motif("Cas9")
    return hash_len <= 16 &&
        length_noPAM(dbi.motif) == length_noPAM(cas9) &&
        dbi.motif.extends5 == cas9.extends5 &&
        dbi.motif.pam_loci_fwd == cas9.pam_loci_fwd &&
        dbi.motif.pam_loci_rve == cas9.pam_loci_rve
end

function candidate_prefix_hashes_direct_cas9(
    chrom_seq::LongDNA{4},
    candidate_range::UnitRange{Int64},
    is_antisense::Bool,
    hash_len::Int,
    hash_type::Type{<:Unsigned})

    h = zero(hash_type)
    geometry = CAS9_D3_PREFIX_SCAN_GEOMETRY
    if is_antisense
        start_pos = first(candidate_range) + geometry.pam_bases
        @inbounds for pos in start_pos:(start_pos + hash_len - 1)
            code = dna_twobit_complement_base(chrom_seq[pos])
            code === nothing && return nothing
            h = (h << 2) | convert(hash_type, code)
        end
    else
        start_pos = first(candidate_range) + geometry.guide_bases - 1
        @inbounds for pos in start_pos:-1:(start_pos - hash_len + 1)
            code = dna_twobit_base(chrom_seq[pos])
            code === nothing && return nothing
            h = (h << 2) | convert(hash_type, code)
        end
    end
    return hash_type[h]
end

function append_prefix_hash_scan_guides!(candidate_guides::Vector{Int}, query::Dict, hashes)
    empty!(candidate_guides)
    for h in hashes
        append!(candidate_guides, get(query, h, Int[]))
    end
    isempty(candidate_guides) && return false
    sort!(candidate_guides)
    unique!(candidate_guides)
    return true
end

function prefix_hash_scan_candidate_mask(query::PrefixHashScanBitmaskQuery, hashes)
    mask = zero(UInt64)
    for h in hashes
        mask |= get(query.masks, h, zero(UInt64))
    end
    return mask
end

@inline prefix_hash_scan_candidate_mask(query::PrefixHashScanBitmaskQuery, hash::Unsigned) =
    get(query.masks, hash, zero(UInt64))

function build_prefix_hash_scan_directory(
    query::PrefixHashScanBitmaskQuery,
    hash_len::Int,
    bucket_bases::Int)

    1 <= bucket_bases < hash_len || error("bucket_bases must be in 1:(hash_len - 1).")
    suffix_bits = 2 * (hash_len - bucket_bases)
    suffix_bits <= 16 || error("Directory suffix must fit in UInt16.")

    keys_ = sort!(UInt32.(collect(keys(query.masks))))
    nbuckets = 1 << (2 * bucket_bases)
    offsets = Vector{UInt32}(undef, nbuckets + 1)
    suffixes = Vector{UInt16}(undef, length(keys_))
    masks = Vector{UInt64}(undef, length(keys_))
    suffix_mask = (UInt32(1) << suffix_bits) - UInt32(1)

    key_idx = 1
    @inbounds for bucket in 0:(nbuckets - 1)
        offsets[bucket + 1] = UInt32(key_idx - 1)
        while key_idx <= length(keys_) && Int(keys_[key_idx] >> suffix_bits) == bucket
            key = keys_[key_idx]
            suffixes[key_idx] = UInt16(key & suffix_mask)
            masks[key_idx] = query.masks[key]
            key_idx += 1
        end
    end
    offsets[end] = UInt32(length(keys_))
    return PrefixHashScanDirectory(offsets, suffixes, masks, UInt8(hash_len), UInt8(bucket_bases))
end


function merge_prefix_hash_scan_hash_lists(
    lists::Vector{Vector{T}}) where T <: Unsigned

    positions = ones(Int, length(lists))
    heap = [idx for idx in eachindex(lists) if !isempty(lists[idx])]
    @inline less(left, right) =
        @inbounds lists[left][positions[left]] < lists[right][positions[right]]

    function sift_down!(idx)
        while true
            left = idx << 1
            left > length(heap) && return
            right = left + 1
            child = right <= length(heap) && less(heap[right], heap[left]) ?
                right : left
            less(heap[child], heap[idx]) || return
            heap[idx], heap[child] = heap[child], heap[idx]
            idx = child
        end
    end

    for idx in (length(heap) >> 1):-1:1
        sift_down!(idx)
    end
    capacity = sum(length, lists)
    keys_ = T[]
    masks = UInt64[]
    sizehint!(keys_, capacity)
    sizehint!(masks, capacity)

    while !isempty(heap)
        hash = @inbounds lists[heap[1]][positions[heap[1]]]
        mask = UInt64(0)
        while !isempty(heap) &&
                (@inbounds lists[heap[1]][positions[heap[1]]]) == hash
            guide_idx = heap[1]
            mask |= UInt64(1) << (guide_idx - 1)
            positions[guide_idx] += 1
            if positions[guide_idx] > length(lists[guide_idx])
                heap[1] = heap[end]
                pop!(heap)
            else
                heap[1] = guide_idx
            end
            isempty(heap) || sift_down!(1)
        end
        push!(keys_, hash)
        push!(masks, mask)
    end
    return keys_, masks
end

function build_prefix_hash_scan_directory(
    keys_::Vector{UInt32},
    masks::Vector{UInt64},
    hash_len::Int,
    bucket_bases::Int)

    length(keys_) == length(masks) || error("Hash and mask counts differ.")
    1 <= bucket_bases < hash_len ||
        error("bucket_bases must be in 1:(hash_len - 1).")
    suffix_bits = 2 * (hash_len - bucket_bases)
    suffix_bits <= 16 || error("Directory suffix must fit in UInt16.")

    nbuckets = 1 << (2 * bucket_bases)
    offsets = Vector{UInt32}(undef, nbuckets + 1)
    suffixes = Vector{UInt16}(undef, length(keys_))
    suffix_mask = (UInt32(1) << suffix_bits) - UInt32(1)
    key_idx = 1
    @inbounds for bucket in 0:(nbuckets - 1)
        offsets[bucket + 1] = UInt32(key_idx - 1)
        while key_idx <= length(keys_) &&
                Int(keys_[key_idx] >> suffix_bits) == bucket
            suffixes[key_idx] = UInt16(keys_[key_idx] & suffix_mask)
            key_idx += 1
        end
    end
    offsets[end] = UInt32(length(keys_))
    return PrefixHashScanDirectory(
        offsets, suffixes, masks, UInt8(hash_len), UInt8(bucket_bases))
end

function build_prefix_hash_scan_compact_guide_hashes(
    paths, guide::LongDNA{4}, ::Val{Timed}) where Timed

    format_start = Timed ? time_ns() : UInt64(0)
    formatted = guide_to_template_format(
        guide; alphabet = ALPHABET_TWOBIT)
    format_ns = Timed ? time_ns() - format_start : UInt64(0)

    hashes = Vector{UInt32}(undef, size(paths, 1))
    fold_start = Timed ? time_ns() : UInt64(0)
    fill_prefix_hashes_columnwise!(hashes, paths, formatted)
    fold_ns = Timed ? time_ns() - fold_start : UInt64(0)

    dedup_start = Timed ? time_ns() : UInt64(0)
    sort!(hashes)
    unique_sorted_prefix_hashes!(hashes)
    dedup_ns = Timed ? time_ns() - dedup_start : UInt64(0)
    return hashes, (format_ns, fold_ns, dedup_ns)
end

function build_prefix_hash_scan_compact_lists!(
    lists, timings, paths, guides_, worker_count::Int, ::Val{Parallel},
    timed::Val) where Parallel

    function build_guide!(guide_idx)
        hashes, guide_timings = build_prefix_hash_scan_compact_guide_hashes(
            paths, guides_[guide_idx], timed)
        lists[guide_idx] = hashes
        timings === nothing || (timings[guide_idx] = guide_timings)
        return nothing
    end

    if Parallel && worker_count > 1
        workers = map(1:worker_count) do worker_idx
            Threads.@spawn for guide_idx in worker_idx:worker_count:length(guides_)
                build_guide!(guide_idx)
            end
        end
        fetch.(workers)
    else
        for guide_idx in eachindex(guides_)
            build_guide!(guide_idx)
        end
    end
    return nothing
end

@inline function resolve_prefix_hash_scan_query_build_backend(
    query_build_backend::Symbol,
    worker_count::Int)

    return query_build_backend == :auto && worker_count > 1 ?
        :parallel : (query_build_backend == :auto ? :serial : query_build_backend)
end

function build_prefix_hash_scan_compact_query(
    guides::Vector{LongDNA{4}},
    motif::Motif,
    distance::Int,
    hash_len::Int,
    stats::Union{Nothing, PrefixHashScanStats};
    bucket_bases::Int,
    prefilter_bits::Int,
    query_build_backend::Symbol = :serial,
    query_threads::Int = 1)

    query_build_backend in (:auto, :serial, :parallel) ||
        error("query_build_backend must be :auto, :serial, or :parallel.")
    query_threads >= 1 || error("query_threads must be positive.")

    paths, _ = load_prefix_hash_scan_paths(
        motif, distance, hash_len, stats)
    guides_ = oriented_prefix_hash_scan_guides(guides, motif)
    lists = Vector{Vector{UInt32}}(undef, length(guides_))
    timings = stats === nothing ? nothing :
        Vector{NTuple{3, UInt64}}(undef, length(guides_))
    worker_count = min(query_threads, Threads.nthreads(), length(guides_))
    resolved_backend = resolve_prefix_hash_scan_query_build_backend(
        query_build_backend, worker_count)
    build_prefix_hash_scan_compact_lists!(
        lists, timings, paths, guides_, worker_count,
        Val(resolved_backend == :parallel), Val(stats !== nothing))
    if stats !== nothing
        stats.query_format_ns += sum(first, timings)
        stats.query_fold_ns += sum(x -> x[2], timings)
        stats.query_dedup_ns += sum(last, timings)
        stats.query_hashes += sum(length, lists)
    end

    merge_start = prefix_hash_scan_timer(stats)
    keys_, masks = merge_prefix_hash_scan_hash_lists(lists)
    directory = build_prefix_hash_scan_directory(
        keys_, masks, hash_len, bucket_bases)
    if prefilter_bits != 0
        directory = build_prefix_hash_scan_prefilter(
            directory, keys_, prefilter_bits)
    end
    if stats !== nothing
        stats.query_insert_ns += time_ns() - merge_start
        stats.query_variant = :bitmask64
    end
    return directory, guides_
end

@inline function prefix_hash_scan_candidate_mask(query::PrefixHashScanDirectory, hash::Unsigned)
    suffix_bits = 2 * (Int(query.hash_len) - Int(query.bucket_bases))
    key = UInt32(hash)
    bucket = Int(key >> suffix_bits)
    suffix_mask = (UInt32(1) << suffix_bits) - UInt32(1)
    suffix = UInt16(key & suffix_mask)
    first_idx = Int(@inbounds query.offsets[bucket + 1]) + 1
    last_idx = Int(@inbounds query.offsets[bucket + 2])
    mask = zero(UInt64)
    @inbounds for idx in first_idx:last_idx
        query.suffixes[idx] == suffix && (mask |= query.masks[idx])
    end
    return mask
end


function build_prefix_hash_scan_prefilter(
    directory::PrefixHashScanDirectory,
    hashes,
    prefix_bits::Int)

    prefix_bits in (22, 24, 26) ||
        error("prefilter_bits must be 22, 24, or 26.")
    presence = zeros(UInt64, 1 << (prefix_bits - 6))
    shift = 32 - prefix_bits
    for hash in hashes
        prefix = Int(UInt32(hash) >> shift)
        @inbounds presence[(prefix >> 6) + 1] |= UInt64(1) << (prefix & 63)
    end
    return PrefixHashScanPrefilteredDirectory(
        directory, presence, UInt8(prefix_bits))
end

@inline function prefix_hash_scan_candidate_mask(
    query::PrefixHashScanPrefilteredDirectory,
    hash::Unsigned)

    prefix_hash_scan_prefilter_contains(query, hash) || return UInt64(0)
    return prefix_hash_scan_candidate_mask(query.directory, hash)
end

@inline function prefix_hash_scan_prefilter_contains(
    query::PrefixHashScanPrefilteredDirectory,
    hash::Unsigned)

    shift = 32 - Int(query.prefix_bits)
    prefix = Int(UInt32(hash) >> shift)
    word = @inbounds query.presence[(prefix >> 6) + 1]
    return ((word >> (prefix & 63)) & UInt64(1)) != 0
end
