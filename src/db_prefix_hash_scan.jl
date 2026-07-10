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
end

PrefixHashScanStats() = PrefixHashScanStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, :none, :none)

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
    return stats
end

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

    load_start = time_ns()
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

struct PrefixHashScanHit
    start::Int
    mask::UInt64
end

struct PrefixHashScanFASTARecords{R, C}
    reader::R
    chrom::C
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
    FASTA.identifier(record) == chrom_name || error("FASTA record order differs from GenomeInfo.")
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
    FASTA.identifier(record) == chrom_name || error("FASTA record order differs from GenomeInfo.")
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
    hash_start = time_ns()
    for (guide_idx, guide) in enumerate(guides_)
        format_start = time_ns()
        guide_formatted = guide_to_template_format(guide; alphabet = ALPHABET_TWOBIT)
        if stats !== nothing
            stats.query_format_ns += time_ns() - format_start
        end

        if variant == :baseline
            fold_start = time_ns()
            hashes = prefix_hashes_baseline(paths, guide_formatted, hash_type)
            if stats !== nothing
                stats.query_fold_ns += time_ns() - fold_start
            end
        else
            hashes = Vector{hash_type}(undef, size(paths, 1))
            fold_start = time_ns()
            fill_prefix_hashes_columnwise!(hashes, paths, guide_formatted)
            if stats !== nothing
                stats.query_fold_ns += time_ns() - fold_start
            end

            dedup_start = time_ns()
            sort!(hashes)
            unique_sorted_prefix_hashes!(hashes)
            if stats !== nothing
                stats.query_dedup_ns += time_ns() - dedup_start
            end
        end

        if stats !== nothing
            stats.query_hashes += length(hashes)
        end

        insert_start = time_ns()
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

function materialize_normalized_candidate(
    chrom_seq::LongDNA{4},
    candidate_range::UnitRange{Int64},
    dbi::DBInfo,
    is_antisense::Bool)

    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    ot = removepam(chrom_seq[candidate_range], pam_loci)
    if dbi.motif.distance > 0
        ots = add_extension([ot], [candidate_range], dbi, chrom_seq, is_antisense)
        ot = ots[1]
    end
    ots, norm_pos = normalize_to_PAMseqEXT([ot], [candidate_range], dbi, is_antisense)
    return ots[1], norm_pos[1]
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
    if is_antisense
        start_pos = first(candidate_range) + 3
        @inbounds for pos in start_pos:(start_pos + hash_len - 1)
            code = dna_twobit_complement_base(chrom_seq[pos])
            code === nothing && return nothing
            h = (h << 2) | convert(hash_type, code)
        end
    else
        start_pos = first(candidate_range) + 19
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

@inline function prefix_hash_scan_twobit_nibble(nibble::UInt8)
    nibble == 0x01 && return UInt8(0)
    nibble == 0x02 && return UInt8(1)
    nibble == 0x04 && return UInt8(2)
    nibble == 0x08 && return UInt8(3)
    return UInt8(0xff)
end

function cas9_prefix_scan_bounds(chrom_seq::LongDNA{4}, dbi::DBInfo)
    n = length(chrom_seq)
    n < 23 && return nothing
    seq_start, seq_stop = locate_telomeres(chrom_seq)
    plus_first = max(seq_start - dbi.motif.distance, 1)
    plus_last = seq_stop - 22
    minus_first = seq_start
    minus_last = min(seq_stop + dbi.motif.distance, n) - 22
    firsts = Int[]
    lasts = Int[]
    if plus_first <= plus_last
        push!(firsts, plus_first)
        push!(lasts, plus_last)
    end
    if minus_first <= minus_last
        push!(firsts, minus_first)
        push!(lasts, minus_last)
    end
    isempty(firsts) && return nothing
    return minimum(firsts), maximum(lasts), plus_first, plus_last, minus_first, minus_last
end

function scan_cas9_prefix_hits_range(
    chrom_seq::LongDNA{4},
    query,
    hash_len::Int,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int)

    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    candidate_first > candidate_last && return plus_hits, minus_hits, 0

    hash_bits = 2 * hash_len
    hash_mask = (UInt64(1) << hash_bits) - UInt64(1)
    window_mask = (UInt64(1) << 46) - UInt64(1)
    hash_shift = 2 * (20 - hash_len)
    fwd_window = zero(UInt64)
    rev_window = zero(UInt64)
    valid_run = 0
    previous_code = UInt8(0xff)
    motif_candidates = 0

    @inbounds for pos in candidate_first:(candidate_first + 21)
        nibble = UInt8(BioSequences.extract_encoded_element(chrom_seq, pos))
        code = prefix_hash_scan_twobit_nibble(nibble)
        if code == 0xff
            valid_run = 0
            code = 0x00
        else
            valid_run += 1
        end
        fwd_window = ((fwd_window << 2) | UInt64(code)) & window_mask
        rev_window = (rev_window >> 2) | (UInt64(code) << 44)
        previous_code = code
    end

    @inbounds for pos in (candidate_first + 22):(candidate_last + 22)
        nibble = UInt8(BioSequences.extract_encoded_element(chrom_seq, pos))
        code = prefix_hash_scan_twobit_nibble(nibble)
        if code == 0xff
            valid_run = 0
            code = 0x00
        else
            valid_run += 1
        end

        fwd_window = ((fwd_window << 2) | UInt64(code)) & window_mask
        rev_window = (rev_window >> 2) | (UInt64(code) << 44)
        candidate_start = pos - 22
        if valid_run >= 23
            if plus_first <= candidate_start <= plus_last && code == 0x02 && previous_code == 0x02
                motif_candidates += 1
                hash = (rev_window >> hash_shift) & hash_mask
                mask = prefix_hash_scan_candidate_mask(query, hash)
                mask != 0 && push!(plus_hits, PrefixHashScanHit(candidate_start, mask))
            end
            if minus_first <= candidate_start <= minus_last &&
                    ((fwd_window >> 44) & UInt64(0x03)) == UInt64(0x01) &&
                    ((fwd_window >> 42) & UInt64(0x03)) == UInt64(0x01)
                motif_candidates += 1
                hash = xor((fwd_window >> hash_shift) & hash_mask, hash_mask)
                mask = prefix_hash_scan_candidate_mask(query, hash)
                mask != 0 && push!(minus_hits, PrefixHashScanHit(candidate_start, mask))
            end
        end
        previous_code = code
    end
    return plus_hits, minus_hits, motif_candidates
end

function scan_cas9_prefix_hits(
    chrom_seq::LongDNA{4},
    dbi::DBInfo,
    query,
    hash_len::Int,
    stats::Union{Nothing, PrefixHashScanStats} = nothing;
    scan_threads::Int = Threads.nthreads())

    bounds = cas9_prefix_scan_bounds(chrom_seq, dbi)
    bounds === nothing && return PrefixHashScanHit[], PrefixHashScanHit[]
    candidate_first, candidate_last, plus_first, plus_last, minus_first, minus_last = bounds
    candidate_count = candidate_last - candidate_first + 1
    thread_count = min(max(scan_threads, 1), candidate_count)

    if thread_count == 1
        plus_hits, minus_hits, motif_candidates = scan_cas9_prefix_hits_range(
            chrom_seq, query, hash_len, candidate_first, candidate_last,
            plus_first, plus_last, minus_first, minus_last)
    else
        chunk_size = cld(candidate_count, thread_count)
        ranges = [
            first:min(first + chunk_size - 1, candidate_last)
            for first in candidate_first:chunk_size:candidate_last
        ]
        tasks = map(ranges) do range
            Threads.@spawn scan_cas9_prefix_hits_range(
                chrom_seq, query, hash_len, first(range), last(range),
                plus_first, plus_last, minus_first, minus_last)
        end
        plus_hits = PrefixHashScanHit[]
        minus_hits = PrefixHashScanHit[]
        motif_candidates = 0
        for task in tasks
            local_plus, local_minus, local_candidates = fetch(task)
            append!(plus_hits, local_plus)
            append!(minus_hits, local_minus)
            motif_candidates += local_candidates
        end
    end

    if stats !== nothing
        stats.motif_candidates += motif_candidates
    end
    return plus_hits, minus_hits
end

function materialize_normalized_candidate_cas9(
    chrom_seq::LongDNA{4},
    candidate_start::Int,
    dbi::DBInfo,
    is_antisense::Bool)

    distance = dbi.motif.distance
    if is_antisense
        guide_start = candidate_start + 3
        extension_end = candidate_start + 22 + distance
        if extension_end <= length(chrom_seq)
            ot = complement(chrom_seq[guide_start:extension_end])
        else
            guide = chrom_seq[guide_start:(candidate_start + 22)]
            extension = getExt3(
                chrom_seq, length(chrom_seq), candidate_start + 23, distance)
            ot = complement(guide * extension)
        end
        return ot, candidate_start
    end

    extension_start = candidate_start - distance
    if extension_start >= 1
        ot = reverse(chrom_seq[extension_start:(candidate_start + 19)])
    else
        extension = getExt5(chrom_seq, candidate_start - 1, distance)
        guide = chrom_seq[candidate_start:(candidate_start + 19)]
        ot = reverse(extension * guide)
    end
    return ot, candidate_start + 22
end

function verify_prefix_hash_scan_bitmask_candidate!(
    out,
    chrom_seq::LongDNA{4},
    candidate_range::UnitRange{Int},
    dbi::DBInfo,
    is_antisense::Bool,
    candidate_mask::UInt64,
    guides::Vector{LongDNA{4}},
    guides_::Vector{LongDNA{4}},
    chrom_name::String,
    distance::Int,
    early_stopping::Vector{Int},
    es_acc,
    is_es,
    seen,
    stats::Union{Nothing, PrefixHashScanStats};
    distance_first::Bool = false)

    materialize_start = time_ns()
    ot, pos = materialize_normalized_candidate_cas9(
        chrom_seq, first(candidate_range), dbi, is_antisense)
    if stats !== nothing
        stats.candidate_materialize_ns += time_ns() - materialize_start
    end

    strand = is_antisense ? "-" : "+"
    verify_start = time_ns()
    mask = candidate_mask
    while mask != 0
        guide_idx = trailing_zeros(mask) + 1
        mask &= mask - 1
        is_es[guide_idx] && continue

        if stats !== nothing
            stats.alignment_calls += 1
        end
        align_start = time_ns()
        if distance_first
            if stats !== nothing
                stats.distance_calls += 1
            end
            dist = levenshtein(guides_[guide_idx], ot, distance, iscompatible)
            if dist > distance
                if stats !== nothing
                    stats.align_ns += time_ns() - align_start
                end
                continue
            end
        end
        if stats !== nothing
            stats.traceback_calls += 1
        end
        aln = align(guides_[guide_idx], ot, distance, iscompatible)
        if stats !== nothing
            stats.align_ns += time_ns() - align_start
        end
        aln.dist > distance && continue

        if dbi.motif.extends5
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
    if stats !== nothing
        stats.verify_ns += time_ns() - verify_start
    end
    return nothing
end

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
    verify_variant::Symbol = :auto,
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
    verify_variant in (:auto, :align, :distance_first) ||
        error("verify_variant must be :auto, :align, or :distance_first.")
    if any(isambig.(guides))
        error("search_prefixHashScan does not support ambiguous query guides.")
    end

    if stats !== nothing
        reset!(stats)
    end

    dbi = DBInfo(genome_path, "prefixHashScan", motif)
    hash_type = smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
    resolved_query_variant = resolve_prefix_hash_scan_query_variant(query_variant, length(guides))
    query = nothing
    query_start = time_ns()
    if resolved_query_variant == :bruteforce
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
    scan_backend in (:auto, :legacy, :fused_dict, :fused_directory) ||
        error("scan_backend must be :auto, :legacy, :fused_dict, or :fused_directory.")
    supports_fused = distance == 3 &&
        resolved_query_variant == :bitmask64 &&
        is_cas9_prefix_hash_candidate(dbi, hash_len)
    resolved_scan_backend = scan_backend == :auto ?
        (supports_fused ? :fused_directory : :legacy) : scan_backend
    if resolved_scan_backend != :legacy && !supports_fused
        error("Fused scan backends require Cas9 distance 3, hash_len <= 16, and at most 64 guides.")
    end
    if resolved_scan_backend == :fused_directory
        query = build_prefix_hash_scan_directory(query, hash_len, bucket_bases)
    end
    if stats !== nothing
        stats.query_build_ns += time_ns() - query_start
    end

    es_acc = zeros(Int, length(guides), length(early_stopping))
    is_es = falses(length(guides))
    seen = [Set{Tuple{String, Int, String, Int, String, String, String}}() for _ in guides]
    candidate_guides = Int[]
    use_bruteforce_query = resolved_query_variant == :bruteforce
    use_bitmask_query = query isa PrefixHashScanBitmaskQuery
    use_direct_cas9_hash = is_cas9_prefix_hash_candidate(dbi, hash_len)
    use_fused_scan = resolved_scan_backend in (:fused_dict, :fused_directory)
    use_distance_first = verify_variant == :distance_first ||
        (verify_variant == :auto && use_fused_scan)
    mkpath(dirname(output_file))
    open(output_file, "w") do out
        write(out, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")

        ref = open(dbi.gi.filepath, "r")
        try
            reader = dbi.gi.is_fa ? FASTA.Reader(ref) : TwoBit.Reader(ref)
            records = dbi.gi.is_fa ?
                PrefixHashScanFASTARecords(reader, dbi.gi.chrom) :
                PrefixHashScanIndexedRecords(reader, dbi.gi.chrom)
            scan_start = time_ns()
            for (chrom_name, record, record_io_ns) in records
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
                            candidate_range = hit.start:(hit.start + 22)
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
