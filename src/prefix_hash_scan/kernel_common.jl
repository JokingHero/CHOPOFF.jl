# Geometry-neutral scalar and compact-lookup helpers. ISA primitives live in isa.jl.

@inline function prefix_hash_scan_twobit_nibble(nibble::UInt8)
    nibble == 0x01 && return UInt8(0)
    nibble == 0x02 && return UInt8(1)
    nibble == 0x04 && return UInt8(2)
    nibble == 0x08 && return UInt8(3)
    return UInt8(0xff)
end

@inline function prefix_hash_scan_reverse_codes(hash::UInt32)
    reversed = bitreverse(hash)
    return ((reversed & UInt32(0xaaaaaaaa)) >> 1) |
        ((reversed & UInt32(0x55555555)) << 1)
end

@inline function prefix_hash_scan_raw_code(base::UInt8)
    (base == UInt8('A') || base == UInt8('a')) && return UInt8(0)
    (base == UInt8('C') || base == UInt8('c')) && return UInt8(1)
    (base == UInt8('G') || base == UInt8('g')) && return UInt8(2)
    (base == UInt8('T') || base == UInt8('t')) && return UInt8(3)
    return UInt8(0xff)
end

@inline prefix_hash_scan_hash_type(
    ::PrefixHashScanBitmaskQuery{H}) where {H <: Unsigned} = H
@inline prefix_hash_scan_hash_type(::PrefixHashScanDirectory) = UInt32
@inline prefix_hash_scan_hash_type(::PrefixHashScanPrefilteredDirectory) = UInt32

@inline function prefix_hash_scan_ambiguous_window_mask(
    exact::UInt128, window_bases::Int, ::Val{K}) where K

    one_or_more = UInt64(0)
    two_or_more = UInt64(0)
    three_or_more = UInt64(0)
    four_or_more = UInt64(0)
    ambiguous = ~exact
    @inbounds for offset in 0:(window_bases - 1)
        present = UInt64(
            (ambiguous >> offset) & UInt128(typemax(UInt64)))
        four_or_more |= three_or_more & present
        three_or_more |= two_or_more & present
        two_or_more |= one_or_more & present
        one_or_more |= present
    end
    too_many = K == 1 ? two_or_more : K == 2 ? three_or_more : four_or_more
    return one_or_more & ~too_many
end

@inline function prefix_hash_scan_exact_block(
    raw::AbstractVector{UInt8}, start_pos::Int, required_bases::Int,
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    if start_pos + 127 <= length(raw)
        a0, c0, g0, t0 = prefix_hash_scan_raw_profile64(
            raw, start_pos, simd_backend)
        a1, c1, g1, t1 = prefix_hash_scan_raw_profile64(
            raw, start_pos + 64, simd_backend)
        return UInt128(a0 | c0 | g0 | t0) |
            (UInt128(a1 | c1 | g1 | t1) << 64)
    end
    exact = UInt128(0)
    @inbounds for offset in 0:(required_bases - 1)
        prefix_hash_scan_raw_code(raw[start_pos + offset]) != 0xff &&
            (exact |= UInt128(1) << offset)
    end
    return exact
end

@inline function prefix_hash_scan_exact_block(
    chrom_seq::LongDNA{4}, start_pos::Int, required_bases::Int,
    ::Val = Val(:avx2))

    exact = UInt128(0)
    @inbounds for offset in 0:(required_bases - 1)
        nibble = UInt8(BioSequences.extract_encoded_element(
            chrom_seq, start_pos + offset))
        prefix_hash_scan_twobit_nibble(nibble) != 0xff &&
            (exact |= UInt128(1) << offset)
    end
    return exact
end

@inline function prefix_hash_scan_candidate_sequence(
    raw::AbstractVector{UInt8}, candidate_start::Int, candidate_bases::Int)

    @inbounds for pos in candidate_start:(candidate_start + candidate_bases - 1)
        prefix_hash_scan_iupac_mask(raw[pos]) == 0 && return nothing
    end
    return LongDNA{4}(
        @view raw[candidate_start:(candidate_start + candidate_bases - 1)])
end
@inline prefix_hash_scan_candidate_sequence(
    chrom_seq::LongDNA{4}, candidate_start::Int, candidate_bases::Int) =
    chrom_seq[candidate_start:(candidate_start + candidate_bases - 1)]

function prefix_hash_scan_window_prefix(
    candidate::LongDNA{4}, motif::Motif, is_antisense::Bool, hash_len::Int)

    pam_loci = is_antisense ? motif.pam_loci_rve : motif.pam_loci_fwd
    guide = removepam(candidate, pam_loci)
    if motif.extends5 && is_antisense
        guide = complement(guide)
    elseif motif.extends5 && !is_antisense
        guide = reverse(guide)
    elseif !motif.extends5 && is_antisense
        guide = reverse_complement(guide)
    end
    return guide[1:hash_len]
end

function prefix_hash_scan_ambiguous_candidate_mask(
    candidate::LongDNA{4}, motif::Motif, is_antisense::Bool,
    hash_len::Int, query)

    pattern = is_antisense ? motif.rve : motif.fwd
    @inbounds for idx in eachindex(pattern)
        iscompatible(pattern[idx], candidate[idx]) || return false, UInt64(0), false
    end
    prefix = prefix_hash_scan_window_prefix(
        candidate, motif, is_antisense, hash_len)
    hashes = candidate_prefix_hashes(
        prefix, prefix_hash_scan_hash_type(query), nothing)
    mask = UInt64(0)
    @inbounds for hash in hashes
        mask |= prefix_hash_scan_candidate_mask(query, hash)
    end
    return true, mask, isambig(prefix)
end

function scan_ambiguous_prefix_hits_range!(
    plus_hits::Vector{PrefixHashScanHit},
    minus_hits::Vector{PrefixHashScanHit},
    source,
    geometry::PrefixScanGeometry,
    dbi::DBInfo,
    query,
    hash_len::Int,
    bounds::PrefixScanBounds,
    ::Val{K},
    stats::Union{Nothing, PrefixHashScanStats} = nothing,
    simd_backend::Val = default_prefix_hash_scan_simd_backend()) where K

    candidate_first, candidate_last = first(bounds.all), last(bounds.all)
    candidate_first > candidate_last && return nothing
    candidate_bases = prefix_scan_candidate_bases(geometry)
    for block_start in candidate_first:64:candidate_last
        count = min(64, candidate_last - block_start + 1)
        required_bases = count + candidate_bases - 1
        exact = prefix_hash_scan_exact_block(
            source, block_start, required_bases, simd_backend)
        required_mask = (UInt128(1) << required_bases) - UInt128(1)
        exact & required_mask == required_mask && continue
        starts = prefix_hash_scan_ambiguous_window_mask(
            exact, candidate_bases, Val(K))
        count < 64 && (starts &= (UInt64(1) << count) - UInt64(1))
        while starts != 0
            bit = trailing_zeros(starts)
            starts &= starts - 1
            candidate_start = block_start + bit
            candidate = prefix_hash_scan_candidate_sequence(
                source, candidate_start, candidate_bases)
            candidate === nothing && continue
            for (is_antisense, hits, strand_range) in (
                    (false, plus_hits, bounds.plus),
                    (true, minus_hits, bounds.minus))
                candidate_start in strand_range || continue
                pattern = is_antisense ? dbi.motif.rve : dbi.motif.fwd
                isempty(pattern) && continue
                matched, mask, ambiguous_prefix =
                    prefix_hash_scan_ambiguous_candidate_mask(
                        candidate, dbi.motif, is_antisense, hash_len, query)
                matched || continue
                if stats !== nothing
                    stats.motif_candidates += 1
                    stats.ambiguous_prefixes += ambiguous_prefix
                end
                mask == 0 || push!(
                    hits, PrefixHashScanHit(candidate_start, mask))
            end
        end
    end
    sort!(plus_hits; by = hit -> hit.start, alg = QuickSort)
    sort!(minus_hits; by = hit -> hit.start, alg = QuickSort)
    return nothing
end

# Whole-chromosome LongDNA scan for the :fused_directory fallback. Splits the
# candidate span across tasks, each running the geometry's range kernel.
function scan_prefix_hits(
    geometry::PrefixScanGeometry,
    chrom_seq::LongDNA{4},
    dbi::DBInfo,
    query,
    hash_len::Int,
    stats::Union{Nothing, PrefixHashScanStats} = nothing;
    scan_threads::Int = Threads.nthreads())

    bounds = prefix_scan_bounds(geometry, chrom_seq, dbi)
    bounds === nothing && return PrefixHashScanHit[], PrefixHashScanHit[]
    candidate_first, candidate_last = first(bounds.all), last(bounds.all)
    thread_count = min(max(scan_threads, 1), length(bounds.all))
    chunk_size = cld(length(bounds.all), thread_count)
    tasks = map(candidate_first:chunk_size:candidate_last) do first_
        part = PrefixScanBounds(
            first_:min(first_ + chunk_size - 1, candidate_last),
            bounds.plus, bounds.minus)
        Threads.@spawn scan_prefix_hits_range(
            geometry, chrom_seq, query, hash_len, part)
    end
    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    motif_candidates = 0
    for task in tasks
        local_plus, local_minus, local_count = fetch(task)
        append!(plus_hits, local_plus)
        append!(minus_hits, local_minus)
        motif_candidates += local_count
    end
    if dbi.motif.ambig_max > 0
        scan_ambiguous_prefix_hits_range!(
            plus_hits, minus_hits, chrom_seq, geometry, dbi, query,
            hash_len, bounds, Val(dbi.motif.ambig_max), stats)
    end
    stats === nothing || (stats.motif_candidates += motif_candidates)
    return plus_hits, minus_hits
end

@inline function prefix_hash_scan_record_candidate!(
    hits::Vector{PrefixHashScanHit},
    ::Nothing,
    query,
    candidate_start::Int,
    hash::Unsigned,
    ::Val{false})

    mask = prefix_hash_scan_candidate_mask(query, hash)
    mask != 0 && push!(hits, PrefixHashScanHit(candidate_start, mask))
    return nothing
end

@inline function prefix_hash_scan_record_candidate!(
    hits::Vector{PrefixHashScanHit},
    candidates::Vector{UInt64},
    query::PrefixHashScanPrefilteredDirectory,
    candidate_start::Int,
    hash::Unsigned,
    ::Val{true})

    prefix_hash_scan_prefilter_contains(query, hash) || return nothing
    packed = (UInt64(UInt32(hash)) << 32) | UInt64(UInt32(candidate_start))
    push!(candidates, packed)
    return nothing
end

function prefix_hash_scan_radix_pass!(
    destination::Vector{UInt64},
    source::Vector{UInt64},
    counts::Vector{Int},
    shift::Int,
    mask::UInt64)

    fill!(counts, 0)
    @inbounds for value in source
        counts[Int((value >> shift) & mask) + 1] += 1
    end
    next_idx = 1
    @inbounds for idx in eachindex(counts)
        count = counts[idx]
        counts[idx] = next_idx
        next_idx += count
    end
    @inbounds for value in source
        digit = Int((value >> shift) & mask) + 1
        destination[counts[digit]] = value
        counts[digit] += 1
    end
    return destination
end

function prefix_hash_scan_radix_order!(
    candidates::Vector{UInt64},
    scratch::Vector{UInt64},
    counts::Vector{Int})

    resize!(scratch, length(candidates))
    isempty(candidates) && return scratch
    prefix_hash_scan_radix_pass!(
        scratch, candidates, counts, 32, UInt64(0x03ff))
    prefix_hash_scan_radix_pass!(
        candidates, scratch, counts, 42, UInt64(0x07ff))
    prefix_hash_scan_radix_pass!(
        scratch, candidates, counts, 53, UInt64(0x07ff))
    return scratch
end

function resolve_prefix_hash_scan_bucketed_hits!(
    hits::Vector{PrefixHashScanHit},
    candidates::Vector{UInt64},
    scratch::Vector{UInt64},
    counts::Vector{Int},
    query::PrefixHashScanPrefilteredDirectory)

    ordered = prefix_hash_scan_radix_order!(candidates, scratch, counts)
    previous_hash = UInt32(0)
    candidate_mask = UInt64(0)
    have_previous = false
    @inbounds for packed in ordered
        hash = UInt32(packed >> 32)
        if !have_previous || hash != previous_hash
            previous_hash = hash
            candidate_mask = prefix_hash_scan_candidate_mask(
                query.directory, hash)
            have_previous = true
        end
        candidate_mask == 0 && continue
        candidate_start = Int(
            UInt32(packed & UInt64(typemax(UInt32))))
        push!(hits, PrefixHashScanHit(candidate_start, candidate_mask))
    end
    sort!(hits; by = hit -> hit.start, alg = QuickSort)
    return hits
end
