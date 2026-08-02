# Motif-specialized generic scanner. Motif decisions live in the matcher type,
# so generated block operations contain no runtime motif branches.

struct PrefixScanMatcher{Spec} end

@inline prefix_scan_matcher_spec(::PrefixScanMatcher{Spec}) where Spec = Spec

function prefix_scan_pattern_constraints(pattern::LongDNA{4})
    bases = (DNA_A, DNA_C, DNA_G, DNA_T)
    constraints = Tuple{Int, UInt8}[]
    for (idx, base) in enumerate(pattern)
        mask = UInt8(0)
        for (code, concrete) in enumerate(bases)
            iscompatible(base, concrete) &&
                (mask |= UInt8(1) << (code - 1))
        end
        mask == 0x0f || push!(constraints, (idx - 1, mask))
    end
    return Tuple(constraints)
end

function prefix_scan_normalized_offsets(
    span::Int, pam_loci::UnitRange{<:Integer},
    is_antisense::Bool, extends5::Bool)

    offsets = collect(0:(span - 1))
    isempty(pam_loci) || deleteat!(offsets, pam_loci)
    extends5 == is_antisense || reverse!(offsets)
    return Tuple(offsets)
end

function resolve_generic_prefix_scan_geometry(
    motif::Motif, distance::Int, hash_len::Int)

    distance in 0:4 || return nothing
    hash_len == 16 || return nothing
    guide_bases = length_noPAM(motif)
    span = length(motif)
    16 <= guide_bases <= 64 || return nothing
    guide_bases - distance >= hash_len || return nothing
    1 <= span <= 65 || return nothing
    fwd_enabled = !isempty(motif.fwd)
    rev_enabled = !isempty(motif.rve)
    (fwd_enabled || rev_enabled) || return nothing
    fwd_enabled && length(motif.fwd) != span && return nothing
    rev_enabled && length(motif.rve) != span && return nothing
    fwd_enabled && length(motif.fwd) - length(motif.pam_loci_fwd) != guide_bases &&
        return nothing
    rev_enabled && length(motif.rve) - length(motif.pam_loci_rve) != guide_bases &&
        return nothing

    fwd_offsets = fwd_enabled ? prefix_scan_normalized_offsets(
        span, motif.pam_loci_fwd, false, motif.extends5) : ()
    rev_offsets = rev_enabled ? prefix_scan_normalized_offsets(
        span, motif.pam_loci_rve, true, motif.extends5) : ()
    spec = (
        span = span,
        fwd_enabled = fwd_enabled,
        rev_enabled = rev_enabled,
        fwd_constraints = fwd_enabled ?
            prefix_scan_pattern_constraints(motif.fwd) : (),
        rev_constraints = rev_enabled ?
            prefix_scan_pattern_constraints(motif.rve) : (),
        fwd_offsets = fwd_offsets,
        rev_offsets = rev_offsets,
        fwd_step = motif.extends5 ? -1 : 1,
        rev_step = motif.extends5 ? 1 : -1,
        fwd_pos_offset = motif.extends5 ? span - 1 : 0,
        rev_pos_offset = motif.extends5 ? 0 : span - 1,
    )
    matcher = PrefixScanMatcher{spec}()
    return PrefixScanGeometry{:generic, typeof(matcher)}(
        guide_bases, span - guide_bases, hash_len, distance, matcher)
end

function candidate_prefix_hashes_direct(
    geometry::PrefixScanGeometry{:generic},
    chrom_seq::LongDNA{4}, candidate_range::UnitRange{Int64},
    is_antisense::Bool, hash_len::Int,
    hash_type::Type{<:Unsigned})

    hash_len == geometry.prefix_bases || return nothing
    spec = prefix_scan_matcher_spec(geometry.matcher)
    offsets = is_antisense ? spec.rev_offsets : spec.fwd_offsets
    hash = prefix_hash_scan_generic_hash_scalar(
        chrom_seq, first(candidate_range), offsets, is_antisense)
    hash === nothing && return nothing
    return hash_type[convert(hash_type, hash)]
end

@generated function prefix_hash_scan_generic_valid(
    exact::UInt128, ::PrefixScanMatcher{Spec}) where Spec

    span = Spec.span
    statements = Expr[]
    masks = Dict{Int, Any}(1 => :exact)
    width = 1
    while width * 2 <= span
        previous = masks[width]
        current = gensym(:valid)
        push!(statements, :($current = $previous & ($previous >> $width)))
        width *= 2
        masks[width] = current
    end
    offset = 0
    result = nothing
    power = width
    while power >= 1
        if span & power != 0
            term = offset == 0 ? masks[power] : :($(masks[power]) >> $offset)
            result = result === nothing ? term : :($result & $term)
            offset += power
        end
        power >>= 1
    end
    return quote
        $(statements...)
        UInt64(($result) & UInt128(typemax(UInt64)))
    end
end

@generated function prefix_hash_scan_generic_pattern_mask(
    a::UInt128, c::UInt128, g::UInt128, t::UInt128,
    ::PrefixScanMatcher{Spec}, ::Val{Anti}) where {Spec, Anti}

    enabled = Anti ? Spec.rev_enabled : Spec.fwd_enabled
    enabled || return :(UInt64(0))
    constraints = Anti ? Spec.rev_constraints : Spec.fwd_constraints
    result = :(typemax(UInt64))
    profiles = (:a, :c, :g, :t)
    for (offset, allowed) in constraints
        terms = Any[profiles[idx] for idx in 1:4 if allowed & (UInt8(1) << (idx - 1)) != 0]
        compatible = isempty(terms) ? :(UInt128(0)) : foldl(
            (left, right) -> :($left | $right), terms)
        shifted = :(UInt64(($compatible >> $offset) & UInt128(typemax(UInt64))))
        result = :($result & $shifted)
    end
    return result
end

@inline function prefix_hash_scan_reverse_lowbits(value::UInt64, ::Val{N}) where N
    return bitreverse(value) >> (64 - N)
end

function prefix_scan_hash_bits_expr(profile, offsets)
    runs = UnitRange{Int}[]
    first_idx = 1
    for idx in 2:length(offsets)
        if abs(offsets[idx] - offsets[idx - 1]) != 1
            push!(runs, first_idx:(idx - 1))
            first_idx = idx
        end
    end
    push!(runs, first_idx:length(offsets))

    result = :(UInt64(0))
    consumed = 0
    for run in runs
        run_offsets = offsets[run]
        count = length(run_offsets)
        minimum_offset = minimum(run_offsets)
        mask = (UInt64(1) << count) - UInt64(1)
        bits = :(UInt64(($profile >> (bit + $minimum_offset)) & UInt128($mask)))
        if count > 1 && run_offsets[2] < run_offsets[1]
            bits = :(prefix_hash_scan_reverse_lowbits($bits, Val($count)))
        end
        result = :($result | ($bits << $consumed))
        consumed += count
    end
    return result
end

@generated function prefix_hash_scan_generic_hash(
    low::UInt128, high::UInt128, bit::Int,
    ::PrefixScanMatcher{Spec}, ::Val{Anti}) where {Spec, Anti}

    offsets = collect((Anti ? Spec.rev_offsets : Spec.fwd_offsets)[1:16])
    low_bits = prefix_scan_hash_bits_expr(:low, offsets)
    high_bits = prefix_scan_hash_bits_expr(:high, offsets)
    complement_expr = Anti ? :(xor(hash, typemax(UInt32))) : :hash
    return quote
        low16 = $low_bits
        high16 = $high_bits
        hash = prefix_hash_scan_reverse_codes(
            prefix_hash_scan_pack_codes(low16, high16))
        $complement_expr
    end
end

@inline function prefix_hash_scan_generic_matches(
    raw::AbstractVector{UInt8}, candidate_start::Int, constraints)

    @inbounds for (offset, allowed) in constraints
        code = prefix_hash_scan_raw_code(raw[candidate_start + offset])
        code == 0xff && return false
        allowed & (UInt8(1) << code) != 0 || return false
    end
    return true
end

@inline function prefix_hash_scan_generic_matches(
    chrom_seq::LongDNA{4}, candidate_start::Int, constraints)

    @inbounds for (offset, allowed) in constraints
        code = prefix_hash_scan_twobit_nibble(UInt8(
            BioSequences.extract_encoded_element(chrom_seq, candidate_start + offset)))
        code == 0xff && return false
        allowed & (UInt8(1) << code) != 0 || return false
    end
    return true
end

@inline function prefix_hash_scan_generic_hash_scalar(
    source, candidate_start::Int, offsets, is_antisense::Bool)

    hash = UInt32(0)
    @inbounds for offset in offsets[1:16]
        code = source isa AbstractVector{UInt8} ?
            prefix_hash_scan_raw_code(source[candidate_start + offset]) :
            prefix_hash_scan_twobit_nibble(UInt8(
                BioSequences.extract_encoded_element(
                    source, candidate_start + offset)))
        code == 0xff && return nothing
        is_antisense && (code = UInt8(3) - code)
        hash = (hash << 2) | UInt32(code)
    end
    return hash
end

function prefix_scan_generic_bounds(
    n::Int, seq_start::Int, seq_stop::Int,
    geometry::PrefixScanGeometry{:generic}, dbi::DBInfo)

    n < prefix_scan_candidate_bases(geometry) && return nothing
    spec = prefix_scan_matcher_spec(geometry.matcher)
    last_offset = prefix_scan_candidate_last_offset(geometry)
    function strand_bounds(is_antisense, enabled)
        enabled || return (1, 0)
        if xor(is_antisense, dbi.motif.extends5)
            return max(seq_start - dbi.motif.distance, 1), seq_stop - last_offset
        end
        return seq_start, min(seq_stop + dbi.motif.distance, n) - last_offset
    end
    plus_first, plus_last = strand_bounds(false, spec.fwd_enabled)
    minus_first, minus_last = strand_bounds(true, spec.rev_enabled)
    firsts = Int[]
    lasts = Int[]
    plus_first <= plus_last && (push!(firsts, plus_first); push!(lasts, plus_last))
    minus_first <= minus_last && (push!(firsts, minus_first); push!(lasts, minus_last))
    isempty(firsts) && return nothing
    return minimum(firsts), maximum(lasts), plus_first, plus_last, minus_first, minus_last
end

function generic_prefix_scan_bounds(raw::AbstractVector{UInt8}, geometry, dbi)
    seq_start = 1
    seq_stop = length(raw)
    @inbounds while seq_start <= seq_stop &&
            (raw[seq_start] == UInt8('N') || raw[seq_start] == UInt8('n'))
        seq_start += 1
    end
    @inbounds while seq_stop > 0 &&
            (raw[seq_stop] == UInt8('N') || raw[seq_stop] == UInt8('n'))
        seq_stop -= 1
    end
    return prefix_scan_generic_bounds(
        length(raw), seq_start, seq_stop, geometry, dbi)
end

function generic_prefix_scan_bounds(chrom_seq::LongDNA{4}, geometry, dbi)
    seq_start, seq_stop = locate_telomeres(chrom_seq)
    return prefix_scan_generic_bounds(
        length(chrom_seq), seq_start, seq_stop, geometry, dbi)
end

function scan_generic_prefix_hits_raw_range_impl!(
    plus_hits::Vector{PrefixHashScanHit},
    minus_hits::Vector{PrefixHashScanHit},
    plus_candidates, minus_candidates, plus_radix_scratch,
    minus_radix_scratch, radix_counts,
    raw::AbstractVector{UInt8}, query,
    geometry::PrefixScanGeometry{:generic},
    candidate_first::Int, candidate_last::Int,
    plus_first::Int, plus_last::Int,
    minus_first::Int, minus_last::Int,
    ::Val{Bucketed}) where Bucketed

    empty!(plus_hits)
    empty!(minus_hits)
    if Bucketed
        empty!(plus_candidates)
        empty!(minus_candidates)
    end
    candidate_first > candidate_last && return 0
    matcher = geometry.matcher
    spec = prefix_scan_matcher_spec(matcher)
    motif_candidates = 0
    block_start = candidate_first
    n = length(raw)
    if block_start + 127 <= n && block_start + 63 <= candidate_last
        a0, c0, g0, t0 = prefix_hash_scan_raw_profile64(raw, block_start)
    end
    while block_start + 127 <= n && block_start + 63 <= candidate_last
        a1, c1, g1, t1 = prefix_hash_scan_raw_profile64(raw, block_start + 64)
        a = UInt128(a0) | (UInt128(a1) << 64)
        c = UInt128(c0) | (UInt128(c1) << 64)
        g = UInt128(g0) | (UInt128(g1) << 64)
        t = UInt128(t0) | (UInt128(t1) << 64)
        valid = prefix_hash_scan_generic_valid(a | c | g | t, matcher)
        count = min(64, candidate_last - block_start + 1)
        count < 64 && (valid &= (UInt64(1) << count) - UInt64(1))
        plus_mask = valid & prefix_hash_scan_generic_pattern_mask(
            a, c, g, t, matcher, Val(false))
        minus_mask = valid & prefix_hash_scan_generic_pattern_mask(
            a, c, g, t, matcher, Val(true))
        low = c | t
        high = g | t
        while plus_mask != 0
            bit = trailing_zeros(plus_mask)
            plus_mask &= plus_mask - 1
            candidate_start = block_start + bit
            plus_first <= candidate_start <= plus_last || continue
            motif_candidates += 1
            hash = prefix_hash_scan_generic_hash(
                low, high, bit, matcher, Val(false))
            prefix_hash_scan_record_candidate!(
                plus_hits, plus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
        while minus_mask != 0
            bit = trailing_zeros(minus_mask)
            minus_mask &= minus_mask - 1
            candidate_start = block_start + bit
            minus_first <= candidate_start <= minus_last || continue
            motif_candidates += 1
            hash = prefix_hash_scan_generic_hash(
                low, high, bit, matcher, Val(true))
            prefix_hash_scan_record_candidate!(
                minus_hits, minus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
        block_start += 64
        a0, c0, g0, t0 = a1, c1, g1, t1
    end

    @inbounds for candidate_start in block_start:candidate_last
        valid = true
        for offset in 0:(spec.span - 1)
            if prefix_hash_scan_raw_code(raw[candidate_start + offset]) == 0xff
                valid = false
                break
            end
        end
        valid || continue
        if plus_first <= candidate_start <= plus_last && spec.fwd_enabled &&
                prefix_hash_scan_generic_matches(
                    raw, candidate_start, spec.fwd_constraints)
            motif_candidates += 1
            hash = prefix_hash_scan_generic_hash_scalar(
                raw, candidate_start, spec.fwd_offsets, false)
            prefix_hash_scan_record_candidate!(
                plus_hits, plus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
        if minus_first <= candidate_start <= minus_last && spec.rev_enabled &&
                prefix_hash_scan_generic_matches(
                    raw, candidate_start, spec.rev_constraints)
            motif_candidates += 1
            hash = prefix_hash_scan_generic_hash_scalar(
                raw, candidate_start, spec.rev_offsets, true)
            prefix_hash_scan_record_candidate!(
                minus_hits, minus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
    end
    if Bucketed
        resolve_prefix_hash_scan_bucketed_hits!(
            plus_hits, plus_candidates, plus_radix_scratch, radix_counts, query)
        resolve_prefix_hash_scan_bucketed_hits!(
            minus_hits, minus_candidates, minus_radix_scratch, radix_counts, query)
    end
    return motif_candidates
end

function scan_generic_prefix_hits_raw_range!(
    plus_hits, minus_hits, raw, query, geometry, bounds...)
    return scan_generic_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, nothing, nothing, nothing, nothing, nothing,
        raw, query, geometry, bounds..., Val(false))
end

function scan_generic_prefix_hits_raw_range_bucketed!(
    plus_hits, minus_hits, plus_candidates, minus_candidates,
    plus_radix_scratch, minus_radix_scratch, radix_counts,
    raw, query, geometry, bounds...)
    return scan_generic_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, plus_candidates, minus_candidates,
        plus_radix_scratch, minus_radix_scratch, radix_counts,
        raw, query, geometry, bounds..., Val(true))
end

function scan_generic_prefix_hits_raw_range(raw, query, geometry, bounds...)
    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    count = scan_generic_prefix_hits_raw_range!(
        plus_hits, minus_hits, raw, query, geometry, bounds...)
    return plus_hits, minus_hits, count
end

function scan_generic_prefix_hits_range(
    chrom_seq::LongDNA{4}, query, geometry::PrefixScanGeometry{:generic},
    candidate_first::Int, candidate_last::Int,
    plus_first::Int, plus_last::Int,
    minus_first::Int, minus_last::Int)

    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    spec = prefix_scan_matcher_spec(geometry.matcher)
    motif_candidates = 0
    @inbounds for candidate_start in candidate_first:candidate_last
        valid = true
        for offset in 0:(spec.span - 1)
            code = prefix_hash_scan_twobit_nibble(UInt8(
                BioSequences.extract_encoded_element(
                    chrom_seq, candidate_start + offset)))
            if code == 0xff
                valid = false
                break
            end
        end
        valid || continue
        if plus_first <= candidate_start <= plus_last && spec.fwd_enabled &&
                prefix_hash_scan_generic_matches(
                    chrom_seq, candidate_start, spec.fwd_constraints)
            motif_candidates += 1
            hash = prefix_hash_scan_generic_hash_scalar(
                chrom_seq, candidate_start, spec.fwd_offsets, false)
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || push!(plus_hits, PrefixHashScanHit(candidate_start, mask))
        end
        if minus_first <= candidate_start <= minus_last && spec.rev_enabled &&
                prefix_hash_scan_generic_matches(
                    chrom_seq, candidate_start, spec.rev_constraints)
            motif_candidates += 1
            hash = prefix_hash_scan_generic_hash_scalar(
                chrom_seq, candidate_start, spec.rev_offsets, true)
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || push!(minus_hits, PrefixHashScanHit(candidate_start, mask))
        end
    end
    return plus_hits, minus_hits, motif_candidates
end

function scan_generic_prefix_hits_raw(
    raw, dbi, query, geometry::PrefixScanGeometry{:generic},
    stats = nothing; scan_threads::Int = Threads.nthreads())

    bounds = generic_prefix_scan_bounds(raw, geometry, dbi)
    bounds === nothing && return PrefixHashScanHit[], PrefixHashScanHit[]
    candidate_first, candidate_last, plus_first, plus_last, minus_first, minus_last = bounds
    candidate_count = candidate_last - candidate_first + 1
    thread_count = min(max(scan_threads, 1), candidate_count)
    chunk_size = cld(candidate_count, thread_count)
    tasks = map(candidate_first:chunk_size:candidate_last) do first_
        last_ = min(first_ + chunk_size - 1, candidate_last)
        Threads.@spawn scan_generic_prefix_hits_raw_range(
            raw, query, geometry, first_, last_, plus_first, plus_last,
            minus_first, minus_last)
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
            plus_hits, minus_hits, raw, geometry, dbi, query,
            geometry.prefix_bases, candidate_first, candidate_last,
            plus_first, plus_last, minus_first, minus_last,
            Val(dbi.motif.ambig_max), stats)
    end
    stats === nothing || (stats.motif_candidates += motif_candidates)
    return plus_hits, minus_hits
end

function scan_generic_prefix_hits(
    chrom_seq::LongDNA{4}, dbi, query, geometry::PrefixScanGeometry{:generic},
    stats = nothing; scan_threads::Int = Threads.nthreads())

    bounds = generic_prefix_scan_bounds(chrom_seq, geometry, dbi)
    bounds === nothing && return PrefixHashScanHit[], PrefixHashScanHit[]
    candidate_first, candidate_last, plus_first, plus_last, minus_first, minus_last = bounds
    candidate_count = candidate_last - candidate_first + 1
    thread_count = min(max(scan_threads, 1), candidate_count)
    chunk_size = cld(candidate_count, thread_count)
    tasks = map(candidate_first:chunk_size:candidate_last) do first_
        last_ = min(first_ + chunk_size - 1, candidate_last)
        Threads.@spawn scan_generic_prefix_hits_range(
            chrom_seq, query, geometry, first_, last_, plus_first, plus_last,
            minus_first, minus_last)
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
            geometry.prefix_bases, candidate_first, candidate_last,
            plus_first, plus_last, minus_first, minus_last,
            Val(dbi.motif.ambig_max), stats)
    end
    stats === nothing || (stats.motif_candidates += motif_candidates)
    return plus_hits, minus_hits
end

scan_prefix_hits(
    geometry::PrefixScanGeometry{:generic}, chrom_seq, dbi, query, hash_len,
    stats = nothing; kwargs...) =
    scan_generic_prefix_hits(
        chrom_seq, dbi, query, geometry, stats; kwargs...)

scan_prefix_hits_raw(
    geometry::PrefixScanGeometry{:generic}, raw, dbi, query, stats = nothing;
    kwargs...) =
    scan_generic_prefix_hits_raw(raw, dbi, query, geometry, stats; kwargs...)

scan_prefix_hits_raw_range!(
    geometry::PrefixScanGeometry{:generic}, plus_hits, minus_hits,
    raw, query, bounds...) =
    scan_generic_prefix_hits_raw_range!(
        plus_hits, minus_hits, raw, query, geometry, bounds...)

scan_prefix_hits_raw_range_bucketed!(
    geometry::PrefixScanGeometry{:generic}, plus_hits, minus_hits,
    plus_candidates, minus_candidates, plus_radix_scratch,
    minus_radix_scratch, radix_counts, raw, query, bounds...) =
    scan_generic_prefix_hits_raw_range_bucketed!(
        plus_hits, minus_hits, plus_candidates, minus_candidates,
        plus_radix_scratch, minus_radix_scratch, radix_counts,
        raw, query, geometry, bounds...)

scan_prefix_hits_raw_range(
    geometry::PrefixScanGeometry{:generic}, raw, query, bounds...) =
    scan_generic_prefix_hits_raw_range(raw, query, geometry, bounds...)

function scan_verify_prefix_raw_range!(
    geometry::PrefixScanGeometry{:generic}, plus, minus, raw, query,
    candidate_first, candidate_last, plus_first, plus_last,
    minus_first, minus_last, global_offset, dbi, guides_, myers_profiles,
    distance, stats)

    plus_hits, minus_hits, motif_candidates = scan_generic_prefix_hits_raw_range(
        raw, query, geometry, candidate_first, candidate_last,
        plus_first, plus_last, minus_first, minus_last)
    stats === nothing || (stats.motif_candidates += motif_candidates)
    evaluate_prefix_hash_scan_hits!(
        plus, raw, geometry, plus_hits, global_offset, dbi, false,
        guides_, myers_profiles, distance, stats)
    evaluate_prefix_hash_scan_hits!(
        minus, raw, geometry, minus_hits, global_offset, dbi, true,
        guides_, myers_profiles, distance, stats)
    return plus, minus
end
