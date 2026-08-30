# Specialized scalar and AVX2/BMI2 Cas12a scan kernels.
# Cas12a geometry literals remain compile-time constants in these hot loops.

function cas12a_prefix_scan_bounds(chrom_seq::LongDNA{4}, dbi::DBInfo)
    n = length(chrom_seq)
    n < 25 && return nothing
    seq_start, seq_stop = locate_telomeres(chrom_seq)
    plus_first = seq_start
    plus_last = min(seq_stop + dbi.motif.distance, n) - 24
    minus_first = max(seq_start - dbi.motif.distance, 1)
    minus_last = seq_stop - 24
    firsts = Int[]
    lasts = Int[]
    plus_first <= plus_last && (push!(firsts, plus_first); push!(lasts, plus_last))
    minus_first <= minus_last && (push!(firsts, minus_first); push!(lasts, minus_last))
    isempty(firsts) && return nothing
    return minimum(firsts), maximum(lasts), plus_first, plus_last, minus_first, minus_last
end

function scan_cas12a_prefix_hits_range(
    chrom_seq::LongDNA{4},
    query,
    hash_len::Int,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int)

    hash_len == 16 || error("The specialized Cas12a scan requires hash_len=16.")
    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    candidate_first > candidate_last && return plus_hits, minus_hits, 0

    window_mask = (UInt64(1) << 50) - UInt64(1)
    fwd_window = zero(UInt64)
    valid_run = 0
    motif_candidates = 0

    @inbounds for pos in candidate_first:(candidate_first + 23)
        code = prefix_hash_scan_twobit_nibble(
            UInt8(BioSequences.extract_encoded_element(chrom_seq, pos)))
        if code == 0xff
            valid_run = 0
            code = 0x00
        else
            valid_run += 1
        end
        fwd_window = ((fwd_window << 2) | UInt64(code)) & window_mask
    end

    @inbounds for pos in (candidate_first + 24):(candidate_last + 24)
        code = prefix_hash_scan_twobit_nibble(
            UInt8(BioSequences.extract_encoded_element(chrom_seq, pos)))
        if code == 0xff
            valid_run = 0
            code = 0x00
        else
            valid_run += 1
        end
        fwd_window = ((fwd_window << 2) | UInt64(code)) & window_mask
        candidate_start = pos - 24
        valid_run >= 25 || continue

        if plus_first <= candidate_start <= plus_last &&
                ((fwd_window >> 48) & 0x03) == 0x03 &&
                ((fwd_window >> 46) & 0x03) == 0x03 &&
                ((fwd_window >> 44) & 0x03) == 0x03 &&
                ((fwd_window >> 42) & 0x03) != 0x03
            motif_candidates += 1
            hash = UInt32((fwd_window >> 10) & UInt64(typemax(UInt32)))
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask != 0 && push!(plus_hits, PrefixHashScanHit(candidate_start, mask))
        end
        if minus_first <= candidate_start <= minus_last &&
                ((fwd_window >> 6) & 0x03) != 0x00 &&
                ((fwd_window >> 4) & 0x03) == 0x00 &&
                ((fwd_window >> 2) & 0x03) == 0x00 &&
                (fwd_window & 0x03) == 0x00
            motif_candidates += 1
            hash = xor(
                prefix_hash_scan_reverse_codes(
                    UInt32((fwd_window >> 8) & UInt64(typemax(UInt32)))),
                typemax(UInt32),
            )
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask != 0 && push!(minus_hits, PrefixHashScanHit(candidate_start, mask))
        end
    end
    return plus_hits, minus_hits, motif_candidates
end

@inline function prefix_hash_scan_valid25(exact::UInt128)
    valid2 = exact & (exact >> 1)
    valid4 = valid2 & (valid2 >> 2)
    valid8 = valid4 & (valid4 >> 4)
    valid16 = valid8 & (valid8 >> 8)
    return valid16 & (valid8 >> 16) & (exact >> 24)
end

@inline function prefix_hash_scan_raw_hash_cas12a(
    raw::AbstractVector{UInt8}, candidate_start::Int, is_antisense::Bool)

    hash = UInt32(0)
    positions = is_antisense ?
        ((candidate_start + 20):-1:(candidate_start + 5)) :
        ((candidate_start + 4):(candidate_start + 19))
    @inbounds for pos in positions
        code = prefix_hash_scan_raw_code(raw[pos])
        code == 0xff && return nothing
        is_antisense && (code = UInt8(3) - code)
        hash = (hash << 2) | UInt32(code)
    end
    return hash
end

function cas12a_prefix_scan_bounds_raw(raw::AbstractVector{UInt8}, dbi::DBInfo)
    n = length(raw)
    n < 25 && return nothing
    seq_start = 1
    seq_stop = n
    @inbounds while seq_start <= seq_stop &&
            (raw[seq_start] == UInt8('N') || raw[seq_start] == UInt8('n'))
        seq_start += 1
    end
    @inbounds while seq_stop > 0 &&
            (raw[seq_stop] == UInt8('N') || raw[seq_stop] == UInt8('n'))
        seq_stop -= 1
    end
    plus_first = seq_start
    plus_last = min(seq_stop + dbi.motif.distance, n) - 24
    minus_first = max(seq_start - dbi.motif.distance, 1)
    minus_last = seq_stop - 24
    firsts = Int[]
    lasts = Int[]
    plus_first <= plus_last && (push!(firsts, plus_first); push!(lasts, plus_last))
    minus_first <= minus_last && (push!(firsts, minus_first); push!(lasts, minus_last))
    isempty(firsts) && return nothing
    return minimum(firsts), maximum(lasts), plus_first, plus_last, minus_first, minus_last
end

function scan_cas12a_prefix_hits_raw_range_impl!(
    plus_hits::Vector{PrefixHashScanHit},
    minus_hits::Vector{PrefixHashScanHit},
    plus_candidates,
    minus_candidates,
    plus_radix_scratch,
    minus_radix_scratch,
    radix_counts,
    raw::AbstractVector{UInt8},
    query,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int,
    ::Val{Bucketed},
    simd_backend::Val = default_prefix_hash_scan_simd_backend()) where Bucketed

    empty!(plus_hits)
    empty!(minus_hits)
    if Bucketed
        empty!(plus_candidates)
        empty!(minus_candidates)
    end
    motif_candidates = 0
    candidate_first > candidate_last && return motif_candidates
    n = length(raw)
    block_start = candidate_first

    if block_start + 127 <= n && block_start + 63 <= candidate_last
        a0, c0, g0, t0 = prefix_hash_scan_raw_profile64(
            raw, block_start, simd_backend)
    end
    while block_start + 127 <= n && block_start + 63 <= candidate_last
        a1, c1, g1, t1 = prefix_hash_scan_raw_profile64(
            raw, block_start + 64, simd_backend)
        a = UInt128(a0) | (UInt128(a1) << 64)
        c = UInt128(c0) | (UInt128(c1) << 64)
        g = UInt128(g0) | (UInt128(g1) << 64)
        t = UInt128(t0) | (UInt128(t1) << 64)
        exact = a | c | g | t
        valid = UInt64(prefix_hash_scan_valid25(exact) & UInt128(typemax(UInt64)))
        count = min(64, candidate_last - block_start + 1)
        count_mask = count == 64 ? typemax(UInt64) : (UInt64(1) << count) - 1
        valid &= count_mask
        plus_mask = valid &
            UInt64(t & UInt128(typemax(UInt64))) &
            UInt64((t >> 1) & UInt128(typemax(UInt64))) &
            UInt64((t >> 2) & UInt128(typemax(UInt64))) &
            UInt64(((a | c | g) >> 3) & UInt128(typemax(UInt64)))
        minus_mask = valid &
            UInt64(((c | g | t) >> 21) & UInt128(typemax(UInt64))) &
            UInt64((a >> 22) & UInt128(typemax(UInt64))) &
            UInt64((a >> 23) & UInt128(typemax(UInt64))) &
            UInt64((a >> 24) & UInt128(typemax(UInt64)))
        low = c | t
        high = g | t

        while plus_mask != 0
            bit = trailing_zeros(plus_mask)
            plus_mask &= plus_mask - 1
            candidate_start = block_start + bit
            plus_first <= candidate_start <= plus_last || continue
            motif_candidates += 1
            low16 = UInt64((low >> (bit + 4)) & UInt128(0xffff))
            high16 = UInt64((high >> (bit + 4)) & UInt128(0xffff))
            hash = prefix_hash_scan_reverse_codes(
                prefix_hash_scan_pack_codes(low16, high16))
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
            low16 = UInt64((low >> (bit + 5)) & UInt128(0xffff))
            high16 = UInt64((high >> (bit + 5)) & UInt128(0xffff))
            hash = xor(
                prefix_hash_scan_pack_codes(low16, high16), typemax(UInt32))
            prefix_hash_scan_record_candidate!(
                minus_hits, minus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
        block_start += 64
        a0, c0, g0, t0 = a1, c1, g1, t1
    end

    @inbounds for candidate_start in block_start:candidate_last
        valid = true
        for pos in candidate_start:(candidate_start + 24)
            if prefix_hash_scan_raw_code(raw[pos]) == 0xff
                valid = false
                break
            end
        end
        valid || continue
        if plus_first <= candidate_start <= plus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start]) == 3 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 1]) == 3 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 2]) == 3 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 3]) != 3
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash_cas12a(raw, candidate_start, false)
            prefix_hash_scan_record_candidate!(
                plus_hits, plus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
        if minus_first <= candidate_start <= minus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start + 21]) != 0 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 22]) == 0 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 23]) == 0 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 24]) == 0
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash_cas12a(raw, candidate_start, true)
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

function scan_cas12a_prefix_hits_raw_range!(
    plus_hits::Vector{PrefixHashScanHit},
    minus_hits::Vector{PrefixHashScanHit},
    raw::AbstractVector{UInt8},
    query,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int,
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    return scan_cas12a_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, nothing, nothing, nothing, nothing, nothing,
        raw, query, candidate_first, candidate_last, plus_first, plus_last,
        minus_first, minus_last, Val(false), simd_backend)
end

function scan_cas12a_prefix_hits_raw_range_bucketed!(
    plus_hits::Vector{PrefixHashScanHit},
    minus_hits::Vector{PrefixHashScanHit},
    plus_candidates::Vector{UInt64},
    minus_candidates::Vector{UInt64},
    plus_radix_scratch::Vector{UInt64},
    minus_radix_scratch::Vector{UInt64},
    radix_counts::Vector{Int},
    raw::AbstractVector{UInt8},
    query::PrefixHashScanPrefilteredDirectory,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int,
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    return scan_cas12a_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, plus_candidates, minus_candidates,
        plus_radix_scratch, minus_radix_scratch, radix_counts, raw, query,
        candidate_first, candidate_last, plus_first, plus_last, minus_first,
        minus_last, Val(true), simd_backend)
end

function scan_cas12a_prefix_hits_raw_range(
    raw::AbstractVector{UInt8},
    query,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int,
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    motif_candidates = scan_cas12a_prefix_hits_raw_range!(
        plus_hits, minus_hits, raw, query, candidate_first, candidate_last,
        plus_first, plus_last, minus_first, minus_last, simd_backend)
    return plus_hits, minus_hits, motif_candidates
end

function scan_cas12a_prefix_hits_raw(
    raw::AbstractVector{UInt8},
    dbi::DBInfo,
    query,
    stats::Union{Nothing, PrefixHashScanStats} = nothing;
    scan_threads::Int = Threads.nthreads(),
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    bounds = cas12a_prefix_scan_bounds_raw(raw, dbi)
    bounds === nothing && return PrefixHashScanHit[], PrefixHashScanHit[]
    candidate_first, candidate_last, plus_first, plus_last, minus_first, minus_last = bounds
    candidate_count = candidate_last - candidate_first + 1
    thread_count = min(max(scan_threads, 1), candidate_count)
    chunk_size = cld(candidate_count, thread_count)
    ranges = [
        first:min(first + chunk_size - 1, candidate_last)
        for first in candidate_first:chunk_size:candidate_last
    ]
    tasks = map(ranges) do range
        Threads.@spawn scan_cas12a_prefix_hits_raw_range(
            raw, query, first(range), last(range),
            plus_first, plus_last, minus_first, minus_last, simd_backend)
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
    if dbi.motif.ambig_max > 0
        scan_ambiguous_prefix_hits_range!(
            plus_hits, minus_hits, raw, CAS12A_D3_PREFIX_SCAN_GEOMETRY,
            dbi, query, 16, candidate_first, candidate_last,
            plus_first, plus_last, minus_first, minus_last,
            Val(dbi.motif.ambig_max), stats, simd_backend)
    end
    stats !== nothing && (stats.motif_candidates += motif_candidates)
    return plus_hits, minus_hits
end

function scan_cas12a_prefix_hits(
    chrom_seq::LongDNA{4},
    dbi::DBInfo,
    query,
    hash_len::Int,
    stats::Union{Nothing, PrefixHashScanStats} = nothing;
    scan_threads::Int = Threads.nthreads())

    bounds = cas12a_prefix_scan_bounds(chrom_seq, dbi)
    bounds === nothing && return PrefixHashScanHit[], PrefixHashScanHit[]
    candidate_first, candidate_last, plus_first, plus_last, minus_first, minus_last = bounds
    candidate_count = candidate_last - candidate_first + 1
    thread_count = min(max(scan_threads, 1), candidate_count)
    if thread_count == 1
        plus_hits, minus_hits, motif_candidates = scan_cas12a_prefix_hits_range(
            chrom_seq, query, hash_len, candidate_first, candidate_last,
            plus_first, plus_last, minus_first, minus_last)
    else
        chunk_size = cld(candidate_count, thread_count)
        ranges = [
            first:min(first + chunk_size - 1, candidate_last)
            for first in candidate_first:chunk_size:candidate_last
        ]
        tasks = map(ranges) do range
            Threads.@spawn scan_cas12a_prefix_hits_range(
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
    if dbi.motif.ambig_max > 0
        scan_ambiguous_prefix_hits_range!(
            plus_hits, minus_hits, chrom_seq, CAS12A_D3_PREFIX_SCAN_GEOMETRY,
            dbi, query, hash_len, candidate_first, candidate_last,
            plus_first, plus_last, minus_first, minus_last,
            Val(dbi.motif.ambig_max), stats)
    end
    stats !== nothing && (stats.motif_candidates += motif_candidates)
    return plus_hits, minus_hits
end

function scan_verify_cas12a_prefix_raw_range!(
    plus::Vector{PrefixHashScanVerifiedHit},
    minus::Vector{PrefixHashScanVerifiedHit},
    raw::AbstractVector{UInt8},
    query,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int,
    global_offset::Int,
    dbi::DBInfo,
    guides_::Vector{LongDNA{4}},
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    stats::S,
    simd_backend::Val = default_prefix_hash_scan_simd_backend()) where {S <: Union{Nothing, PrefixHashScanStats}}

    motif_candidates = 0
    candidate_first > candidate_last && return plus, minus
    geometry = PrefixScanGeometry{:cas12a}(21, 4, 16, distance)
    n = length(raw)
    block_start = candidate_first

    if block_start + 127 <= n && block_start + 63 <= candidate_last
        a0, c0, g0, t0 = prefix_hash_scan_raw_profile64(
            raw, block_start, simd_backend)
    end
    while block_start + 127 <= n && block_start + 63 <= candidate_last
        a1, c1, g1, t1 = prefix_hash_scan_raw_profile64(
            raw, block_start + 64, simd_backend)
        a = UInt128(a0) | (UInt128(a1) << 64)
        c = UInt128(c0) | (UInt128(c1) << 64)
        g = UInt128(g0) | (UInt128(g1) << 64)
        t = UInt128(t0) | (UInt128(t1) << 64)
        exact = a | c | g | t
        valid = UInt64(prefix_hash_scan_valid25(exact) & UInt128(typemax(UInt64)))
        count = min(64, candidate_last - block_start + 1)
        count_mask = count == 64 ? typemax(UInt64) : (UInt64(1) << count) - 1
        valid &= count_mask
        plus_mask = valid &
            UInt64(t & UInt128(typemax(UInt64))) &
            UInt64((t >> 1) & UInt128(typemax(UInt64))) &
            UInt64((t >> 2) & UInt128(typemax(UInt64))) &
            UInt64(((a | c | g) >> 3) & UInt128(typemax(UInt64)))
        minus_mask = valid &
            UInt64(((c | g | t) >> 21) & UInt128(typemax(UInt64))) &
            UInt64((a >> 22) & UInt128(typemax(UInt64))) &
            UInt64((a >> 23) & UInt128(typemax(UInt64))) &
            UInt64((a >> 24) & UInt128(typemax(UInt64)))
        low = c | t
        high = g | t

        while plus_mask != 0
            bit = trailing_zeros(plus_mask)
            plus_mask &= plus_mask - 1
            candidate_start = block_start + bit
            plus_first <= candidate_start <= plus_last || continue
            motif_candidates += 1
            low16 = UInt64((low >> (bit + 4)) & UInt128(0xffff))
            high16 = UInt64((high >> (bit + 4)) & UInt128(0xffff))
            hash = prefix_hash_scan_reverse_codes(
                prefix_hash_scan_pack_codes(low16, high16))
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                plus, raw, geometry,
                candidate_start, mask, global_offset, dbi, false,
                guides_, myers_profiles, distance, stats)
        end

        while minus_mask != 0
            bit = trailing_zeros(minus_mask)
            minus_mask &= minus_mask - 1
            candidate_start = block_start + bit
            minus_first <= candidate_start <= minus_last || continue
            motif_candidates += 1
            low16 = UInt64((low >> (bit + 5)) & UInt128(0xffff))
            high16 = UInt64((high >> (bit + 5)) & UInt128(0xffff))
            hash = xor(
                prefix_hash_scan_pack_codes(low16, high16), typemax(UInt32))
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                minus, raw, geometry,
                candidate_start, mask, global_offset, dbi, true,
                guides_, myers_profiles, distance, stats)
        end
        block_start += 64
        a0, c0, g0, t0 = a1, c1, g1, t1
    end

    @inbounds for candidate_start in block_start:candidate_last
        valid = true
        for pos in candidate_start:(candidate_start + 24)
            if prefix_hash_scan_raw_code(raw[pos]) == 0xff
                valid = false
                break
            end
        end
        valid || continue
        if plus_first <= candidate_start <= plus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start]) == 3 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 1]) == 3 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 2]) == 3 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 3]) != 3
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash_cas12a(raw, candidate_start, false)
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                plus, raw, geometry,
                candidate_start, mask, global_offset, dbi, false,
                guides_, myers_profiles, distance, stats)
        end
        if minus_first <= candidate_start <= minus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start + 21]) != 0 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 22]) == 0 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 23]) == 0 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 24]) == 0
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash_cas12a(raw, candidate_start, true)
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                minus, raw, geometry,
                candidate_start, mask, global_offset, dbi, true,
                guides_, myers_profiles, distance, stats)
        end
    end
    stats !== nothing && (stats.motif_candidates += motif_candidates)
    return plus, minus
end

scan_prefix_hits(
    ::PrefixScanGeometry{:cas9}, args...; kwargs...) =
    scan_cas9_prefix_hits(args...; kwargs...)
scan_prefix_hits(
    ::PrefixScanGeometry{:cas12a}, args...; kwargs...) =
    scan_cas12a_prefix_hits(args...; kwargs...)

scan_prefix_hits_raw(
    ::PrefixScanGeometry{:cas9}, args...; kwargs...) =
    scan_cas9_prefix_hits_raw(args...; kwargs...)
scan_prefix_hits_raw(
    ::PrefixScanGeometry{:cas12a}, args...; kwargs...) =
    scan_cas12a_prefix_hits_raw(args...; kwargs...)

scan_prefix_hits_raw_range!(
    ::PrefixScanGeometry{:cas9}, args...) =
    scan_cas9_prefix_hits_raw_range!(args...)
scan_prefix_hits_raw_range!(
    ::PrefixScanGeometry{:cas12a}, args...) =
    scan_cas12a_prefix_hits_raw_range!(args...)

scan_prefix_hits_raw_range_bucketed!(
    ::PrefixScanGeometry{:cas9}, args...) =
    scan_cas9_prefix_hits_raw_range_bucketed!(args...)
scan_prefix_hits_raw_range_bucketed!(
    ::PrefixScanGeometry{:cas12a}, args...) =
    scan_cas12a_prefix_hits_raw_range_bucketed!(args...)

scan_prefix_hits_raw_range(
    ::PrefixScanGeometry{:cas9}, args...) =
    scan_cas9_prefix_hits_raw_range(args...)
scan_prefix_hits_raw_range(
    ::PrefixScanGeometry{:cas12a}, args...) =
    scan_cas12a_prefix_hits_raw_range(args...)

scan_verify_prefix_raw_range!(
    ::PrefixScanGeometry{:cas9}, args...) =
    scan_verify_cas9_prefix_raw_range!(args...)
scan_verify_prefix_raw_range!(
    ::PrefixScanGeometry{:cas12a}, args...) =
    scan_verify_cas12a_prefix_raw_range!(args...)
