# Specialized scalar and AVX2/BMI2 Cas12a scan kernels.
# Cas12a geometry literals remain compile-time constants in these hot loops.

function scan_cas12a_prefix_hits_range(
    chrom_seq::LongDNA{4},
    query,
    hash_len::Int,
    bounds::PrefixScanBounds)

    hash_len == 16 || error("The specialized Cas12a scan requires hash_len=16.")
    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    candidate_first, candidate_last = first(bounds.all), last(bounds.all)
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

        if candidate_start in bounds.plus &&
                ((fwd_window >> 48) & 0x03) == 0x03 &&
                ((fwd_window >> 46) & 0x03) == 0x03 &&
                ((fwd_window >> 44) & 0x03) == 0x03 &&
                ((fwd_window >> 42) & 0x03) != 0x03
            motif_candidates += 1
            hash = UInt32((fwd_window >> 10) & UInt64(typemax(UInt32)))
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask != 0 && push!(plus_hits, PrefixHashScanHit(candidate_start, mask))
        end
        if candidate_start in bounds.minus &&
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
    bounds::PrefixScanBounds,
    ::Val{Bucketed},
    simd_backend::Val = default_prefix_hash_scan_simd_backend()) where Bucketed

    empty!(plus_hits)
    empty!(minus_hits)
    if Bucketed
        empty!(plus_candidates)
        empty!(minus_candidates)
    end
    motif_candidates = 0
    candidate_first, candidate_last = first(bounds.all), last(bounds.all)
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
            candidate_start in bounds.plus || continue
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
            candidate_start in bounds.minus || continue
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
        if candidate_start in bounds.plus &&
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
        if candidate_start in bounds.minus &&
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
    bounds::PrefixScanBounds,
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    return scan_cas12a_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, nothing, nothing, nothing, nothing, nothing,
        raw, query, bounds, Val(false), simd_backend)
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
    bounds::PrefixScanBounds,
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    return scan_cas12a_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, plus_candidates, minus_candidates,
        plus_radix_scratch, minus_radix_scratch, radix_counts, raw, query,
        bounds, Val(true), simd_backend)
end

function scan_cas12a_prefix_hits_raw_range(
    raw::AbstractVector{UInt8},
    query,
    bounds::PrefixScanBounds,
    simd_backend::Val = default_prefix_hash_scan_simd_backend())

    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    motif_candidates = scan_cas12a_prefix_hits_raw_range!(
        plus_hits, minus_hits, raw, query, bounds, simd_backend)
    return plus_hits, minus_hits, motif_candidates
end

scan_prefix_hits_range(
    ::PrefixScanGeometry{:cas12a}, chrom_seq, query, hash_len, bounds) =
    scan_cas12a_prefix_hits_range(chrom_seq, query, hash_len, bounds)

scan_prefix_hits_raw_range!(
    ::PrefixScanGeometry{:cas12a}, args...) =
    scan_cas12a_prefix_hits_raw_range!(args...)

scan_prefix_hits_raw_range_bucketed!(
    ::PrefixScanGeometry{:cas12a}, args...) =
    scan_cas12a_prefix_hits_raw_range_bucketed!(args...)

scan_prefix_hits_raw_range(
    ::PrefixScanGeometry{:cas12a}, args...) =
    scan_cas12a_prefix_hits_raw_range(args...)
