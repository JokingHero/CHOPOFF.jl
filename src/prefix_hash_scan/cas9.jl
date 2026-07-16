# Specialized scalar and AVX2/BMI2 Cas9 scan kernels.
# Cas9 geometry literals remain compile-time constants in these hot loops.

@inline function prefix_hash_scan_twobit_nibble(nibble::UInt8)
    nibble == 0x01 && return UInt8(0)
    nibble == 0x02 && return UInt8(1)
    nibble == 0x04 && return UInt8(2)
    nibble == 0x08 && return UInt8(3)
    return UInt8(0xff)
end

function cas9_prefix_scan_bounds(chrom_seq::LongDNA{4}, dbi::DBInfo)
    n = length(chrom_seq)
    geometry = CAS9_D3_PREFIX_SCAN_GEOMETRY
    candidate_bases = prefix_scan_candidate_bases(geometry)
    candidate_last_offset = prefix_scan_candidate_last_offset(geometry)
    n < candidate_bases && return nothing
    seq_start, seq_stop = locate_telomeres(chrom_seq)
    plus_first = max(seq_start - dbi.motif.distance, 1)
    plus_last = seq_stop - candidate_last_offset
    minus_first = seq_start
    minus_last = min(seq_stop + dbi.motif.distance, n) - candidate_last_offset
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


const PREFIX_HASH_SCAN_V32U8 = Vec{32, UInt8}
const PREFIX_HASH_SCAN_AVX2_FEATURE = UInt32(32 * 2 + 5)
const PREFIX_HASH_SCAN_BMI2_FEATURE = UInt32(32 * 2 + 8)
const PREFIX_HASH_SCAN_PDEP_EVEN = UInt64(0x5555555555555555)
const PREFIX_HASH_SCAN_PDEP_ODD = UInt64(0xaaaaaaaaaaaaaaaa)

@inline function can_use_prefix_hash_scan_simd()
    (Sys.ARCH === :x86_64 || Sys.ARCH === :i686) || return false
    return ccall(:jl_test_cpu_feature, Bool, (UInt32,), PREFIX_HASH_SCAN_AVX2_FEATURE) &&
        ccall(:jl_test_cpu_feature, Bool, (UInt32,), PREFIX_HASH_SCAN_BMI2_FEATURE)
end

@inline function prefix_hash_scan_movemask(v::PREFIX_HASH_SCAN_V32U8)
    Base.llvmcall(
        ("""
         declare i32 @llvm.x86.avx2.pmovmskb(<32 x i8>)
         define i32 @entry(<32 x i8> %0) #0 {
             %res = call i32 @llvm.x86.avx2.pmovmskb(<32 x i8> %0)
             ret i32 %res
         }
         attributes #0 = { "target-features"="+avx2" }
         """, "entry"),
        UInt32, Tuple{PREFIX_HASH_SCAN_V32U8}, v)
end

@inline function prefix_hash_scan_pdep(value::UInt64, mask::UInt64)
    Base.llvmcall(
        ("""
         declare i64 @llvm.x86.bmi.pdep.64(i64, i64)
         define i64 @entry(i64 %0, i64 %1) #0 {
             %res = call i64 @llvm.x86.bmi.pdep.64(i64 %0, i64 %1)
             ret i64 %res
         }
         attributes #0 = { "target-features"="+bmi2" }
         """, "entry"),
        UInt64, Tuple{UInt64, UInt64}, value, mask)
end

@inline function prefix_hash_scan_ascii_mask(
    folded::PREFIX_HASH_SCAN_V32U8, upper::UInt8)

    bytes = vifelse(
        folded == PREFIX_HASH_SCAN_V32U8(upper),
        PREFIX_HASH_SCAN_V32U8(0xff),
        PREFIX_HASH_SCAN_V32U8(0x00),
    )
    return UInt64(prefix_hash_scan_movemask(bytes))
end

@inline function prefix_hash_scan_raw_profile64(
    raw::AbstractVector{UInt8}, start_pos::Int)

    case_mask = PREFIX_HASH_SCAN_V32U8(0xdf)
    chunk0 = vload(PREFIX_HASH_SCAN_V32U8, pointer(raw, start_pos)) & case_mask
    chunk1 = vload(PREFIX_HASH_SCAN_V32U8, pointer(raw, start_pos + 32)) & case_mask
    profile(upper) =
        prefix_hash_scan_ascii_mask(chunk0, upper) |
        (prefix_hash_scan_ascii_mask(chunk1, upper) << 32)
    return (
        profile(UInt8('A')),
        profile(UInt8('C')),
        profile(UInt8('G')),
        profile(UInt8('T')),
    )
end

@inline function prefix_hash_scan_valid23(exact::UInt128)
    valid2 = exact & (exact >> 1)
    valid4 = valid2 & (valid2 >> 2)
    valid8 = valid4 & (valid4 >> 4)
    valid16 = valid8 & (valid8 >> 8)
    return valid16 & (UInt128(valid4) >> 16) &
        (UInt128(valid2) >> 20) & (exact >> 22)
end

@inline function prefix_hash_scan_pack_codes(
    low_bits::UInt64, high_bits::UInt64)

    return UInt32(
        prefix_hash_scan_pdep(low_bits, PREFIX_HASH_SCAN_PDEP_EVEN) |
        prefix_hash_scan_pdep(high_bits, PREFIX_HASH_SCAN_PDEP_ODD))
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

@inline function prefix_hash_scan_raw_hash(
    raw::AbstractVector{UInt8}, candidate_start::Int, is_antisense::Bool)

    hash = UInt32(0)
    positions = if is_antisense
        (candidate_start + 3):(candidate_start + 18)
    else
        (candidate_start + 19):-1:(candidate_start + 4)
    end
    @inbounds for pos in positions
        code = prefix_hash_scan_raw_code(raw[pos])
        code == 0xff && return nothing
        is_antisense && (code = UInt8(3) - code)
        hash = (hash << 2) | UInt32(code)
    end
    return hash
end

function cas9_prefix_scan_bounds_raw(raw::AbstractVector{UInt8}, dbi::DBInfo)
    n = length(raw)
    geometry = CAS9_D3_PREFIX_SCAN_GEOMETRY
    candidate_bases = prefix_scan_candidate_bases(geometry)
    candidate_last_offset = prefix_scan_candidate_last_offset(geometry)
    n < candidate_bases && return nothing
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
    plus_first = max(seq_start - dbi.motif.distance, 1)
    plus_last = seq_stop - candidate_last_offset
    minus_first = seq_start
    minus_last = min(seq_stop + dbi.motif.distance, n) - candidate_last_offset
    firsts = Int[]
    lasts = Int[]
    plus_first <= plus_last && (push!(firsts, plus_first); push!(lasts, plus_last))
    minus_first <= minus_last && (push!(firsts, minus_first); push!(lasts, minus_last))
    isempty(firsts) && return nothing
    return minimum(firsts), maximum(lasts), plus_first, plus_last, minus_first, minus_last
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

function scan_cas9_prefix_hits_raw_range_impl!(
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
    ::Val{Bucketed}) where Bucketed

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
        a0, c0, g0, t0 = prefix_hash_scan_raw_profile64(raw, block_start)
    end
    while block_start + 127 <= n && block_start + 63 <= candidate_last
        a1, c1, g1, t1 = prefix_hash_scan_raw_profile64(raw, block_start + 64)
        a = UInt128(a0) | (UInt128(a1) << 64)
        c = UInt128(c0) | (UInt128(c1) << 64)
        g = UInt128(g0) | (UInt128(g1) << 64)
        t = UInt128(t0) | (UInt128(t1) << 64)
        valid = UInt64(prefix_hash_scan_valid23(a | c | g | t) & UInt128(typemax(UInt64)))
        count = min(64, candidate_last - block_start + 1)
        count_mask = count == 64 ? typemax(UInt64) : (UInt64(1) << count) - 1
        valid &= count_mask
        plus_mask = valid & UInt64((g >> 21) & UInt128(typemax(UInt64))) &
            UInt64((g >> 22) & UInt128(typemax(UInt64)))
        minus_mask = valid & UInt64(c & UInt128(typemax(UInt64))) &
            UInt64((c >> 1) & UInt128(typemax(UInt64)))
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
            hash = prefix_hash_scan_pack_codes(low16, high16)
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
            low16 = UInt64((low >> (bit + 3)) & UInt128(0xffff))
            high16 = UInt64((high >> (bit + 3)) & UInt128(0xffff))
            hash = xor(
                prefix_hash_scan_reverse_codes(
                    prefix_hash_scan_pack_codes(low16, high16)),
                typemax(UInt32),
            )
            prefix_hash_scan_record_candidate!(
                minus_hits, minus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
        block_start += 64
        a0, c0, g0, t0 = a1, c1, g1, t1
    end

    @inbounds for candidate_start in block_start:candidate_last
        valid = true
        for pos in candidate_start:(candidate_start + 22)
            if prefix_hash_scan_raw_code(raw[pos]) == 0xff
                valid = false
                break
            end
        end
        valid || continue
        if plus_first <= candidate_start <= plus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start + 21]) == 2 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 22]) == 2
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash(raw, candidate_start, false)
            prefix_hash_scan_record_candidate!(
                plus_hits, plus_candidates, query, candidate_start, hash,
                Val(Bucketed))
        end
        if minus_first <= candidate_start <= minus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start]) == 1 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 1]) == 1
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash(raw, candidate_start, true)
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

function scan_cas9_prefix_hits_raw_range!(
    plus_hits::Vector{PrefixHashScanHit},
    minus_hits::Vector{PrefixHashScanHit},
    raw::AbstractVector{UInt8},
    query,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int)

    return scan_cas9_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, nothing, nothing, nothing, nothing, nothing,
        raw, query, candidate_first, candidate_last, plus_first, plus_last,
        minus_first, minus_last, Val(false))
end

function scan_cas9_prefix_hits_raw_range_bucketed!(
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
    minus_last::Int)

    return scan_cas9_prefix_hits_raw_range_impl!(
        plus_hits, minus_hits, plus_candidates, minus_candidates,
        plus_radix_scratch, minus_radix_scratch, radix_counts, raw, query,
        candidate_first, candidate_last, plus_first, plus_last, minus_first,
        minus_last, Val(true))
end

function scan_cas9_prefix_hits_raw_range(
    raw::AbstractVector{UInt8},
    query,
    candidate_first::Int,
    candidate_last::Int,
    plus_first::Int,
    plus_last::Int,
    minus_first::Int,
    minus_last::Int)

    plus_hits = PrefixHashScanHit[]
    minus_hits = PrefixHashScanHit[]
    motif_candidates = scan_cas9_prefix_hits_raw_range!(
        plus_hits, minus_hits, raw, query, candidate_first, candidate_last,
        plus_first, plus_last, minus_first, minus_last)
    return plus_hits, minus_hits, motif_candidates
end

function scan_cas9_prefix_hits_raw(
    raw::AbstractVector{UInt8},
    dbi::DBInfo,
    query,
    stats::Union{Nothing, PrefixHashScanStats} = nothing;
    scan_threads::Int = Threads.nthreads())

    bounds = cas9_prefix_scan_bounds_raw(raw, dbi)
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
        Threads.@spawn scan_cas9_prefix_hits_raw_range(
            raw, query, first(range), last(range),
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
    if stats !== nothing
        stats.motif_candidates += motif_candidates
    end
    return plus_hits, minus_hits
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

function scan_verify_cas9_prefix_raw_range!(
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
    stats::S) where {S <: Union{Nothing, PrefixHashScanStats}}

    motif_candidates = 0
    candidate_first > candidate_last && return plus, minus
    n = length(raw)
    block_start = candidate_first

    if block_start + 127 <= n && block_start + 63 <= candidate_last
        a0, c0, g0, t0 = prefix_hash_scan_raw_profile64(raw, block_start)
    end
    while block_start + 127 <= n && block_start + 63 <= candidate_last
        a1, c1, g1, t1 = prefix_hash_scan_raw_profile64(raw, block_start + 64)
        a = UInt128(a0) | (UInt128(a1) << 64)
        c = UInt128(c0) | (UInt128(c1) << 64)
        g = UInt128(g0) | (UInt128(g1) << 64)
        t = UInt128(t0) | (UInt128(t1) << 64)
        valid = UInt64(prefix_hash_scan_valid23(a | c | g | t) & UInt128(typemax(UInt64)))
        count = min(64, candidate_last - block_start + 1)
        count_mask = count == 64 ? typemax(UInt64) : (UInt64(1) << count) - 1
        valid &= count_mask
        plus_mask = valid & UInt64((g >> 21) & UInt128(typemax(UInt64))) &
            UInt64((g >> 22) & UInt128(typemax(UInt64)))
        minus_mask = valid & UInt64(c & UInt128(typemax(UInt64))) &
            UInt64((c >> 1) & UInt128(typemax(UInt64)))
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
            hash = prefix_hash_scan_pack_codes(low16, high16)
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                plus, raw, candidate_start, mask, global_offset, dbi, false,
                guides_, myers_profiles, distance, stats)
        end

        while minus_mask != 0
            bit = trailing_zeros(minus_mask)
            minus_mask &= minus_mask - 1
            candidate_start = block_start + bit
            minus_first <= candidate_start <= minus_last || continue
            motif_candidates += 1
            low16 = UInt64((low >> (bit + 3)) & UInt128(0xffff))
            high16 = UInt64((high >> (bit + 3)) & UInt128(0xffff))
            hash = xor(
                prefix_hash_scan_reverse_codes(
                    prefix_hash_scan_pack_codes(low16, high16)),
                typemax(UInt32),
            )
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                minus, raw, candidate_start, mask, global_offset, dbi, true,
                guides_, myers_profiles, distance, stats)
        end
        block_start += 64
        a0, c0, g0, t0 = a1, c1, g1, t1
    end

    @inbounds for candidate_start in block_start:candidate_last
        valid = true
        for pos in candidate_start:(candidate_start + 22)
            if prefix_hash_scan_raw_code(raw[pos]) == 0xff
                valid = false
                break
            end
        end
        valid || continue
        if plus_first <= candidate_start <= plus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start + 21]) == 2 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 22]) == 2
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash(raw, candidate_start, false)
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                plus, raw, candidate_start, mask, global_offset, dbi, false,
                guides_, myers_profiles, distance, stats)
        end
        if minus_first <= candidate_start <= minus_last &&
                prefix_hash_scan_raw_code(raw[candidate_start]) == 1 &&
                prefix_hash_scan_raw_code(raw[candidate_start + 1]) == 1
            motif_candidates += 1
            hash = prefix_hash_scan_raw_hash(raw, candidate_start, true)
            mask = prefix_hash_scan_candidate_mask(query, hash)
            mask == 0 || evaluate_prefix_hash_scan_candidate!(
                minus, raw, candidate_start, mask, global_offset, dbi, true,
                guides_, myers_profiles, distance, stats)
        end
    end
    if stats !== nothing
        stats.motif_candidates += motif_candidates
    end
    return plus, minus
end
