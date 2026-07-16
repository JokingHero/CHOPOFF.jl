# Geometry-neutral scalar, AVX2/BMI2, and compact-lookup helpers.

@inline function prefix_hash_scan_twobit_nibble(nibble::UInt8)
    nibble == 0x01 && return UInt8(0)
    nibble == 0x02 && return UInt8(1)
    nibble == 0x04 && return UInt8(2)
    nibble == 0x08 && return UInt8(3)
    return UInt8(0xff)
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
             %res = call i64 @llvm.x86.bmi.pdep.64(i64 %0, i64 %1) #0
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
