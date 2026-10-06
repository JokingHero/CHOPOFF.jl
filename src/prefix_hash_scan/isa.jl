# Instruction-set layer: CPU feature detection, backend resolution, and the
# primitives whose implementation depends on the backend (`Val(:avx2)`,
# `Val(:avx512)`, `Val(:portable)`). All x86 intrinsics live in this file;
# `:portable` methods are plain Julia and run on any CPU.

const PREFIX_HASH_SCAN_V32U8 = Vec{32, UInt8}
const PREFIX_HASH_SCAN_V64U8 = Vec{64, UInt8}
const PREFIX_HASH_SCAN_AVX2_FEATURE = UInt32(32 * 2 + 5)
const PREFIX_HASH_SCAN_BMI2_FEATURE = UInt32(32 * 2 + 8)
const PREFIX_HASH_SCAN_AVX512F_FEATURE = UInt32(32 * 2 + 16)
const PREFIX_HASH_SCAN_AVX512BW_FEATURE = UInt32(32 * 2 + 30)
const PREFIX_HASH_SCAN_PDEP_EVEN = UInt64(0x5555555555555555)
const PREFIX_HASH_SCAN_PDEP_ODD = UInt64(0xaaaaaaaaaaaaaaaa)
const PREFIX_HASH_SCAN_SIMD_BACKENDS = (:auto, :avx512, :avx2, :portable)
const PREFIX_HASH_SCAN_AVX512_AUTO_TARGETS = (
    ("skylake-avx512", :cas9),
    ("skylake-avx512", :cas12a),
)

@inline function can_use_prefix_hash_scan_avx2()
    (Sys.ARCH === :x86_64 || Sys.ARCH === :i686) || return false
    return ccall(:jl_test_cpu_feature, Bool, (UInt32,), PREFIX_HASH_SCAN_AVX2_FEATURE) &&
        ccall(:jl_test_cpu_feature, Bool, (UInt32,), PREFIX_HASH_SCAN_BMI2_FEATURE)
end

@inline function can_use_prefix_hash_scan_avx512()
    (Sys.ARCH === :x86_64 || Sys.ARCH === :i686) || return false
    return ccall(:jl_test_cpu_feature, Bool, (UInt32,), PREFIX_HASH_SCAN_AVX512F_FEATURE) &&
        ccall(:jl_test_cpu_feature, Bool, (UInt32,), PREFIX_HASH_SCAN_AVX512BW_FEATURE) &&
        ccall(:jl_test_cpu_feature, Bool, (UInt32,), PREFIX_HASH_SCAN_BMI2_FEATURE)
end

function resolve_prefix_hash_scan_simd_backend(
    requested::Symbol = :auto;
    scan_kind::Symbol = :generic,
    cpu_name::AbstractString = Sys.CPU_NAME,
    avx2::Bool = can_use_prefix_hash_scan_avx2(),
    avx512::Bool = can_use_prefix_hash_scan_avx512())

    requested in PREFIX_HASH_SCAN_SIMD_BACKENDS ||
        error("simd_backend must be :auto, :avx2, :avx512, or :portable.")
    if requested == :auto
        if avx512 && (lowercase(cpu_name), scan_kind) in
                PREFIX_HASH_SCAN_AVX512_AUTO_TARGETS
            return :avx512
        end
        return avx2 ? :avx2 : :portable
    elseif requested == :avx512
        avx512 || error("simd_backend=:avx512 requires AVX-512F, AVX-512BW, and BMI2.")
    elseif requested == :avx2
        avx2 || error("simd_backend=:avx2 requires AVX2 and BMI2.")
    end
    return requested
end

@inline default_prefix_hash_scan_simd_backend() =
    Val(resolve_prefix_hash_scan_simd_backend(:auto))

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

@inline function prefix_hash_scan_ascii_mask(
    folded::PREFIX_HASH_SCAN_V64U8, upper::UInt8)

    Base.llvmcall(
        ("""
         define i64 @entry(<64 x i8> %0, <64 x i8> %1) #0 {
             %cmp = icmp eq <64 x i8> %0, %1
             %mask = bitcast <64 x i1> %cmp to i64
             ret i64 %mask
         }
         attributes #0 = { alwaysinline "target-features"="+avx512f,+avx512bw" }
         """, "entry"),
        UInt64,
        Tuple{PREFIX_HASH_SCAN_V64U8, PREFIX_HASH_SCAN_V64U8},
        folded,
        PREFIX_HASH_SCAN_V64U8(upper),
    )
end

@inline function prefix_hash_scan_raw_profile64(
    raw::AbstractVector{UInt8}, start_pos::Int, ::Val{:avx2})

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

@inline function prefix_hash_scan_raw_profile64(
    raw::AbstractVector{UInt8}, start_pos::Int, ::Val{:avx512})

    folded = vload(PREFIX_HASH_SCAN_V64U8, pointer(raw, start_pos)) &
        PREFIX_HASH_SCAN_V64U8(0xdf)
    return (
        prefix_hash_scan_ascii_mask(folded, UInt8('A')),
        prefix_hash_scan_ascii_mask(folded, UInt8('C')),
        prefix_hash_scan_ascii_mask(folded, UInt8('G')),
        prefix_hash_scan_ascii_mask(folded, UInt8('T')),
    )
end

const PREFIX_HASH_SCAN_SWAR_LOW7 = 0x7f7f7f7f7f7f7f7f
const PREFIX_HASH_SCAN_SWAR_CASE_MASK = 0xdfdfdfdfdfdfdfdf

# Bit j set when byte j of `word` equals the byte broadcast in `pattern`.
@inline function prefix_hash_scan_swar_eq_mask(word::UInt64, pattern::UInt64)
    x = xor(word, pattern)
    # High bit of each byte set exactly when that byte of x is zero; the
    # masked add cannot carry across bytes.
    zero_bytes = ~(((x & PREFIX_HASH_SCAN_SWAR_LOW7) + PREFIX_HASH_SCAN_SWAR_LOW7) |
        x | PREFIX_HASH_SCAN_SWAR_LOW7)
    # Gather bit 8j to bit 56 + j; the partial products occupy distinct bits.
    return ((zero_bytes >> 7) * 0x0102040810204080) >> 56
end

# Same masks as the SIMD profiles, from eight little-endian 64-bit loads.
@inline function prefix_hash_scan_raw_profile64(
    raw::AbstractVector{UInt8}, start_pos::Int, ::Val{:portable})

    a = c = g = t = UInt64(0)
    GC.@preserve raw begin
        ptr = Ptr{UInt64}(pointer(raw, start_pos))
        for k in 0:7
            word = ltoh(unsafe_load(ptr, k + 1)) & PREFIX_HASH_SCAN_SWAR_CASE_MASK
            lane = 8 * k
            a |= prefix_hash_scan_swar_eq_mask(word, 0x4141414141414141) << lane
            c |= prefix_hash_scan_swar_eq_mask(word, 0x4343434343434343) << lane
            g |= prefix_hash_scan_swar_eq_mask(word, 0x4747474747474747) << lane
            t |= prefix_hash_scan_swar_eq_mask(word, 0x5454545454545454) << lane
        end
    end
    return a, c, g, t
end

@inline function prefix_hash_scan_pack_codes(
    low_bits::UInt64, high_bits::UInt64)

    return UInt32(
        prefix_hash_scan_pdep(low_bits, PREFIX_HASH_SCAN_PDEP_EVEN) |
        prefix_hash_scan_pdep(high_bits, PREFIX_HASH_SCAN_PDEP_ODD))
end

@inline prefix_hash_scan_pack_codes(low_bits::UInt64, high_bits::UInt64, ::Val) =
    prefix_hash_scan_pack_codes(low_bits, high_bits)

# Moves bit i of the low 16 bits to bit 2i.
@inline function prefix_hash_scan_spread_even(x::UInt64)
    x &= UInt64(0xffff)
    x = (x | (x << 8)) & UInt64(0x00ff00ff)
    x = (x | (x << 4)) & UInt64(0x0f0f0f0f)
    x = (x | (x << 2)) & UInt64(0x33333333)
    return (x | (x << 1)) & UInt64(0x55555555)
end

@inline prefix_hash_scan_pack_codes(
    low_bits::UInt64, high_bits::UInt64, ::Val{:portable}) =
    UInt32(prefix_hash_scan_spread_even(low_bits) |
        (prefix_hash_scan_spread_even(high_bits) << 1))
