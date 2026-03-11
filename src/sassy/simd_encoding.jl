using SIMD

# AVX2-accelerated text block encoding
# Mirrors Rust encode_ref (profiles/iupac.rs:70-130)
#
# Uses vpshufb for parallel IUPAC table lookup on 32 bytes at a time,
# producing per-base UInt64 bitmasks ~25x faster than scalar encoding.

const V32U8 = Vec{32, UInt8}

# --- CPU Feature Detection ---

const JL_X86_AVX2_FEATURE_ID = UInt32(32 * 2 + 5)  # from Julia features_h.jl
const AVX2_CACHE = Ref{Union{Nothing, Bool}}(nothing)

@inline function can_use_avx2()
    cached = AVX2_CACHE[]
    if cached !== nothing
        return cached::Bool
    end

    enabled = false
    if Sys.ARCH === :x86_64 || Sys.ARCH === :i686
        enabled = ccall(:jl_test_cpu_feature, Bool, (UInt32,), JL_X86_AVX2_FEATURE_ID)
    end

    AVX2_CACHE[] = enabled
    return enabled
end

# --- AVX2 Intrinsics ---

@inline function vpshufb(table::V32U8, idx::V32U8)
    Base.llvmcall(
        ("""
         declare <32 x i8> @llvm.x86.avx2.pshuf.b(<32 x i8>, <32 x i8>)
         define <32 x i8> @entry(<32 x i8> %0, <32 x i8> %1) #0 {
             %res = call <32 x i8> @llvm.x86.avx2.pshuf.b(<32 x i8> %0, <32 x i8> %1)
             ret <32 x i8> %res
         }
         attributes #0 = { "target-features"="+avx2" }
         """, "entry"),
        V32U8, Tuple{V32U8, V32U8}, table, idx)
end

@inline function vpmovmskb(v::V32U8)
    Base.llvmcall(
        ("""
         declare i32 @llvm.x86.avx2.pmovmskb(<32 x i8>)
         define i32 @entry(<32 x i8> %0) #0 {
             %res = call i32 @llvm.x86.avx2.pmovmskb(<32 x i8> %0)
             ret i32 %res
         }
         attributes #0 = { "target-features"="+avx2" }
         """, "entry"),
        UInt32, Tuple{V32U8}, v)
end

# --- PACKED_NIBBLES Table ---
# Encodes the IUPAC table in nibble-packed format for vpshufb lookup.
# Low nibble = IUPAC mask for index 0-15, high nibble = mask for index 16-31.
# Duplicated to fill 32-byte register (vpshufb operates on 16-byte halves).

const PACKED_NIBBLES_VEC = let
    packed = zeros(UInt8, 16)
    for i in 0:15
        lo = IUPAC_SASSY[i + 1] & 0x0F
        hi = (i + 16 < length(IUPAC_SASSY)) ? (IUPAC_SASSY[i + 16 + 1] & 0x0F) : UInt8(0)
        packed[i + 1] = (hi << 4) | lo
    end
    V32U8(ntuple(i -> VecElement(packed[((i - 1) % 16) + 1]), Val(32)))
end

# --- Precomputed Base Mask Splats ---

const BASE_MASKS_ACGT = let
    masks = UInt8[0x01, 0x02, 0x04, 0x08]  # A, C, T, G
    ntuple(i -> V32U8(masks[i]), Val(4))
end

# --- Core AVX2 Encoding ---

# Preprocess one 32-byte text chunk into IUPAC nibble masks via vpshufb.
@inline function _avx2_nibbles(chunk::V32U8)
    MASK4 = V32U8(0x0F)
    MASK5 = V32U8(0x1F)
    FIFTEEN = V32U8(15)

    low4 = chunk & MASK4
    idx5 = chunk & MASK5
    is_hi = idx5 > FIFTEEN

    shuffled = vpshufb(PACKED_NIBBLES_VEC, low4)
    lo_nib = shuffled & MASK4
    hi_nib = (shuffled >>> 4) & MASK4

    return vifelse(is_hi, hi_nib, lo_nib)
end

# Extract a UInt64 bitmask for one base from preprocessed nibble halves.
@inline function _avx2_extract_base(nib0::V32U8, nib1::V32U8, base_splat::V32U8)
    ZERO = V32U8(0x00)
    FF = V32U8(0xFF)

    anded0 = nib0 & base_splat
    anded1 = nib1 & base_splat
    ne0 = vifelse(~(anded0 == ZERO), FF, ZERO)
    ne1 = vifelse(~(anded1 == ZERO), FF, ZERO)

    low = UInt64(vpmovmskb(ne0))
    high = UInt64(vpmovmskb(ne1))
    return (high << 32) | low
end

"""
    encode_block_avx2!(out, lane, text, start_pos, n_bases, bases)

Encode a full 64-byte text block into per-base UInt64 bitmasks using AVX2.
Writes results into `out[lane, 1:n_bases]`.

For the standard 4-base case (ACGT), uses precomputed splat vectors.
For extra IUPAC bases in the pattern (N, R, Y, ...), dynamically creates splats.
"""
@inline function encode_block_avx2!(
    out::Matrix{UInt64}, lane::Int,
    text::AbstractVector{UInt8}, start_pos::Int,
    n_bases::Int, bases::Vector{UInt8}
)
    # Load 64 bytes as 2x32 and preprocess
    chunk0 = vload(V32U8, pointer(text, start_pos))
    chunk1 = vload(V32U8, pointer(text, start_pos + 32))
    nib0 = _avx2_nibbles(chunk0)
    nib1 = _avx2_nibbles(chunk1)

    # Standard ACGT bases (unrolled with precomputed splats)
    @inbounds out[lane, 1] = _avx2_extract_base(nib0, nib1, BASE_MASKS_ACGT[1])
    @inbounds out[lane, 2] = _avx2_extract_base(nib0, nib1, BASE_MASKS_ACGT[2])
    @inbounds out[lane, 3] = _avx2_extract_base(nib0, nib1, BASE_MASKS_ACGT[3])
    @inbounds out[lane, 4] = _avx2_extract_base(nib0, nib1, BASE_MASKS_ACGT[4])

    # Extra IUPAC bases (if pattern contains N, R, Y, etc.)
    for j in 5:n_bases
        @inbounds bm_splat = V32U8(get_iupac_mask(bases[j]))
        @inbounds out[lane, j] = _avx2_extract_base(nib0, nib1, bm_splat)
    end
end
