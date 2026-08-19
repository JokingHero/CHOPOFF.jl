#!/usr/bin/env julia

using CHOPOFF
using SIMD
using InteractiveUtils

function native_code(f, types)
    io = IOBuffer()
    code_native(io, f, types; debuginfo = :none, syntax = :intel)
    return lowercase(String(take!(io)))
end

function verify_avx512()
    CHOPOFF.Sassy.can_use_avx512() || error("AVX-512F/BW is unavailable")
    VecT = Vec{8, UInt64}
    core = native_code(CHOPOFF.Sassy.compute_block, Tuple{VecT, VecT, VecT, VecT, VecT})
    encoder = native_code(
        CHOPOFF.Sassy.encode_block_avx512!,
        Tuple{Matrix{UInt64}, Int, Vector{UInt8}, Int, Int, Vector{UInt8}},
    )
    occursin("zmm", core) || error("Vec{8,UInt64} core has no ZMM operations")
    occursin("zmm", encoder) || error("AVX-512 encoder has no ZMM operations")
    println("AVX-512 codegen verified: ZMM Myers core and 64-byte encoder")
end

function verify_avx2_safe()
    VecT = Vec{4, UInt64}
    core = native_code(CHOPOFF.Sassy.compute_block, Tuple{VecT, VecT, VecT, VecT, VecT})
    encoder = native_code(
        CHOPOFF.Sassy.encode_block_avx2!,
        Tuple{Matrix{UInt64}, Int, Vector{UInt8}, Int, Int, Vector{UInt8}},
    )
    minima = native_code(
        CHOPOFF.Sassy.prefix_min_val,
        Tuple{UInt64, UInt64, Val{false}},
    )
    occursin("ymm", core) || error("Vec{4,UInt64} core has no YMM operations")
    occursin("ymm", encoder) || error("AVX2 encoder has no YMM operations")
    occursin("zmm", core * encoder * minima) && error("Safe AVX2 path contains ZMM operations")
    occursin("pext", minima) && error("Safe minima path contains PEXT")
    println("AVX2-safe codegen verified: YMM core/encoder and no PEXT or ZMM")
end

mode = isempty(ARGS) ? "auto" : ARGS[1]
if mode == "auto"
    verify_avx2_safe()
    CHOPOFF.Sassy.can_use_avx512() && verify_avx512()
elseif mode == "avx512"
    verify_avx512()
elseif mode == "avx2_safe"
    verify_avx2_safe()
else
    error("Usage: verify_simd_codegen.jl [auto|avx512|avx2_safe]")
end
