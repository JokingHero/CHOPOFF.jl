# Verify LLVM emits SIMD instructions for compute_block with Vec{4, UInt64}
#
# Usage: julia --project=. scripts/verify_simd_codegen.jl

using CHOPOFF
using SIMD
using InteractiveUtils

function check_vectorization()
    VecT = Vec{4, UInt64}
    args = Tuple{VecT, VecT, VecT, VecT, VecT}

    io = IOBuffer()
    code_llvm(io, CHOPOFF.Sassy.compute_block, args; debuginfo=:none)
    ir = String(take!(io))

    has_vec4 = occursin("<4 x i64>", ir)

    println("=== SIMD Codegen Verification ===")
    println("Function: compute_block(Vec{4,UInt64} × 5)")
    println("Vectorized (<4 x i64>): ", has_vec4 ? "YES" : "NO")

    if has_vec4
        # Count vector operations
        vec_ops = count(r"<4 x i64>", ir)
        println("Vector operation count: $vec_ops")
    else
        println("WARNING: scalar fallback detected — inner DP loop may be 4x slower than expected")
        println()
        println("--- LLVM IR (first 60 lines) ---")
        for (i, line) in enumerate(split(ir, '\n'))
            i > 60 && break
            println(line)
        end
    end

    println()
    println("--- Full LLVM IR ---")
    println(ir)

    return has_vec4
end

check_vectorization()
