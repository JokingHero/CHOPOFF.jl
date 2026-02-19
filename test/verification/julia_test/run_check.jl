using JSON
using Test

# Point this to where you saved your split files
include("../../../src/sassy/core.jl") 

# Helper to format numbers as hex for debugging
hex(x) = "0x" * string(x, base=16, pad=16)

println("Loading Test Vectors...")
const vectors = JSON.parsefile("../test_vectors.json")

println("Running Friction Test on compute_block...")

@testset "Sassy Logic Verification" begin
    for (i, case) in enumerate(vectors)
        # 1. Parse Inputs
        # JSON numbers are Int64 by default in Julia, we strictly cast to UInt64
        inp = case["input"]
        hp_in = reinterpret(UInt64, Int64(inp["hp_in"]))
        hm_in = reinterpret(UInt64, Int64(inp["hm_in"]))
        vp    = reinterpret(UInt64, Int64(inp["vp"]))
        vm    = reinterpret(UInt64, Int64(inp["vm"]))
        eq    = reinterpret(UInt64, Int64(inp["eq"]))

        # 2. Parse Expected Outputs
        out = case["output"]
        exp_hp_out = reinterpret(UInt64, Int64(out["hp_out"]))
        exp_hm_out = reinterpret(UInt64, Int64(out["hm_out"]))
        exp_vp     = reinterpret(UInt64, Int64(out["vp"]))
        exp_vm     = reinterpret(UInt64, Int64(out["vm"]))

        # 3. Run Julia Function
        (act_hp_out, act_hm_out, act_vp, act_vm) = compute_block(
            hp_in, hm_in, vp, vm, eq
        )

        # 4. Compare
        # We perform individual checks to give granular error messages
        if act_hp_out != exp_hp_out
            println("\n[FAIL] Case #$i Mismatch on HP_OUT")
            println("Input: hp_in=$(hex(hp_in)) eq=$(hex(eq)) vp=$(hex(vp))")
            println("Expected: $(hex(exp_hp_out))")
            println("Actual:   $(hex(act_hp_out))")
        end
        @test act_hp_out == exp_hp_out

        if act_hm_out != exp_hm_out
            println("\n[FAIL] Case #$i Mismatch on HM_OUT")
            println("Expected: $(hex(exp_hm_out))")
            println("Actual:   $(hex(act_hm_out))")
        end
        @test act_hm_out == exp_hm_out

        if act_vp != exp_vp
            println("\n[FAIL] Case #$i Mismatch on VP")
            println("Expected: $(hex(exp_vp))")
            println("Actual:   $(hex(act_vp))")
        end
        @test act_vp == exp_vp

        if act_vm != exp_vm
            println("\n[FAIL] Case #$i Mismatch on VM")
            println("Expected: $(hex(exp_vm))")
            println("Actual:   $(hex(act_vm))")
        end
        @test act_vm == exp_vm
    end
end

println("Done.")
