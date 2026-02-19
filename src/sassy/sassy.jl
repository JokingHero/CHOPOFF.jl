module Sassy

# --- External Dependencies ---
using BioSequences
using ThreadsX
using FASTX
using TwoBit

# --- Parent Dependencies ---
import ..CHOPOFF: 
    Motif, 
    DBInfo, 
    Offtarget, 
    Loc, 
    align, 
    decode, 
    length_noPAM, 
    locate_telomeres, 
    getchromseq,
    insert_offtarget!

# --- Sub-components ---
# Order matters here due to internal dependencies.

# 1. Constants: IUPAC tables and block sizes
include("constants.jl")

# 2. Encoding: Logic to convert DNA/Text to bitmasks
include("encoding.jl")

# 3. Core: The pure Myers logic and the scalar-emulated SIMD loop
#    (This is the file you will fuzz-test against Rust)
include("core.jl")

# 4. Interface: The high-level functions that connect Sassy logic 
#    to CHOPOFF's structs (Motif, DBInfo) and file I/O.
include("interface.jl")

# --- Exports ---

# The main entry point for the user
export search_sassy

# Useful if you want to run single-guide searches directly
export search_sassy_guide

# Exported specifically for "Friction Test" (Unit Testing)
# You need these visible to run the JSON comparison tests.
export compute_block
export search_sassy_impl
export encode_text_profile
export encode_pattern_sassy

end 