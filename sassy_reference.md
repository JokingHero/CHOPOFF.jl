To replicate the specific Horizontal Bit-Parallelism (Transposed Myers) + SIMD Tiling architecture used in Sassy for CRISPR search, you can inspect these files in "sassy" folder with the original sassy implementation in rust.
These files contain the exact mathematical logic for the "Transposed" step, the memory layout for the text chunks, and the specific scoring logic used to extract edit distances from the bit-vectors.

1. src/bitpacking.rs
Why: This contains the mathematical core.
Function: compute_block_simd
Significance: This implements the "Transposed Myers" step. Unlike standard Myers, it updates vertical deltas (v) and horizontal deltas (h) for a SIMD vector of text blocks. This is the equation the Julia code must replicate exactly.

2. src/search.rs
Why: This contains the search engine and tiling logic.
Structs: Searcher, LaneState
Functions: search_internal, find_minima_with_overhang
Significance:
search_internal: Shows the main loop structure. It iterates through the text in chunks (SIMD Lanes), and then iterates through the Pattern characters (Horizontal parallelism).
LaneState: Shows how text is buffered into 64-byte blocks.
find_minima_with_overhang: Shows the complex logic required to convert the final bit-vectors (Vp, Vm) back into an actual integer edit distance (score) for every position in the text.

3. src/profiles/iupac.rs
Why: This contains the SIMD-optimized text profiling.
Function: encode_ref
Significance: Because Sassy iterates over the Pattern, it must pre-process the Text chunks. This file shows how a 64-byte text block is converted into bitmasks (e.g., "which bits in this u64 correspond to 'A'?"). It handles the IUPAC logic (N, R, Y, etc.) efficiently.

4. bin/crispr.rs
Why: This contains the CRISPR-specific high-level logic.
Significance: It shows how to invoke the searcher with a filter_fn. This is crucial for CRISPR because you need to check the PAM sequence (e.g., NGG) strictly after finding a candidate match. It shows how to combine the approximate search (edit distance) with the strict suffix check.