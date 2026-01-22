### Algorithm Reference: Sassy's Transposed SIMD-Myers

**Core Concept:**
Instead of computing the edit distance column-by-column (scanning text), we compute it row-by-row (scanning pattern).

1.  **State Vectors:** `Hp` (Horizontal Positive) and `Hm` (Horizontal Negative). These are `UInt64` vectors where each bit corresponds to a position in the **Text**.
2.  **SIMD Tiling:** We process `L` lanes simultaneously. If `L=4` (AVX2 equivalent), we process 4 chunks of text ($4 \times 64 = 256$ characters) in parallel per instruction.
3.  **Data Flow:** We iterate through the **Pattern** characters ($j=1..m$). In each step, we update the `Hp/Hm` vectors for the entire text block and carry the Vertical Deltas ($Vp, Vm$) to the next pattern row.

#### 1. The SIMD Backend (Pseudocode/Julia-like)

We need a structure to hold the vertical deltas (`V`) which flow _down_ the rows, and the horizontal deltas (`H`) which stay within the text block.

```julia
# Concept: A SIMD vector of 4 UInt64s
# In Julia: Use tuples or SIMD.jl (Vec{4, UInt64})
struct SimdState
    # Vertical deltas (passed from pattern row j to j+1)
    Vp::Vec4_UInt64
    Vm::Vec4_UInt64

    # Horizontal deltas (state of the current text block)
    Hp::Vec4_UInt64
    Hm::Vec4_UInt64
end
```

#### 2. Text Profiling (The "Transposed" Pre-computation)

Because we iterate through the _pattern_, we cannot look up the character match mask for the pattern (like standard Myers). Instead, we must convert the **Text** chunks into profile blocks so we can quickly ask: "Which positions in this 64-char text block match 'A'?"

```python
# Reference: src/profiles/iupac.rs -> encode_ref
# Input: A 64-byte slice of the Genome
# Output: An array of bitmasks, one for each IUPAC code.
function encode_text_block(text_slice_64_bytes):
    # This effectively "splats" the text.
    # If text is "ACGT...",
    # masks['A'] will have bit 0 set.
    # masks['C'] will have bit 1 set.

    masks = Dict() # Maps char_code -> UInt64

    # Sassy optimizes this using SIMD shuffling (pshufb),
    # but logically it is:
    for i in 0..63:
        char = text_slice_64_bytes[i]
        # Set bit i for every base that matches 'char'
        # e.g. if char is 'N', it sets bits in masks['A'], masks['C']...
        update_masks(masks, char, bit_index=i)

    return masks
```

#### 3. The Transposed Block Step (The Math)

This is the logic from `src/bitpacking.rs` (`compute_block_simd`). Note the shift direction: `Hp << 1`. In horizontal bit-parallelism, shifting moves us to the next text character.

```python
# Inputs:
#   Hp_in, Hm_in: Horizontal deltas from previous pattern row
#   Vp, Vm: Vertical deltas accumulated so far
#   Eq: Bitmask of matches for the current Pattern character against the 64 Text chars
# Returns:
#   (Hp_out, Hm_out, Vp_out, Vm_out)
function step_transposed_myers(Hp_in, Hm_in, Vp, Vm, Eq):
    # 1. Update Vertical Deltas (Match logic)
    Xv = Eq | Vm

    # 2. Horizontal Update Logic (The "Folding" magic)
    # Note: Hp_in/Hm_in act as the "carry" in this orientation

    # Standard Myers Eq adjustment for horizontal input
    # In Sassy: Eq = Eq | Hm_in (conceptually handling gaps)

    # Core Logic
    t1 = (Eq & Vp) + Vp   # Arithmetic add handles the carry chain
    Xh = (t1 ^ Vp) | Eq

    Ph = Vm | ~(Xh | Vp)
    Mh = Vp & Xh

    # 3. Shift Horizontal Deltas (Move to next text pos)
    # We shift in a 1 (gap penalty) from the left
    Ph_shifted = (Ph << 1) | 1
    Mh_shifted = (Mh << 1)

    # 4. Update Horizontal State (for next pattern char)
    Hp_out = Mh_shifted | ~(Xv | Ph_shifted)
    Hm_out = Ph_shifted & Xv

    # 5. Update Vertical State (to pass to next pattern char)
    # In transposed mode, we keep track of Vp/Vm to affect the NEXT row
    # Ideally, Vp/Vm are actually implicit or used to compute the next H

    # Sassy implementation specifically updates H in place and returns new V
    # See compute_block_simd logic:
    Hp_next = Hp_in # (Updated via logic above, effectively Ph_shifted?)
    # ...

    # LET'S USE THE EXACT SASSY LOGIC from `compute_block_simd`:
    Xv = Eq | Vm
    Eq_local = Eq | Hm_in

    Xh = (((Eq_local & Vp) + Vp) ^ Vp) | Eq_local

    Ph = Vm | ~(Xh | Vp)
    Mh = Vp & Xh

    # Shift horizontal deltas (moving right in text)
    # Hp_in serves as the "carry in" bit if we were chaining 64-bit blocks
    # But usually we assume independent overlapped blocks.
    Ph_shifted = (Ph << 1) | 1
    Mh_shifted = (Mh << 1)

    Hp_out = Ph_shifted
    Hm_out = Mh_shifted

    # Calculate new Vertical deltas
    # (These flow down to the next pattern row)
    Pv_out = Mh | ~(Xv | Ph)
    Mv_out = Ph & Xv

    return (Hp_out, Hm_out, Pv_out, Mv_out)
```

#### 4. The Main Search Function (Tiling)

This replicates `src/search.rs`.

```julia
function search_genome_simd(guide, genome, k)
    # 1. Chunking
    # Divide genome into chunks of size (LANES * 64).
    # e.g., 256 characters per chunk.
    # Ensure chunks overlap by (Length(guide) + k) to catch boundary matches.

    m = length(guide)

    # 2. Encode Guide
    # Convert guide string to array of indices (0=A, 1=C...)
    guide_idx = encode_pattern(guide)

    # 3. Main Loop over Text Chunks
    for chunk in genome_chunks
        # A. Pre-compute Text Masks (The "LaneState")
        # Load 4 blocks of 64 chars.
        # Create masks[Char] -> Vec4_UInt64
        # e.g., text_masks[0] contains bits for 'A' across all 4 lanes
        lane_profiles = profile_text_chunk_simd(chunk)

        # B. Initialize State
        # Hp, Hm initialized (usually Hp=1s, Hm=0s for global alignment start)
        Hp = Vec4(All_Ones)
        Hm = Vec4(Zeros)

        # Vp, Vm initialized (0s)
        Vp = Vec4(Zeros)
        Vm = Vec4(Zeros)

        Score = Vec4(Zeros) # Or keeping track of current edit distance

        # C. Iterate Pattern (Rows)
        for j in 1..m
            char_code = guide_idx[j]

            # 1. Get Equality Mask for this pattern char across all 4 text lanes
            Eq = lane_profiles[char_code]

            # 2. Run Transposed Myers Step
            (Hp, Hm, Vp, Vm) = step_transposed_myers(Hp, Hm, Vp, Vm, Eq)

            # 3. Update Score (Vertical movement)
            # In standard Myers, score is tracked at the last column.
            # In transposed, we track score at the bottom of the column (current row).
            # Score += Vp - Vm
            Score = Score + Vp - Vm
        end

        # D. Check Results
        # 'Score' now contains the edit distance for the alignment ending at
        # the current text position for every lane.
        # Wait, Score is a Vec4. Each element represents the score at the specific
        # text offset inside that 64-bit word?
        # NO. In transposed Myers, the score is distributed.

        # CORRECTION ON SCORING:
        # In transposed Myers, we don't get a simple "Score" register.
        # We have the final H vector (Hp, Hm) after processing the whole pattern.
        # We calculate the score for every bit position `b` (0..63) by summing H deltas?
        # NO. Sassy tracks score diagonally or re-calculates it.

        # Let's look at `search.rs` -> `find_minima_with_overhang`.
        # It calculates the exact score for every bit position `b` by walking `V`.
        # Score[b] = Initial_Score + count_ones(Vp & mask) - count_ones(Vm & mask).

        # Simplified Sassy Scoring Strategy:
        # 1. We track a base `Dist` (scalar) for the lane.
        # 2. After the loop, we iterate bits `b` from 0..63.
        # 3. Current_Cost = Dist + (Vp >> b & 1) - (Vm >> b & 1).
        # 4. If Current_Cost <= k, we have a match ending at `ChunkStart + b`.
    end
end
```

### Prompt for LLM

Use this prompt to get the exact architecture you need:

> "I need to implement a high-performance CRISPR search in Julia by porting the specific architecture of the Rust `sassy` library.
>
> **The Architecture:**
> Unlike standard Myers (which packs the pattern into a bitvector), Sassy uses **Horizontal Bit-Parallelism (Transposed Myers)** combined with **SIMD Tiling**.
>
> **Requirements:**
>
> 1.  **Data Layout:** Use `UInt64` to represent a block of 64 **Genome** characters (not pattern characters).
> 2.  **SIMD:** Process 4 of these 64-char blocks in parallel using `SIMD.jl` (i.e., `Vec{4, UInt64}`). This means processing 256 genome characters per iteration.
> 3.  **Pattern Loop:** Iterate through the **Guide/Pattern** characters one by one.
>     - For each pattern char, load the pre-computed match mask for the 4 text blocks.
>     - Update the Horizontal Deltas (`Hp`, `Hm`) and Vertical Deltas (`Vp`, `Vm`) using the bitwise Myers formula.
> 4.  **Text Profiling:** Implement a function that takes a 256-byte text chunk and converts it into IUPAC bitmasks (A, C, G, T, N, etc.) usable by the search loop.
> 5.  **Scoring:** Since the bitvectors represent text positions, the edit distance is distributed. Implement the logic to extract the specific edit distance for each of the 64 positions in the block after the pattern loop finishes (using `popcount` logic on the V-vectors).
> 6.  **PAM Filter:** Once a bit position `i` is identified as a match (score $\le k$), verify the PAM (NGG) manually.
>
> Please generate:
>
> 1.  The `step_transposed_myers` function (SIMD-compatible).
> 2.  The text-chunk profiling logic.
> 3.  The main search loop structure handling overlaps and scoring."
