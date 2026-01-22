# Sassy SIMD-Optimized Implementation
# Port of Rust Sassy's transposed SIMD-Myers algorithm
#
# Key architecture:
# - Horizontal bit-parallelism: bits represent TEXT positions (not pattern)
# - 4-lane SIMD: processes 256 genome characters per step (4 × 64)
# - Iterate through PATTERN characters, lookup TEXT masks
#
# SIMD is used for fast candidate filtering. Final alignment uses CHOPOFF's align()
# to ensure consistency with linearDB output format.

using SIMD: Vec, vload, vstore
using BioSequences: reverse_complement

# --- Constants ---
const LANES = 4
const BLOCK_SIZE = 64
const CHUNK_SIZE = LANES * BLOCK_SIZE  # 256 characters per SIMD step

# IUPAC encoding: A=1, C=2, T=4, G=8 (matches Sassy Rust)
const IUPAC_SASSY = let
    enc = zeros(UInt8, 32)
    A = UInt8(1 << 0)
    C = UInt8(1 << 1)
    T = UInt8(1 << 2)
    G = UInt8(1 << 3)

    # Map using lower 5 bits of ASCII (& 0x1F)
    enc[Int('A') & 0x1F + 1] = A
    enc[Int('C') & 0x1F + 1] = C
    enc[Int('T') & 0x1F + 1] = T
    enc[Int('U') & 0x1F + 1] = T  # Uracil = Thymine
    enc[Int('G') & 0x1F + 1] = G
    enc[Int('N') & 0x1F + 1] = A | C | T | G

    # IUPAC ambiguity codes
    enc[Int('R') & 0x1F + 1] = A | G
    enc[Int('Y') & 0x1F + 1] = C | T
    enc[Int('S') & 0x1F + 1] = G | C
    enc[Int('W') & 0x1F + 1] = A | T
    enc[Int('K') & 0x1F + 1] = G | T
    enc[Int('M') & 0x1F + 1] = A | C
    enc[Int('B') & 0x1F + 1] = C | G | T
    enc[Int('D') & 0x1F + 1] = A | G | T
    enc[Int('H') & 0x1F + 1] = A | C | T
    enc[Int('V') & 0x1F + 1] = A | C | G
    enc[Int('X') & 0x1F + 1] = 0  # Gap/unknown

    enc
end

@inline function get_iupac_mask(c::UInt8)::UInt8
    @inbounds IUPAC_SASSY[(c & 0x1F) + 1]
end

# --- Text Profile ---
# Pre-computed match masks for a 64-character text block
struct TextBlockProfile
    masks::Vector{UInt64}  # masks[i] has bit j set if text[j] matches bases[i]
end

"""
    encode_text_block(text_slice::AbstractVector{UInt8}, bases::Vector{UInt8})

Encode a 64-byte text slice into IUPAC match masks.
"""
function encode_text_block(text_slice::AbstractVector{UInt8}, bases::Vector{UInt8})
    n_bases = length(bases)
    masks = zeros(UInt64, max(n_bases, 16))

    for (i, base) in enumerate(bases)
        base_mask = get_iupac_mask(base)
        result = UInt64(0)

        for j in 1:min(64, length(text_slice))
            text_char = text_slice[j]
            text_mask = get_iupac_mask(text_char)
            if (base_mask & text_mask) != 0
                result |= (UInt64(1) << (j - 1))
            end
        end

        masks[i] = result
    end

    TextBlockProfile(masks)
end

# --- SIMD Types ---
const SIMDVec = Vec{4, UInt64}
const SIMD_ONES = SIMDVec((~UInt64(0), ~UInt64(0), ~UInt64(0), ~UInt64(0)))
const SIMD_ZEROS = SIMDVec((UInt64(0), UInt64(0), UInt64(0), UInt64(0)))
const SIMD_ONE = SIMDVec((UInt64(1), UInt64(1), UInt64(1), UInt64(1)))

# --- Core SIMD Myers Step ---
@inline function compute_block_simd_row!(
    hp::MVector{4, UInt64},
    hm::MVector{4, UInt64},
    vp::MVector{4, UInt64},
    vm::MVector{4, UInt64},
    eq::NTuple{4, UInt64}
)
    for lane in 1:4
        eq_lane = eq[lane]
        vp_lane = vp[lane]
        vm_lane = vm[lane]
        hp_lane = hp[lane]
        hm_lane = hm[lane]

        # Myers bit-parallel algorithm
        vx = eq_lane | vm_lane
        eq_local = eq_lane | hm_lane
        t1 = (eq_local & vp_lane) + vp_lane
        hx = (t1 ⊻ vp_lane) | eq_local

        hp_new = vm_lane | (~(hx | vp_lane))
        hm_new = vp_lane & hx

        hpw = hp_new >> 63
        hmw = hm_new >> 63

        hp_shifted = (hp_new << 1) | hp_lane
        hm_shifted = (hm_new << 1) | hm_lane

        hp[lane] = hpw
        hm[lane] = hmw

        vp[lane] = hm_shifted | (~(vx | hp_shifted))
        vm[lane] = hp_shifted & vx
    end

    nothing
end

# --- Encode Pattern ---
function encode_pattern_sassy(pattern::AbstractVector{UInt8})
    bases = UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')]
    indices = Int[]

    for c in pattern
        c_upper = c & ~UInt8(0x20)
        idx = findfirst(==(c_upper), bases)
        if idx === nothing
            push!(bases, c_upper)
            push!(indices, length(bases))
        else
            push!(indices, idx)
        end
    end

    (bases, indices)
end

# --- Main SIMD Search (Candidate Filtering) ---
"""
    search_sassy_simd_chunk(pattern_indices, text_chunk, k, bases)

Fast candidate filtering using transposed SIMD Myers.
Returns vector of (global_position, score) for potential matches with score <= k.
These are candidates that need exact verification with align().
"""
function search_sassy_simd_chunk(
    pattern_indices::Vector{Int},
    text_chunk::AbstractVector{UInt8},
    k::Int,
    bases::Vector{UInt8}
)
    m = length(pattern_indices)
    chunk_len = length(text_chunk)

    # Pad text to 256 chars
    padded = fill(UInt8('X'), CHUNK_SIZE)
    copy_len = min(chunk_len, CHUNK_SIZE)
    padded[1:copy_len] .= text_chunk[1:copy_len]

    # Encode 4 text blocks
    lane_profiles = [
        encode_text_block(view(padded, (i-1)*BLOCK_SIZE+1 : i*BLOCK_SIZE), bases)
        for i in 1:LANES
    ]

    # Transposed Myers state
    # hp stores horizontal propagation from previous rows (carry)
    # hm stores horizontal propagation of deletions
    # Initialize to all zeros - the algorithm will set appropriate bits during processing
    hp = [MVector{4, UInt64}(0, 0, 0, 0) for _ in 1:m]
    hm = [MVector{4, UInt64}(0, 0, 0, 0) for _ in 1:m]
    vp = MVector{4, UInt64}(0, 0, 0, 0)
    vm = MVector{4, UInt64}(0, 0, 0, 0)
    dist_to_start = MVector{4, Int}(0, 0, 0, 0)

    # Iterate through pattern rows
    for j in 1:m
        char_idx = pattern_indices[j]

        for lane in 1:LANES
            dist_to_start[lane] += count_ones(hp[j][lane])
            dist_to_start[lane] -= count_ones(hm[j][lane])
        end

        eq = (
            lane_profiles[1].masks[char_idx],
            lane_profiles[2].masks[char_idx],
            lane_profiles[3].masks[char_idx],
            lane_profiles[4].masks[char_idx]
        )

        compute_block_simd_row!(hp[j], hm[j], vp, vm, eq)
    end

    # Extract candidates from V vectors
    matches = Tuple{Int, Int}[]

    for lane in 1:LANES
        lane_vp = vp[lane]
        lane_vm = vm[lane]
        cost = dist_to_start[lane]

        for bit in 0:63
            if (lane_vp & (UInt64(1) << bit)) != 0
                cost += 1
            end
            if (lane_vm & (UInt64(1) << bit)) != 0
                cost -= 1
            end

            global_pos = (lane - 1) * BLOCK_SIZE + bit + 1

            if global_pos <= chunk_len && cost <= k
                push!(matches, (global_pos, cost))
            end
        end
    end

    return matches
end

# --- Result struct for internal use ---
struct SassyMatch
    pos::Int
    dist::Int
    aln_guide::String
    aln_ref::String
end

# --- PAM Verification ---
"""
    check_pam_match(ref_pam, guide_pam)

Check if reference PAM matches guide PAM using IUPAC compatibility.
"""
function check_pam_match(ref_pam::AbstractString, guide_pam::AbstractString)
    length(ref_pam) == length(guide_pam) || return false

    for (r, g) in zip(ref_pam, guide_pam)
        m1 = get_iupac_mask(UInt8(r))
        m2 = get_iupac_mask(UInt8(g))
        (m1 & m2) != 0 || return false
    end

    return true
end

# --- Main Search with align() for exact results ---

"""
    search_sassy(guides, genome_path, motif, output_path; distance=4, early_stopping=nothing)

Search for off-targets using Sassy algorithm.
- SIMD for fast candidate filtering (main sassy optimization)
- CHOPOFF's align() for exact alignment (matches linearDB output)
- Optional early_stopping to limit results per distance tier

Output format matches linearDB for cross-validation testing.
"""
function search_sassy(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif,
    output_path::String;
    distance::Int = 4,
    early_stopping::Union{Vector{Int}, Nothing} = nothing
)
    dbi = DBInfo(genome_path, "sassy_search", motif)

    # Validate guide lengths
    if any(length_noPAM(motif) .!= length.(guides))
        error("Guide queries are not of the correct length to use with this Motif: " * string(motif))
    end

    # Validate early_stopping
    if early_stopping !== nothing
        if length(early_stopping) != (distance + 1)
            error("early_stopping must have $(distance + 1) elements (one for each distance 0 to $distance)")
        end
    end

    # Default early stopping: no limit if not specified
    if early_stopping === nothing
        early_stopping = fill(typemax(Int), distance + 1)
    end

    # Prepare output file
    mkpath(dirname(output_path))
    outfile = open(output_path, "w")
    write(outfile, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")

    # Read genome
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)

    write_lock = ReentrantLock()
    max_range = length_noPAM(motif)
    g_count = length(guides)

    # Early stopping state
    is_es = falses(g_count)
    es_accumulator = zeros(Int, g_count, distance + 1)
    all_offt = [Vector{Offtarget}() for _ in 1:g_count]
    all_offt_lock = ReentrantLock()

    # Process each chromosome
    for (chrom_idx, chrom) in enumerate(dbi.gi.chrom)
        seq = getchromseq(dbi.gi.is_fa, reader[chrom])

        # Locate telomeres - skip leading/trailing Ns
        (seq_start, seq_stop) = locate_telomeres(seq)
        seq = seq[seq_start:seq_stop]
        seq_str = String(seq)

        # Parallel over guides
        ThreadsX.foreach(enumerate(guides)) do (guide_idx, guide)
            if is_es[guide_idx]
                return  # Skip early-stopped guides
            end

            guide_offtargets = Vector{Offtarget}()

            # Forward search
            results_fwd = search_sassy_guide(guide, seq_str, distance, motif, dbi, false, seq_start - 1)
            for r in results_fwd
                loc = Loc(dbi.gi.chrom_type(chrom_idx), dbi.gi.pos_type(r.pos), true)
                offt = Offtarget(loc, r.dist, r.aln_guide, r.aln_ref)

                old_dist, new_dist = insert_offtarget_sassy!(guide_offtargets, offt, max_range)
                if !isnothing(old_dist)
                    es_accumulator[guide_idx, old_dist + 1] -= 1
                end
                if !isnothing(new_dist)
                    es_accumulator[guide_idx, new_dist + 1] += 1
                    if es_accumulator[guide_idx, new_dist + 1] >= early_stopping[new_dist + 1]
                        is_es[guide_idx] = true
                        break  # Early stop
                    end
                end
            end

            # Early stop if triggered
            if is_es[guide_idx]
                return
            end

            # Reverse complement search
            guide_rc = reverse_complement(guide)
            results_rc = search_sassy_guide(guide_rc, seq_str, distance, motif, dbi, true, seq_start - 1)
            for r in results_rc
                loc = Loc(dbi.gi.chrom_type(chrom_idx), dbi.gi.pos_type(r.pos), false)
                offt = Offtarget(loc, r.dist, r.aln_guide, r.aln_ref)

                old_dist, new_dist = insert_offtarget_sassy!(guide_offtargets, offt, max_range)
                if !isnothing(old_dist)
                    es_accumulator[guide_idx, old_dist + 1] -= 1
                end
                if !isnothing(new_dist)
                    es_accumulator[guide_idx, new_dist + 1] += 1
                    if es_accumulator[guide_idx, new_dist + 1] >= early_stopping[new_dist + 1]
                        is_es[guide_idx] = true
                        break  # Early stop
                    end
                end
            end

            # Store results for this guide
            lock(all_offt_lock) do
                append!(all_offt[guide_idx], guide_offtargets)
            end
        end
    end

    # Write all results
    for i in 1:g_count
        for offt in all_offt[i]
            # Apply orientation transformations to match linearDB output
            # For extends5 (Cas9): reverse the alignment strings
            # For extends5=false (Cas12a): reverse complement antisense alignments
            if motif.extends5
                aln_guide = string(reverse(offt.aln_guide))
                aln_ref = string(reverse(offt.aln_ref))
            elseif !offt.loc.isplus
                # Antisense with extends5=false: RC back to original guide orientation
                aln_guide = string(reverse_complement(LongDNA{4}(offt.aln_guide)))
                aln_ref = string(reverse_complement(LongDNA{4}(offt.aln_ref)))
            else
                aln_guide = offt.aln_guide
                aln_ref = offt.aln_ref
            end

            noloc = string(guides[i]) * "," * aln_guide * "," *
                    aln_ref * "," * string(offt.dist) * ","
            write(outfile, noloc * decode(offt.loc, dbi) * "\n")
        end
    end

    close(ref)
    close(outfile)
end

"""
    search_sassy_guide(guide_seq, genome_str, k, motif, dbi, is_antisense, telomere_offset)

Search for one guide in genome.
Uses SIMD for fast candidate filtering, then verifies with align().
Returns Vector{SassyMatch}.

# Arguments
- `telomere_offset`: Number of bases skipped at start of chromosome (for position correction)
"""
function search_sassy_guide(
    guide_seq::LongDNA{4},
    genome_str::String,
    k::Int,
    motif::Motif,
    dbi::DBInfo,
    is_antisense::Bool,
    telomere_offset::Int = 0  # Offset from skipped telomeres at chromosome start
)
    # Build full pattern with PAM
    # For antisense, PAM is at the END (because it's on reverse strand)
    # For forward, PAM is at the START (for extends5=false like Cas12a)
    if motif.extends5
        full_pattern_seq = appendPAM_forward(guide_seq, motif)
        pam_loci = is_antisense ? motif.pam_loci_rve : motif.pam_loci_fwd
    else
        # For extends5=false, PAM is at 5' end
        # Forward: PAM at start of pattern
        # Antisense: PAM at end of pattern (because it's RC on reverse strand)
        if is_antisense
            # PAM at end: guide + rve_pam
            full_pattern_seq = guide_seq * motif.rve[motif.pam_loci_rve]
            pam_loci = (length(guide_seq) + 1):(length(guide_seq) + length(motif.pam_loci_rve))
        else
            full_pattern_seq = appendPAM_forward(guide_seq, motif)
            pam_loci = motif.pam_loci_fwd
        end
    end

    pattern_bytes = Vector{UInt8}(String(full_pattern_seq))
    m = length(pattern_bytes)

    # Encode pattern for SIMD
    (bases, pattern_indices) = encode_pattern_sassy(pattern_bytes)

    pam_len = length(pam_loci)

    # Convert genome to bytes
    genome_bytes = Vector{UInt8}(genome_str)
    n = length(genome_bytes)

    results = SassyMatch[]

    # Overlap size for chunk boundaries - must be at least m + k
    # This ensures we don't miss matches that span chunk boundaries
    overlap = m + k + 1  # Extra +1 for safety at boundaries

    # Process genome in overlapping chunks
    chunk_start = 1
    while chunk_start <= n
        chunk_end = min(chunk_start + CHUNK_SIZE - 1, n)
        chunk = view(genome_bytes, chunk_start:chunk_end)

        # Use SIMD for fast candidate filtering
        chunk_matches = search_sassy_simd_chunk(pattern_indices, chunk, k, bases)

        # Verify each candidate with exact alignment
        for (local_pos, _score) in chunk_matches
            global_end = chunk_start + local_pos - 1
            global_start = max(1, global_end - m + 1)

            if global_end <= n && global_start >= 1
                # Extract and verify PAM
                # For antisense, we need to calculate PAM position based on genome orientation
                # and pattern construction (where PAM is in the pattern)
                actual_pam_start, actual_pam_end = if is_antisense
                    if motif.extends5
                        # PAM is at end of match (in genome coordinates for antisense)
                        (global_end - pam_len + 1, global_end)
                    else
                        # For antisense with extends5=false, PAM is at END of pattern
                        # Pattern: guide + PAM, so PAM is at end
                        (global_end - pam_len + 1, global_end)
                    end
                else
                    if motif.extends5
                        (global_end - pam_len + 1, global_end)
                    else
                        # Forward with extends5=false, PAM is at start
                        (global_start, global_start + pam_len - 1)
                    end
                end

                ref_pam_str, guide_pam_str = extract_pam_strings(
                    genome_bytes, global_start, global_end,
                    full_pattern_seq, actual_pam_start, actual_pam_end, pam_len, motif.extends5, is_antisense
                )

                if !isnothing(ref_pam_str) && check_pam_match(ref_pam_str, guide_pam_str)
                    # For alignment, we need to align the guide (without PAM) against the reference (without PAM)
                    # The output format expects alignment of the guide only, not including PAM

                    # Extract the guide without PAM from the full pattern
                    # For Cas12a (extends5=false):
                    #   Forward: PAM at start, remove from start
                    #   Antisense: PAM at end, remove from end
                    # For Cas9 (extends5=true): PAM at end, remove from end
                    guide_only = if motif.extends5
                        removepam(full_pattern_seq, length(full_pattern_seq) - length(pam_loci) + 1 : length(full_pattern_seq))
                    elseif is_antisense
                        # Antisense with extends5=false: PAM at end
                        removepam(full_pattern_seq, length(full_pattern_seq) - length(pam_loci) + 1 : length(full_pattern_seq))
                    else
                        # Forward with extends5=false: PAM at start
                        removepam(full_pattern_seq, 1 : length(pam_loci))
                    end

                    # Extract reference match without PAM
                    ref_match_seq = LongDNA{4}(genome_str[global_start:global_end])
                    if is_antisense && !motif.extends5
                        # Antisense Cas12a: PAM at end of pattern (and end of reference match)
                        ref_spacer_only = removepam(ref_match_seq, length(ref_match_seq) - length(pam_loci) + 1 : length(ref_match_seq))
                    else
                        # PAM at start of reference match
                        ref_spacer_only = removepam(ref_match_seq, 1 : length(pam_loci))
                    end

                    # Skip if lengths don't match (can happen at boundaries)
                    if length(guide_only) != length(ref_spacer_only)
                        continue
                    end

                    # Use CHOPOFF's align() for exact alignment (guide vs reference spacer)
                    aln_res = align(guide_only, ref_spacer_only, k)

                    if aln_res.dist <= k
                        # Check ambiguity filter (on the spacer, not including PAM)
                        if count(!iscertain, ref_spacer_only) <= motif.ambig_max
                            # Position convention matching linearDB
                            # Add telomere_offset to get actual chromosome position
                            pos = calculate_position(global_start, global_end, motif.extends5, is_antisense)
                            pos += telomere_offset
                            push!(results, SassyMatch(pos, aln_res.dist, aln_res.guide, aln_res.ref))
                        end
                    end
                end
            end
        end

        # Move to next chunk with overlap
        chunk_start += CHUNK_SIZE - overlap
        if chunk_start < 1
            chunk_start = 1
        end
    end

    return results
end

"""
    extract_pam_strings(genome_bytes, global_start, global_end, full_pattern, pam_start, pam_end, pam_len, extends5, is_antisense)

Extract PAM strings from reference and pattern for verification.
Returns (ref_pam, guide_pam) or (nothing, nothing) if out of bounds.
"""
function extract_pam_strings(
    genome_bytes::Vector{UInt8},
    global_start::Int,
    global_end::Int,
    full_pattern::LongDNA{4},
    pam_start::Int,
    pam_end::Int,
    pam_len::Int,
    extends5::Bool,
    is_antisense::Bool
)
    n = length(genome_bytes)

    if pam_start < 1 || pam_end > n
        return (nothing, nothing)
    end

    ref_pam = String(genome_bytes[pam_start:pam_end])

    # Get guide PAM from pattern
    if extends5
        # PAM is at end of pattern
        guide_pam = String(full_pattern[end-pam_len+1:end])
    elseif is_antisense
        # Antisense with extends5=false: PAM is at end of pattern
        guide_pam = String(full_pattern[end-pam_len+1:end])
    else
        # Forward with extends5=false: PAM is at start of pattern
        guide_pam = String(full_pattern[1:pam_len])
    end

    return (ref_pam, guide_pam)
end

"""
    calculate_position(global_start, global_end, extends5, is_antisense)

Calculate position following linearDB convention.
- extends5=true (Cas9): antisense uses global_start, sense uses global_end
- extends5=false (Cas12a): antisense uses global_end, sense uses global_start
"""
@inline function calculate_position(global_start::Int, global_end::Int, extends5::Bool, is_antisense::Bool)
    if extends5
        return is_antisense ? global_start : global_end
    else
        return is_antisense ? global_end : global_start
    end
end

# --- CLI Compatibility ---

# Signature for direct CLI call (for storage_dir compatibility - not used for sassy)
function search_sassy(
    storage_dir::String,
    guides::Vector{LongDNA{4}},
    output_file::String;
    distance::Int = 4
)
    error("Direct search_sassy requires genome_path and motif. Use: search_sassy(guides, genome_path, motif, output_path; distance=4)")
end
