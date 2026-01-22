# Sassy SIMD-Optimized Implementation
# Port of Rust Sassy's transposed SIMD-Myers algorithm

using SIMD: Vec, vload, vstore
using BioSequences: reverse_complement, complement

# --- Constants ---
const LANES = 4
const BLOCK_SIZE = 64
const CHUNK_SIZE = LANES * BLOCK_SIZE

# IUPAC encoding
const IUPAC_SASSY = let
    enc = zeros(UInt8, 32)
    A = UInt8(1 << 0)
    C = UInt8(1 << 1)
    T = UInt8(1 << 2)
    G = UInt8(1 << 3)
    enc[Int('A') & 0x1F + 1] = A
    enc[Int('C') & 0x1F + 1] = C
    enc[Int('T') & 0x1F + 1] = T
    enc[Int('U') & 0x1F + 1] = T
    enc[Int('G') & 0x1F + 1] = G
    enc[Int('N') & 0x1F + 1] = A | C | T | G
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
    enc[Int('X') & 0x1F + 1] = 0
    enc
end

@inline get_iupac_mask(c::UInt8) = @inbounds IUPAC_SASSY[(c & 0x1F) + 1]

struct TextBlockProfile
    masks::Vector{UInt64}
end

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

function search_sassy_simd_chunk(
    pattern_indices::Vector{Int},
    text_chunk::AbstractVector{UInt8},
    k::Int,
    bases::Vector{UInt8}
)
    m = length(pattern_indices)
    chunk_len = length(text_chunk)
    padded = fill(UInt8('X'), CHUNK_SIZE)
    copy_len = min(chunk_len, CHUNK_SIZE)
    padded[1:copy_len] .= text_chunk[1:copy_len]

    lane_profiles = [
        encode_text_block(view(padded, (i-1)*BLOCK_SIZE+1 : i*BLOCK_SIZE), bases)
        for i in 1:LANES
    ]

    # Initialize state
    # vp initialized to all 1s (D[i, -1] = i)
    vp = [~UInt64(0) for _ in 1:LANES]
    vm = [UInt64(0) for _ in 1:LANES]
    
    # Last row outputs for score reconstruction
    last_row_hp = zeros(UInt64, LANES)
    last_row_hm = zeros(UInt64, LANES)
    
    for j in 1:m
        char_idx = pattern_indices[j]
        
        # Carry into the first block of this chunk.
        # Boundary condition for semi-global: D[0, j] = 0.
        # This implies HP carry (D[j, 0] - D[j, -1]) should be related.
        # Rust initializes HP carries to 1.
        carry_hp = UInt64(1)
        carry_hm = UInt64(0)
        
        for lane in 1:LANES
            eq_mask = lane_profiles[lane].masks[char_idx]
            
            # 1. Update Eq with Carry HM
            eq_c = eq_mask | carry_hm
            
            # 2. Compute HX (Horizontal X)
            # (Eq & VP) + VP. NO CARRY IN ADDITION.
            term = (eq_c & vp[lane]) + vp[lane]
            hx = (term ⊻ vp[lane]) | eq_c
            
            # 3. Compute HP/HM (Internal)
            hp_internal = vm[lane] | (~(hx | vp[lane]))
            hm_internal = vp[lane] & hx
            
            # 4. Shift and Inject Carry
            hp_shifted = (hp_internal << 1) | carry_hp
            hm_shifted = (hm_internal << 1) | carry_hm
            
            # 5. Extract Output Carries for next lane
            carry_hp = hp_internal >> 63
            carry_hm = hm_internal >> 63
            
            # 6. Update VP/VM
            # Uses Eq with HM folded in (eq_c)
            vx = eq_c | vm[lane]
            
            vp[lane] = hm_shifted | (~(vx | hp_shifted))
            vm[lane] = hp_shifted & vx
            
            if j == m
                last_row_hp[lane] = hp_internal
                last_row_hm[lane] = hm_internal
            end
        end
    end

    matches = Tuple{Int, Int}[]
    
    # Score Reconstruction
    # D[m, 0] = m.
    current_score = m 
    
    for lane in 1:LANES
        hp_vec = last_row_hp[lane]
        hm_vec = last_row_hm[lane]

        for bit in 0:63
            if (hp_vec & (UInt64(1) << bit)) != 0; current_score += 1; end
            if (hm_vec & (UInt64(1) << bit)) != 0; current_score -= 1; end

            global_pos = (lane - 1) * BLOCK_SIZE + bit + 1
            if global_pos <= chunk_len && current_score <= k
                push!(matches, (global_pos, current_score))
            end
        end
    end
    
    return matches
end

struct SassyMatch
    pos::Int
    dist::Int
    aln_guide::String
    aln_ref::String
end

function check_pam_match(ref_pam::AbstractString, guide_pam::AbstractString)
    length(ref_pam) == length(guide_pam) || return false
    for (r, g) in zip(ref_pam, guide_pam)
        m1 = get_iupac_mask(UInt8(r))
        m2 = get_iupac_mask(UInt8(g))
        (m1 & m2) != 0 || return false
    end
    return true
end

function search_sassy(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif,
    output_path::String;
    distance::Int = 4,
    early_stopping::Union{Vector{Int}, Nothing} = nothing
)
    dbi = DBInfo(genome_path, "sassy_search", motif)
    if any(length_noPAM(motif) .!= length.(guides))
        error("Guide queries are not of the correct length to use with this Motif: " * string(motif))
    end

    use_es = early_stopping !== nothing
    es_limits = use_es ? early_stopping : Int[]

    mkpath(dirname(output_path))
    outfile = open(output_path, "w")
    write(outfile, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")

    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)

    g_count = length(guides)
    is_es = falses(g_count)
    es_accumulator = zeros(Int, g_count, distance + 1)
    all_offt = [Vector{Offtarget}() for _ in 1:g_count]
    all_offt_lock = ReentrantLock()

    for (chrom_idx, chrom) in enumerate(dbi.gi.chrom)
        seq = getchromseq(dbi.gi.is_fa, reader[chrom])
        (seq_start, seq_stop) = locate_telomeres(seq)
        seq = seq[seq_start:seq_stop]
        seq_str = String(seq)

        ThreadsX.foreach(enumerate(guides)) do (guide_idx, guide)
            if is_es[guide_idx]; return; end

            guide_offtargets = Vector{Offtarget}()

            function process_results!(results, is_plus)
                for r in results
                    loc = Loc(dbi.gi.chrom_type(chrom_idx), dbi.gi.pos_type(r.pos), is_plus)
                    offt = Offtarget(loc, r.dist, r.aln_guide, r.aln_ref)
                    push!(guide_offtargets, offt)
                    
                    if use_es && r.dist <= distance
                        idx = r.dist + 1
                        es_accumulator[guide_idx, idx] += 1
                        if es_accumulator[guide_idx, idx] >= es_limits[idx]
                            is_es[guide_idx] = true
                        end
                    end
                end
            end

            results_fwd = search_sassy_guide(guide, seq_str, distance, motif, dbi, false, seq_start - 1)
            process_results!(results_fwd, true)

            if !is_es[guide_idx]
                results_rc = search_sassy_guide(guide, seq_str, distance, motif, dbi, true, seq_start - 1)
                process_results!(results_rc, false)
            end

            lock(all_offt_lock) do
                append!(all_offt[guide_idx], guide_offtargets)
            end
        end
    end

    for i in 1:g_count
        for offt in all_offt[i]
            if motif.extends5 # Cas9
                aln_guide = string(reverse(offt.aln_guide))
                aln_ref = string(reverse(offt.aln_ref))
            elseif !offt.loc.isplus # Cas12a Antisense
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

function search_sassy_guide(
    guide_seq::LongDNA{4},
    genome_str::String,
    k::Int,
    motif::Motif,
    dbi::DBInfo,
    is_antisense::Bool,
    telomere_offset::Int = 0
)
    if motif.extends5
        if !is_antisense
            pam_seq = motif.fwd[motif.pam_loci_fwd]
            full_pattern_seq = guide_seq * pam_seq
            pam_start_idx = length(guide_seq) + 1
            pam_end_idx = length(guide_seq) + length(pam_seq)
        else
            pam_seq = motif.rve[motif.pam_loci_rve] 
            rc_guide = reverse_complement(guide_seq)
            full_pattern_seq = pam_seq * rc_guide
            pam_start_idx = 1
            pam_end_idx = length(pam_seq)
        end
    else
        if !is_antisense
            pam_seq = motif.fwd[motif.pam_loci_fwd]
            full_pattern_seq = pam_seq * guide_seq
            pam_start_idx = 1
            pam_end_idx = length(pam_seq)
        else
            pam_seq = motif.rve[motif.pam_loci_rve]
            rc_guide = reverse_complement(guide_seq)
            full_pattern_seq = rc_guide * pam_seq
            pam_start_idx = length(rc_guide) + 1
            pam_end_idx = length(rc_guide) + length(pam_seq)
        end
    end

    pattern_bytes = Vector{UInt8}(String(full_pattern_seq))
    m = length(pattern_bytes)
    (bases, pattern_indices) = encode_pattern_sassy(pattern_bytes)

    genome_bytes = Vector{UInt8}(genome_str)
    n = length(genome_bytes)
    results = SassyMatch[]
    overlap = m + k + 1

    chunk_start = 1
    while chunk_start <= n
        chunk_end = min(chunk_start + CHUNK_SIZE - 1, n)
        chunk = view(genome_bytes, chunk_start:chunk_end)
        chunk_matches = search_sassy_simd_chunk(pattern_indices, chunk, k, bases)

        for (local_pos, _score) in chunk_matches
            global_end = chunk_start + local_pos - 1
            global_start = max(1, global_end - m + 1)

            if global_end <= n && global_start >= 1
                match_seq = LongDNA{4}(genome_str[global_start:global_end])
                
                if pam_end_idx > length(match_seq) || pam_start_idx < 1; continue; end
                
                ref_pam_str = String(match_seq[pam_start_idx:pam_end_idx])
                if pam_start_idx == 1
                    ref_spacer = match_seq[pam_end_idx+1:end]
                else
                    ref_spacer = match_seq[1:pam_start_idx-1]
                end

                expected_pam_str = String(full_pattern_seq[pam_start_idx:pam_end_idx])
                
                if check_pam_match(ref_pam_str, expected_pam_str)
                    if motif.extends5 # Cas9
                        if !is_antisense
                            guide_for_aln = reverse(guide_seq)
                            ref_for_aln = reverse(ref_spacer)
                        else
                            guide_for_aln = reverse(guide_seq)
                            ref_for_aln = complement(ref_spacer)
                        end
                    else # Cas12a
                        if !is_antisense
                            guide_for_aln = guide_seq
                            ref_for_aln = ref_spacer
                        else
                            guide_for_aln = guide_seq
                            ref_for_aln = reverse_complement(ref_spacer)
                        end
                    end

                    aln_res = align(guide_for_aln, ref_for_aln, k)

                    if aln_res.dist <= k
                        if count(!iscertain, ref_spacer) <= motif.ambig_max
                            pos = 0
                            if motif.extends5
                                if !is_antisense
                                    pos = global_end
                                else
                                    pos = global_start
                                end
                            else
                                if !is_antisense
                                    pos = global_start
                                else
                                    pos = global_end
                                end
                            end
                            pos += telomere_offset
                            push!(results, SassyMatch(pos, aln_res.dist, aln_res.guide, aln_res.ref))
                        end
                    end
                end
            end
        end
        chunk_start += CHUNK_SIZE - overlap
        if chunk_start < 1; chunk_start = 1; end
    end
    return results
end

function search_sassy(
    storage_dir::String,
    guides::Vector{LongDNA{4}},
    output_file::String;
    distance::Int = 4
)
    error("Direct search_sassy requires genome_path and motif.")
end