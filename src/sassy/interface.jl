struct SassyMatch
    pos::Int
    dist::Int
    aln_guide::String
    aln_ref::String
end

const LINEAR_TRACE_PREFIX_LEN = 7

function check_pam_match(ref_pam::AbstractString, guide_pam::AbstractString)
    Base.length(ref_pam) == Base.length(guide_pam) || return false
    for (r, g) in zip(ref_pam, guide_pam)
        m1 = Sassy.get_iupac_mask(UInt8(r))
        m2 = Sassy.get_iupac_mask(UInt8(g))
        (m1 & m2) != 0 || return false
    end
    return true
end

@inline function pam_at_start(motif::Motif, is_antisense::Bool)
    # Cas9: antisense has PAM on the left (CCN + guide_rc)
    # Cas12a: forward has PAM on the left (TTTV + guide)
    return motif.extends5 == is_antisense
end

function pam_matches_candidate(
    genome_bytes::AbstractVector{UInt8},
    guide_end::Int,
    guide_len::Int,
    pam_seq::AbstractString,
    pam_on_left::Bool,
)
    pam_len = Base.length(pam_seq)
    pam_start = if pam_on_left
        guide_start = guide_end - guide_len + 1
        guide_start - pam_len
    else
        guide_end + 1
    end
    pam_stop = pam_start + pam_len - 1
    if pam_start < 1 || pam_stop > Base.length(genome_bytes)
        return false
    end
    for i in 1:pam_len
        ref_mask = Sassy.get_iupac_mask(@inbounds genome_bytes[pam_start + i - 1])
        pam_mask = Sassy.get_iupac_mask(UInt8(@inbounds pam_seq[i]))
        if (ref_mask & pam_mask) == 0
            return false
        end
    end
    return true
end

@inline function candidate_bounds(
    match_end::Int,
    guide_len::Int,
    pam_len::Int,
    pam_on_left::Bool,
    strict_pam::Bool,
)
    if strict_pam
        guide_end = match_end
        guide_start = guide_end - guide_len + 1
        motif_start = pam_on_left ? (guide_start - pam_len) : guide_start
        motif_end = pam_on_left ? guide_end : (guide_end + pam_len)
    else
        motif_end = match_end
        motif_start = motif_end - (guide_len + pam_len) + 1
        if pam_on_left
            guide_start = motif_start + pam_len
            guide_end = motif_end
        else
            guide_start = motif_start
            guide_end = motif_end - pam_len
        end
    end
    return motif_start, motif_end, guide_start, guide_end
end

@inline function extends_right(motif::Motif, is_antisense::Bool)
    # Mirrors CHOPOFF.add_extension dispatch.
    return (motif.extends5 && is_antisense) || (!motif.extends5 && !is_antisense)
end

@inline function normalize_ref_for_linear(
    seq::LongDNA{4},
    motif::Motif,
    is_antisense::Bool,
)
    if motif.extends5 && is_antisense
        return complement(seq)
    elseif motif.extends5 && !is_antisense
        return reverse(seq)
    elseif !motif.extends5 && is_antisense
        return reverse_complement(seq)
    else
        return seq
    end
end

function search_sassy_guide(
    guide_seq::LongDNA{4},
    genome_str::String,
    k::Int,
    motif::Motif,
    dbi::DBInfo,
    is_antisense::Bool,
    telomere_offset::Int = 0,
    impl_func::Function = (idx, txt, k, b) -> search_sassy_impl(idx, txt, k, b, Val(4), Val(true));
    strict_pam::Bool = true,
)
    # 1. Prepare guide orientation and PAM pattern for this strand.
    guide_pattern = is_antisense ? reverse_complement(guide_seq) : guide_seq
    pam_str = ""
    if motif.extends5 # Cas9
        if !is_antisense
            pam_str = String(motif.fwd[motif.pam_loci_fwd])
            full_pattern_seq = guide_pattern * LongDNA{4}(pam_str)
        else
            pam_str = String(motif.rve[motif.pam_loci_rve])
            full_pattern_seq = LongDNA{4}(pam_str) * guide_pattern
        end
    else # Cas12a
        if !is_antisense
            pam_str = String(motif.fwd[motif.pam_loci_fwd])
            full_pattern_seq = LongDNA{4}(pam_str) * guide_pattern
        else
            pam_str = String(motif.rve[motif.pam_loci_rve])
            full_pattern_seq = guide_pattern * LongDNA{4}(pam_str)
        end
    end

    search_pattern = strict_pam ? guide_pattern : full_pattern_seq
    pattern_bytes = Vector{UInt8}(String(search_pattern))
    search_len = Base.length(pattern_bytes)
    guide_len = Base.length(guide_pattern)
    pam_len = Base.length(pam_str)
    (bases, pattern_indices) = encode_pattern_sassy(pattern_bytes)

    genome_bytes = Vector{UInt8}(genome_str)
    genome_seq = LongDNA{4}(genome_str)
    n = Base.length(genome_bytes)
    pam_on_left = pam_at_start(motif, is_antisense)
    results = SassyMatch[]
    seen = Set{Tuple{Int, Int, String, String}}()
    query_for_linear_aln = motif.extends5 ? reverse(guide_seq) : guide_seq

    # 3. Run Sassy Algorithm
    matches = impl_func(pattern_indices, genome_bytes, k, bases)

    # 4. Filter and Verify (Traceback/Alignment)
    for (local_pos, _score) in matches
        global_reported_end = local_pos

        # Search around the reported position and emit all valid alignments <= k.
        for shift in -k:k
            trial_end = global_reported_end + shift
            if trial_end < 1 || trial_end > n; continue; end

            motif_start, motif_end, guide_start, guide_end = candidate_bounds(
                trial_end,
                guide_len,
                pam_len,
                pam_on_left,
                strict_pam,
            )
            if motif_start < 1 || motif_end > n || guide_start < 1 || guide_end > n
                continue
            end

            if strict_pam && !pam_matches_candidate(genome_bytes, guide_end, guide_len, pam_str, pam_on_left)
                continue
            end

            aln_res = if strict_pam
                guide_slice = LongDNA{4}(genome_seq[guide_start:guide_end])
                if extends_right(motif, is_antisense)
                    ext = getExt3(genome_seq, n, guide_end + 1, motif.distance)
                    ref_with_ext = guide_slice * ext
                else
                    ext = getExt5(genome_seq, guide_start - 1, motif.distance)
                    ref_with_ext = ext * guide_slice
                end
                ref_for_aln = normalize_ref_for_linear(ref_with_ext, motif, is_antisense)
                prefix_len = min(LINEAR_TRACE_PREFIX_LEN, Base.length(ref_for_aln))
                prefix = ref_for_aln[1:prefix_len]
                suffix = if prefix_len < Base.length(ref_for_aln)
                    ref_for_aln[(prefix_len + 1):end]
                else
                    LongDNA{4}("")
                end
                pa = prefix_align(query_for_linear_aln, prefix, Base.length(suffix), k)
                suffix_align(suffix, pa)
            else
                align_start = trial_end - search_len + 1
                if align_start < 1
                    continue
                end

                reverse_back = false
                if pam_on_left
                    win_start = align_start
                    win_end = min(n, trial_end + k)
                    if win_start > win_end
                        continue
                    end
                    ref_for_aln = LongDNA{4}(genome_str[win_start:win_end])
                    guide_for_aln = search_pattern
                else
                    win_start = max(1, align_start - k)
                    win_end = trial_end
                    if win_start > win_end
                        continue
                    end
                    ref_for_aln = reverse(LongDNA{4}(genome_str[win_start:win_end]))
                    guide_for_aln = reverse(search_pattern)
                    reverse_back = true
                end

                aln_tmp = align(guide_for_aln, ref_for_aln, k)
                if reverse_back
                    Aln(reverse(aln_tmp.guide), reverse(aln_tmp.ref), aln_tmp.dist)
                else
                    aln_tmp
                end
            end
            
            if aln_res.dist <= k
                aln_guide = aln_res.guide
                aln_ref = aln_res.ref
                if strict_pam && motif.extends5
                    # Mirrors search_linearDB output rendering for extends5 motifs.
                    aln_guide = reverse(aln_guide)
                    aln_ref = reverse(aln_ref)
                end

                # Match linearDB location normalization (PAM-side anchor convention).
                pos_local = (motif.extends5 == is_antisense) ? motif_start : motif_end
                pos = pos_local + telomere_offset

                key = (pos, aln_res.dist, aln_guide, aln_ref)
                if !(key in seen)
                    push!(results, SassyMatch(pos, aln_res.dist, aln_guide, aln_ref))
                    push!(seen, key)
                end
            end
        end
    end
    return results
end

# --- Main Entry Point (Batched Database Search) ---

function search_sassy(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif,
    output_path::String;
    distance::Int = 4,
    early_stopping::Union{Vector{Int}, Nothing} = nothing,
    use_avx512::Bool = false,
    force_safe_minima::Bool = false,
    strict_pam::Bool = true,
)
    # Dispatch Logic based on flags
    # We create a function barrier or closure to avoid dynamic dispatch in the loop
    # Determine PEXT usage (Default: true, unless forced safe)
    use_pext = !force_safe_minima
    
    _search_impl = if use_avx512
        (idx, txt, k, b) -> search_sassy_impl(idx, txt, k, b, Val(8), Val(use_pext))
    else
        (idx, txt, k, b) -> search_sassy_impl(idx, txt, k, b, Val(4), Val(use_pext))
    end
    dbi = DBInfo(genome_path, "sassy_search", motif)
    if any(length_noPAM(motif) .!= length.(guides))
        error("Guide queries are not of the correct length to use with this Motif.")
    end

    use_es = early_stopping !== nothing
    es_limits = use_es ? early_stopping : Int[]

    mkpath(dirname(output_path))
    outfile = open(output_path, "w")
    write(outfile, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")

    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)

    g_count = Base.length(guides)
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
                    # Convert internal SassyMatch to Offtarget
                    loc = Loc(dbi.gi.chrom_type(chrom_idx), dbi.gi.pos_type(r.pos), is_plus)
                    offt = Offtarget(loc, r.dist, r.aln_guide, r.aln_ref)

                    # Keep all valid alignments (linearDB parity), only lock for shared state.
                    lock(all_offt_lock) do
                        push!(all_offt[guide_idx], offt)

                        if use_es
                            es_accumulator[guide_idx, offt.dist + 1] += 1
                            if es_accumulator[guide_idx, offt.dist + 1] >= es_limits[offt.dist + 1]
                                is_es[guide_idx] = true
                            end
                        end
                    end
                end
            end

            # Forward Search
            results_fwd = search_sassy_guide(
                guide, seq_str, distance, motif, dbi, false, seq_start - 1, _search_impl;
                strict_pam = strict_pam,
            )
            process_results!(results_fwd, true)

            # Reverse Complement Search (only if not early stopped)
            if !is_es[guide_idx]
                results_rc = search_sassy_guide(
                    guide, seq_str, distance, motif, dbi, true, seq_start - 1, _search_impl;
                    strict_pam = strict_pam,
                )
                process_results!(results_rc, false)
            end
        end
    end

    # Write output
    for i in 1:g_count
        for offt in all_offt[i]
            # Adjust alignment strings based on orientation/motif for display
            # Alignments in all_offt are already oriented correctly
            aln_guide = offt.aln_guide
            aln_ref = offt.aln_ref

            noloc = string(guides[i]) * "," * aln_guide * "," *
                    aln_ref * "," * string(offt.dist) * ","
            write(outfile, noloc * decode(offt.loc, dbi) * "\n")
        end
    end

    close(ref)
    close(outfile)
end
