function seq_to_bytes!(buf::Vector{UInt8}, seq::LongDNA{4})
    n = length(seq)
    resize!(buf, n)
    @inbounds for i in 1:n
        nibble = BioSequences.extract_encoded_element(seq, i)
        buf[i] = NIBBLE_TO_ASCII[nibble + 1]
    end
    return buf
end

struct SassyMatch
    pos::Int
    dist::Int
    aln_guide::String
    aln_ref::String
end

const TB_MATCH = UInt8(0x00)
const TB_INS = UInt8(0x01)
const TB_SUB = UInt8(0x02)
const TB_DEL = UInt8(0x03)
const TB_INVALID = UInt8(0xff)
const TB_GAP = UInt8('-')

@inline function tb_idx(i::Int, j::Int, ncols::Int)
    return i * ncols + j + 1
end

@inline function validate_traceback_backend(traceback_backend::Symbol)
    if traceback_backend === :align || traceback_backend === :custom
        return traceback_backend
    end
    throw(
        ArgumentError(
            "Invalid traceback_backend=$(repr(traceback_backend)). Allowed values are :align or :custom.",
        ),
    )
end

@inline function traceback_custom(
    guide_for_aln::LongDNA{4},
    ref_for_aln::LongDNA{4},
    k::Int,
)
    # Fast rejection for non-hits; full traceback is only required for accepted candidates.
    quick_dist = levenshtein(guide_for_aln, ref_for_aln, k)
    quick_dist > k && return Aln("", "", k + 1)

    len1 = Base.length(guide_for_aln)
    len2 = Base.length(ref_for_aln)

    # Strip common prefix, matching CHOPOFF.align semiglobal semantics.
    p = 0
    min_len = min(len1, len2)
    @inbounds for i in 1:min_len
        iscompatible(guide_for_aln[i], ref_for_aln[i]) || break
        p += 1
    end
    if p == len1
        return Aln(String(guide_for_aln), String(ref_for_aln[1:len1]), 0)
    end

    g = guide_for_aln[(p + 1):len1]
    r = ref_for_aln[(p + 1):len2]
    glen = Base.length(g)
    rlen = Base.length(r)
    k_eff = min(k, glen, rlen)

    ncols = glen + 1
    inf_cost = Int16(typemax(Int16))
    prev = Vector{Int16}(undef, ncols)
    cur = Vector{Int16}(undef, ncols)
    parent = fill(TB_INVALID, (rlen + 1) * ncols)

    @inbounds for j in 0:glen
        prev[j + 1] = Int16(j)
    end
    @inbounds for j in 1:glen
        parent[tb_idx(0, j, ncols)] = TB_DEL
    end
    @inbounds for i in 1:rlen
        parent[tb_idx(i, 0, ncols)] = TB_INS
    end

    j_start_max = max(1, glen - k_eff)
    j_end = k_eff

    best_cost = Int(prev[end])
    best_i = 0

    @inbounds for i in 1:rlen
        fill!(cur, inf_cost)
        row_min = glen

        j_start = max(1, i - k_eff)
        if j_start > j_start_max
            j_start = j_start_max
        end
        if j_end < glen
            j_end += 1
        end

        # Left boundary for the active band: only insertion path is valid here.
        cur[j_start] = Int16(i)
        parent[tb_idx(i, j_start - 1, ncols)] = TB_INS
        row_min = min(row_min, i)

        r_ch = r[i]
        for j in j_start:j_end
            diag = Int(prev[j])
            up = Int(prev[j + 1])
            left = Int(cur[j])

            is_match = iscompatible(r_ch, g[j])
            diag_cost = diag + (is_match ? 0 : 1)
            up_cost = up + 1
            left_cost = left + 1

            min_cost = diag_cost
            op = is_match ? TB_MATCH : TB_SUB

            # Rust-style deterministic tie-break order:
            # match > ins > sub > del
            if is_match
                if up_cost < min_cost
                    min_cost = up_cost
                    op = TB_INS
                end
                if left_cost < min_cost
                    min_cost = left_cost
                    op = TB_DEL
                end
            else
                if up_cost <= min_cost
                    min_cost = up_cost
                    op = TB_INS
                end
                if diag_cost < min_cost
                    min_cost = diag_cost
                    op = TB_SUB
                end
                if left_cost < min_cost
                    min_cost = left_cost
                    op = TB_DEL
                end
            end

            cur[j + 1] = Int16(min_cost)
            parent[tb_idx(i, j, ncols)] = op
            row_min = min(row_min, min_cost)

            if j == glen && min_cost < best_cost
                best_cost = min_cost
                best_i = i
            end
        end

        row_min > k_eff && break
        prev, cur = cur, prev
    end

    if best_cost > k_eff
        return Aln("", "", k_eff + 1)
    end

    g_bytes = codeunits(String(g))
    r_bytes = codeunits(String(r))

    aln_g = UInt8[]
    aln_r = UInt8[]
    sizehint!(aln_g, glen + best_i + 8)
    sizehint!(aln_r, glen + best_i + 8)

    i = best_i
    j = glen
    while i > 0 || j > 0
        op = if j == 0
            TB_INS
        elseif i == 0
            TB_DEL
        else
            parent[tb_idx(i, j, ncols)]
        end

        if op == TB_MATCH || op == TB_SUB
            push!(aln_g, g_bytes[j])
            push!(aln_r, r_bytes[i])
            i -= 1
            j -= 1
        elseif op == TB_INS
            push!(aln_g, TB_GAP)
            push!(aln_r, r_bytes[i])
            i -= 1
        elseif op == TB_DEL
            push!(aln_g, g_bytes[j])
            push!(aln_r, TB_GAP)
            j -= 1
        else
            error("Invalid traceback op encountered: $op at (i=$i, j=$j)")
        end
    end

    reverse!(aln_g)
    reverse!(aln_r)

    p_g = p > 0 ? codeunits(String(guide_for_aln[1:p])) : UInt8[]
    p_r = p > 0 ? codeunits(String(ref_for_aln[1:p])) : UInt8[]

    out_g = Vector{UInt8}(undef, Base.length(p_g) + Base.length(aln_g))
    out_r = Vector{UInt8}(undef, Base.length(p_r) + Base.length(aln_r))

    if Base.length(p_g) > 0
        copyto!(out_g, 1, p_g, 1, Base.length(p_g))
        copyto!(out_r, 1, p_r, 1, Base.length(p_r))
    end
    copyto!(out_g, Base.length(p_g) + 1, aln_g, 1, Base.length(aln_g))
    copyto!(out_r, Base.length(p_r) + 1, aln_r, 1, Base.length(aln_r))

    return Aln(String(out_g), String(out_r), best_cost)
end

@inline function run_traceback(
    guide_for_aln::LongDNA{4},
    ref_for_aln::LongDNA{4},
    k::Int,
    traceback_backend::Symbol,
)
    if traceback_backend === :custom
        return traceback_custom(guide_for_aln, ref_for_aln, k)
    end
    return align(guide_for_aln, ref_for_aln, k)
end

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

@inline function ref_ambig_within_limit(
    genome_bytes::AbstractVector{UInt8},
    start_pos::Int,
    end_pos::Int,
    limit::Int,
)
    limit >= end_pos - start_pos + 1 && return true
    n_ambig = 0
    @inbounds for i in start_pos:end_pos
        n_ambig += Sassy.REF_AMBIG_TABLE[genome_bytes[i] + 1]
        n_ambig > limit && return false
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
    genome_bytes::Vector{UInt8},
    genome_seq::LongDNA{4},
    k::Int,
    motif::Motif,
    dbi::DBInfo,
    is_antisense::Bool,
    telomere_offset::Int = 0,
    impl_func::F = (idx, txt, k, b; workspace=nothing) -> search_sassy_impl(idx, txt, k, b, Val(4), Val(true); workspace);
    strict_pam::Bool = true,
    traceback_backend::Symbol = :custom,
    workspace::Union{SassyWorkspace, Nothing} = nothing,
) where {F <: Function}
    traceback_backend = validate_traceback_backend(traceback_backend)
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

    search_pattern = full_pattern_seq
    pattern_bytes = Vector{UInt8}(String(search_pattern))
    search_len = Base.length(pattern_bytes)
    guide_len = Base.length(guide_pattern)
    pam_len = Base.length(pam_str)
    (bases, pattern_indices) = encode_pattern_sassy(pattern_bytes)

    n = Base.length(genome_bytes)
    pam_on_left = pam_at_start(motif, is_antisense)
    results = SassyMatch[]
    seen = Set{Tuple{Int, Int, String, String}}()
    query_for_linear_aln = motif.extends5 ? reverse(guide_seq) : guide_seq

    # 3. Run Sassy Algorithm
    matches = impl_func(pattern_indices, genome_bytes, k, bases; workspace)

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
                false,
            )
            if motif_start < 1 || motif_end > n || guide_start < 1 || guide_end > n
                continue
            end

            if !ref_ambig_within_limit(genome_bytes, motif_start, motif_end, motif.ambig_max)
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
                run_traceback(query_for_linear_aln, ref_for_aln, k, traceback_backend)
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
                    ref_for_aln = genome_seq[win_start:win_end]
                    guide_for_aln = search_pattern
                else
                    win_start = max(1, align_start - k)
                    win_end = trial_end
                    if win_start > win_end
                        continue
                    end
                    ref_for_aln = reverse(genome_seq[win_start:win_end])
                    guide_for_aln = reverse(search_pattern)
                    reverse_back = true
                end

                aln_tmp = run_traceback(guide_for_aln, ref_for_aln, k, traceback_backend)
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

"""
    search_sassy(guides, genome_path, motif, output_path; kwargs...)

Search off-targets directly on a genome sequence using the SASSY reimplementation.

Keyword arguments:

- `distance::Int=4`: maximum edit distance.
- `early_stopping::Union{Vector{Int}, Nothing}=nothing`: per-distance early-stopping
  thresholds.
- `use_avx512::Bool=false`: use 8-lane AVX-512 kernel when available.
- `force_safe_minima::Bool=false`: disable BMI2/PEXT minima path and force safe minima.
- `strict_pam::Bool=true`: enforce strict PAM matching at candidate coordinates.
- `traceback_backend::Symbol=:custom`: traceback backend (`:custom` or `:align`).

The output CSV columns are:
`guide,alignment_guide,alignment_reference,distance,chromosome,start,strand`.
"""
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
    traceback_backend::Symbol = :custom,
)
    traceback_backend = validate_traceback_backend(traceback_backend)
    # Dispatch Logic based on flags
    # We create a function barrier or closure to avoid dynamic dispatch in the loop
    # Default to PEXT on BMI2-capable x86; otherwise use safe minima scanning.
    use_pext = !force_safe_minima && can_use_bmi2_pext()
    
    _search_impl = if use_avx512
        if use_pext
            (idx, txt, k, b; workspace=nothing) -> search_sassy_impl(idx, txt, k, b, Val(8), Val(true); workspace)
        else
            (idx, txt, k, b; workspace=nothing) -> search_sassy_impl(idx, txt, k, b, Val(8), Val(false); workspace)
        end
    else
        if use_pext
            (idx, txt, k, b; workspace=nothing) -> search_sassy_impl(idx, txt, k, b, Val(4), Val(true); workspace)
        else
            (idx, txt, k, b; workspace=nothing) -> search_sassy_impl(idx, txt, k, b, Val(4), Val(false); workspace)
        end
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
    is_es = fill(false, g_count)
    es_accumulator = zeros(Int, g_count, distance + 1)
    all_offt = [Vector{Offtarget}() for _ in 1:g_count]

    # Pre-allocate one workspace per guide (reused across chromosomes)
    m_hint = Base.length(String(motif.fwd))
    lanes_val = use_avx512 ? Val(8) : Val(4)
    workspaces = [SassyWorkspace(m_hint, lanes_val) for _ in 1:g_count]

    genome_buf = Vector{UInt8}(undef, 0)

    for (chrom_idx, chrom) in enumerate(dbi.gi.chrom)
        seq = getchromseq(dbi.gi.is_fa, reader[chrom])
        (seq_start, seq_stop) = locate_telomeres(seq)
        seq = seq[seq_start:seq_stop]
        seq_to_bytes!(genome_buf, seq)

        ThreadsX.foreach(enumerate(guides)) do (guide_idx, guide)
            if is_es[guide_idx]; return; end
            ws = workspaces[guide_idx]

            guide_offtargets = Vector{Offtarget}()

            function process_results!(results, is_plus)
                for r in results
                    # Convert internal SassyMatch to Offtarget
                    loc = Loc(dbi.gi.chrom_type(chrom_idx), dbi.gi.pos_type(r.pos), is_plus)
                    offt = Offtarget(loc, r.dist, r.aln_guide, r.aln_ref)

                    if use_es
                        dist_idx = offt.dist + 1
                        if es_accumulator[guide_idx, dist_idx] >= es_limits[dist_idx]
                            continue  # skip — limit already reached for this distance
                        end
                        es_accumulator[guide_idx, dist_idx] += 1
                        if es_accumulator[guide_idx, dist_idx] >= es_limits[dist_idx]
                            is_es[guide_idx] = true
                        end
                    end

                    push!(all_offt[guide_idx], offt)
                end
            end

            # Forward Search
            results_fwd = search_sassy_guide(
                guide, genome_buf, seq, distance, motif, dbi, false, seq_start - 1, _search_impl;
                strict_pam = strict_pam,
                traceback_backend = traceback_backend,
                workspace = ws,
            )
            process_results!(results_fwd, true)

            # Reverse Complement Search (only if not early stopped)
            if !is_es[guide_idx]
                results_rc = search_sassy_guide(
                    guide, genome_buf, seq, distance, motif, dbi, true, seq_start - 1, _search_impl;
                    strict_pam = strict_pam,
                    traceback_backend = traceback_backend,
                    workspace = ws,
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

            print(outfile, guides[i], ",", aln_guide, ",", aln_ref, ",", offt.dist, ",", decode(offt.loc, dbi), "\n")
        end
    end

    close(ref)
    close(outfile)
end
