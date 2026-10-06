# Candidate verification, alignment, and result commit.

@inline function prefix_hash_scan_iupac_mask(base::UInt8)
    base = base & UInt8(0xdf)
    base == UInt8('A') && return UInt8(0x01)
    base == UInt8('C') && return UInt8(0x02)
    base == UInt8('G') && return UInt8(0x04)
    base == UInt8('T') && return UInt8(0x08)
    base == UInt8('R') && return UInt8(0x05)
    base == UInt8('Y') && return UInt8(0x0a)
    base == UInt8('S') && return UInt8(0x06)
    base == UInt8('W') && return UInt8(0x09)
    base == UInt8('K') && return UInt8(0x0c)
    base == UInt8('M') && return UInt8(0x03)
    base == UInt8('B') && return UInt8(0x0e)
    base == UInt8('D') && return UInt8(0x0d)
    base == UInt8('H') && return UInt8(0x0b)
    base == UInt8('V') && return UInt8(0x07)
    base == UInt8('N') && return UInt8(0x0f)
    return UInt8(0)
end

@inline function prefix_hash_scan_complement_mask(mask::UInt8)
    return ((mask & UInt8(0x01)) << 3) |
        ((mask & UInt8(0x02)) << 1) |
        ((mask & UInt8(0x04)) >> 1) |
        ((mask & UInt8(0x08)) >> 3)
end

function build_prefix_hash_scan_myers_profile(guide::LongDNA{4})
    length(guide) <= 64 || error("Myers profile supports guides up to 64 bases.")
    peq = zeros(UInt64, 4)
    @inbounds for idx in eachindex(guide)
        bit = UInt64(1) << (idx - 1)
        guide[idx] == DNA_A && (peq[1] |= bit)
        guide[idx] == DNA_C && (peq[2] |= bit)
        guide[idx] == DNA_G && (peq[3] |= bit)
        guide[idx] == DNA_T && (peq[4] |= bit)
    end
    eq_by_iupac = ntuple(16) do idx
        mask = UInt8(idx - 1)
        eq = UInt64(0)
        mask & UInt8(0x01) != 0 && (eq |= peq[1])
        mask & UInt8(0x02) != 0 && (eq |= peq[2])
        mask & UInt8(0x04) != 0 && (eq |= peq[3])
        mask & UInt8(0x08) != 0 && (eq |= peq[4])
        eq
    end
    return PrefixHashScanMyersProfile(
        eq_by_iupac, UInt8(length(guide)), UInt64(1) << (length(guide) - 1))
end

function build_prefix_hash_scan_myers_profiles(guides::Vector{LongDNA{4}})
    return build_prefix_hash_scan_myers_profile.(guides)
end

@inline function reserve_prefix_hash_scan_detail_hit!(
    es_acc, is_es, guide_idx::Int, dist::Int,
    early_stopping::Vector{Int})

    dist_idx = dist + 1
    if es_acc[guide_idx, dist_idx] >= early_stopping[dist_idx]
        is_es[guide_idx] = true
        return false
    end
    es_acc[guide_idx, dist_idx] += 1
    return true
end

@inline function prefix_hash_scan_raw_myers_distance(
    geometry::PrefixScanGeometry,
    profile::PrefixHashScanMyersProfile,
    raw::AbstractVector{UInt8},
    candidate_start::Int,
    is_antisense::Bool,
    distance::Int)

    pattern_length = Int(profile.length)
    reference_length = pattern_length + distance
    first_scored_prefix = pattern_length - distance
    pv = typemax(UInt64)
    mv = UInt64(0)
    score = pattern_length
    best = distance + 1
    @inbounds for ref_idx in 1:reference_length
        raw_idx = candidate_start + prefix_scan_reference_offset(
            geometry.matcher, is_antisense, ref_idx)
        mask = 1 <= raw_idx <= length(raw) ?
            prefix_hash_scan_iupac_mask(raw[raw_idx]) : UInt8(0)
        is_antisense && (mask = prefix_hash_scan_complement_mask(mask))
        eq = profile.eq_by_iupac[Int(mask) + 1]
        xv = eq | mv
        xh = xor((eq & pv) + pv, pv) | eq
        ph = mv | ~(xh | pv)
        mh = pv & xh
        ph & profile.final_bit != 0 && (score += 1)
        mh & profile.final_bit != 0 && (score -= 1)
        ph = (ph << 1) | UInt64(1)
        mh <<= 1
        pv = mh | ~(xv | ph)
        mv = ph & xv
        ref_idx >= first_scored_prefix && (best = min(best, score))
    end
    return min(best, distance + 1)
end

function materialize_normalized_candidate(
    chrom_seq::LongDNA{4},
    candidate_range::UnitRange{Int64},
    dbi::DBInfo,
    is_antisense::Bool)

    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd
    ot = removepam(chrom_seq[candidate_range], pam_loci)
    if dbi.motif.distance > 0
        ots = add_extension([ot], [candidate_range], dbi, chrom_seq, is_antisense)
        ot = ots[1]
    end
    ots, norm_pos = normalize_to_PAMseqEXT([ot], [candidate_range], dbi, is_antisense)
    return ots[1], norm_pos[1]
end

@inline prefix_scan_reference_slice(raw::AbstractVector{UInt8}, range) =
    LongDNA{4}(@view raw[range])

@inline prefix_scan_reference_base(raw::AbstractVector{UInt8}, pos::Int) =
    1 <= pos <= length(raw) ? convert(DNA, Char(raw[pos])) : DNA_Gap

# Guide plus `motif.distance` extension bases, oriented like the guide; bases
# beyond the reference ends are gaps.
function materialize_normalized_candidate_specialized(
    geometry::PrefixScanGeometry,
    source::AbstractVector{UInt8}, candidate_start::Int,
    dbi::DBInfo, is_antisense::Bool)

    matcher = geometry.matcher
    reference_length = geometry.guide_bases + dbi.motif.distance
    first_pos = candidate_start +
        prefix_scan_reference_offset(matcher, is_antisense, 1)
    last_pos = candidate_start +
        prefix_scan_reference_offset(matcher, is_antisense, reference_length)
    low, high = minmax(first_pos, last_pos)
    if prefix_scan_reference_is_affine(matcher, is_antisense) &&
            low >= 1 && high <= length(source)
        ot = prefix_scan_reference_slice(source, low:high)
        first_pos > last_pos && reverse!(ot)
    else
        ot = LongDNA{4}([prefix_scan_reference_base(source, candidate_start +
            prefix_scan_reference_offset(matcher, is_antisense, ref_idx))
            for ref_idx in 1:reference_length])
    end
    is_antisense && complement!(ot)
    spec = prefix_scan_matcher_spec(matcher)
    pos_offset = is_antisense ? spec.rev_pos_offset : spec.fwd_pos_offset
    return ot, candidate_start + pos_offset
end

function evaluate_prefix_hash_scan_candidate!(
    output::Vector{PrefixHashScanVerifiedHit},
    raw::AbstractVector{UInt8},
    geometry::PrefixScanGeometry,
    candidate_start::Int,
    candidate_mask::UInt64,
    global_offset::Int,
    dbi::DBInfo,
    is_antisense::Bool,
    guides_::Vector{LongDNA{4}},
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    stats::S,
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    chunk_counts::Union{Nothing, Matrix{Int}} = nothing,
    ) where {S <: Union{Nothing, PrefixHashScanStats}}

    candidate_mask &= prefix_hash_scan_active_mask(early_stop_state)
    candidate_mask == 0 && return output
    if stats !== nothing
        stats.prefix_hits += 1
        stats.guide_pairs += count_ones(candidate_mask)
    end
    verify_start = prefix_hash_scan_timer(stats)
    ot = LongDNA{4}()
    local_pos = 0
    mask = candidate_mask
    while mask != 0
        guide_idx = trailing_zeros(mask) + 1
        mask &= mask - 1
        align_start = prefix_hash_scan_timer(stats)
        if stats !== nothing
            stats.alignment_calls += 1
            stats.distance_calls += 1
        end
        dist = prefix_hash_scan_raw_myers_distance(
            geometry, myers_profiles[guide_idx], raw, candidate_start,
            is_antisense, distance)
        if dist > distance
            if stats !== nothing
                stats.align_ns += time_ns() - align_start
            end
            continue
        end
        reserve_prefix_hash_scan_hit!(
            early_stop_state, chunk_counts, guide_idx, dist) || continue

        if isempty(ot)
            materialize_start = prefix_hash_scan_timer(stats)
            ot, local_pos = materialize_normalized_candidate_specialized(
                geometry, raw, candidate_start, dbi, is_antisense)
            if stats !== nothing
                stats.candidate_materialize_ns +=
                    time_ns() - materialize_start
            end
        end
        if stats !== nothing
            stats.traceback_calls += 1
        end
        aln = align(guides_[guide_idx], ot, distance, iscompatible)
        if stats !== nothing
            stats.align_ns += time_ns() - align_start
        end
        aln.dist > distance && continue

        if dbi.motif.extends5
            aln_guide = reverse(aln.guide)
            aln_ref = reverse(aln.ref)
        else
            aln_guide = aln.guide
            aln_ref = aln.ref
        end
        push!(output, PrefixHashScanVerifiedHit(
            guide_idx,
            local_pos + global_offset,
            aln.dist,
            is_antisense,
            aln_guide,
            aln_ref,
        ))
    end
    if stats !== nothing
        stats.verify_ns += time_ns() - verify_start
    end
    return output
end

function evaluate_prefix_hash_scan_hits!(
    output::Vector{PrefixHashScanVerifiedHit},
    raw::AbstractVector{UInt8},
    geometry::PrefixScanGeometry,
    hits::Vector{PrefixHashScanHit},
    global_offset::Int,
    dbi::DBInfo,
    is_antisense::Bool,
    guides_::Vector{LongDNA{4}},
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    stats::S,
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    chunk_counts::Union{Nothing, Matrix{Int}} = nothing,
    ) where {S <: Union{Nothing, PrefixHashScanStats}}

    for hit in hits
        evaluate_prefix_hash_scan_candidate!(
            output, raw, geometry, hit.start, hit.mask, global_offset, dbi,
            is_antisense, guides_, myers_profiles, distance, stats,
            early_stop_state, chunk_counts)
    end
    return output
end

function evaluate_prefix_hash_scan_count_candidate!(
    counts::Matrix{Int},
    raw::AbstractVector{UInt8},
    geometry::PrefixScanGeometry,
    candidate_start::Int,
    candidate_mask::UInt64,
    is_antisense::Bool,
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    stats::S,
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    chunk_counts::Union{Nothing, Matrix{Int}} = nothing,
    ) where {S <: Union{Nothing, PrefixHashScanStats}}

    candidate_mask &= prefix_hash_scan_active_mask(early_stop_state)
    candidate_mask == 0 && return counts
    if stats !== nothing
        stats.prefix_hits += 1
        stats.guide_pairs += count_ones(candidate_mask)
    end
    verify_start = prefix_hash_scan_timer(stats)
    mask = candidate_mask
    while mask != 0
        guide_idx = trailing_zeros(mask) + 1
        mask &= mask - 1
        align_start = prefix_hash_scan_timer(stats)
        if stats !== nothing
            stats.alignment_calls += 1
            stats.distance_calls += 1
        end
        dist = prefix_hash_scan_raw_myers_distance(
            geometry, myers_profiles[guide_idx], raw, candidate_start,
            is_antisense, distance)
        if stats !== nothing
            stats.align_ns += time_ns() - align_start
        end
        if dist <= distance
            if early_stop_state === nothing
                counts[guide_idx, dist + 1] += 1
            else
                reserve_prefix_hash_scan_hit!(
                    early_stop_state, chunk_counts, guide_idx, dist)
            end
        end
    end
    if stats !== nothing
        stats.verify_ns += time_ns() - verify_start
    end
    return counts
end

function evaluate_prefix_hash_scan_count_hits!(
    counts::Matrix{Int},
    raw::AbstractVector{UInt8},
    geometry::PrefixScanGeometry,
    hits::Vector{PrefixHashScanHit},
    is_antisense::Bool,
    myers_profiles::Vector{PrefixHashScanMyersProfile},
    distance::Int,
    stats::S,
    early_stop_state::Union{Nothing, PrefixHashScanEarlyStopState} = nothing,
    chunk_counts::Union{Nothing, Matrix{Int}} = nothing,
    ) where {S <: Union{Nothing, PrefixHashScanStats}}

    for hit in hits
        evaluate_prefix_hash_scan_count_candidate!(
            counts, raw, geometry, hit.start, hit.mask, is_antisense,
            myers_profiles, distance, stats, early_stop_state, chunk_counts)
    end
    return counts
end

function merge_prefix_hash_scan_worker_stats!(
    stats::PrefixHashScanStats,
    worker_stats::PrefixHashScanStats)

    stats.motif_candidates += worker_stats.motif_candidates
    stats.ambiguous_prefixes += worker_stats.ambiguous_prefixes
    stats.prefix_hits += worker_stats.prefix_hits
    stats.guide_pairs += worker_stats.guide_pairs
    stats.alignment_calls += worker_stats.alignment_calls
    stats.distance_calls += worker_stats.distance_calls
    stats.traceback_calls += worker_stats.traceback_calls
    stats.record_io_ns += worker_stats.record_io_ns
    stats.chrom_load_ns += worker_stats.chrom_load_ns
    stats.candidate_materialize_ns += worker_stats.candidate_materialize_ns
    stats.align_ns += worker_stats.align_ns
    stats.verify_ns += worker_stats.verify_ns
    return stats
end

function commit_prefix_hash_scan_verified!(
    out,
    hit::PrefixHashScanVerifiedHit,
    guide::LongDNA{4},
    chrom_name::String,
    early_stopping::Vector{Int},
    es_acc,
    is_es,
    seen,
    stats::Union{Nothing, PrefixHashScanStats};
    prelimited::Bool = false)

    guide_idx = hit.guide_idx
    !prelimited && is_es[guide_idx] && return nothing
    strand = hit.is_antisense ? "-" : "+"
    key = (
        string(guide),
        hit.dist,
        chrom_name,
        hit.pos,
        strand,
        hit.aln_guide,
        hit.aln_ref,
    )
    key in seen[guide_idx] && return nothing
    push!(seen[guide_idx], key)
    dist_idx = hit.dist + 1
    if prelimited
        es_acc[guide_idx, dist_idx] >= early_stopping[dist_idx] && return nothing
        es_acc[guide_idx, dist_idx] += 1
    else
        reserve_prefix_hash_scan_detail_hit!(
            es_acc, is_es, guide_idx, hit.dist, early_stopping) || return nothing
    end

    emit_start = prefix_hash_scan_timer(stats)
    print(out, guide, ",", hit.aln_guide, ",", hit.aln_ref, ",",
        hit.dist, ",", chrom_name, ",", hit.pos, ",", strand, "\n")
    if stats !== nothing
        stats.emit_ns += time_ns() - emit_start
        stats.emitted_rows += 1
    end
    return nothing
end

# Commits streamed hits in deterministic output order: chromosome, then strand
# (plus before minus), then chunk.
function commit_prefix_hash_scan_chunks!(
    out,
    chunk_results,
    chrom_chunk_ranges,
    chrom_names::Vector{String},
    guides::Vector{LongDNA{4}},
    early_stopping::Vector{Int},
    es_acc,
    is_es,
    seen,
    stats::Union{Nothing, PrefixHashScanStats};
    prelimited::Bool)

    for chrom_idx in eachindex(chrom_chunk_ranges)
        chrom_name = chrom_names[chrom_idx]
        chunk_range = chrom_chunk_ranges[chrom_idx]
        if stats !== nothing
            for chunk_idx in chunk_range
                result_ = chunk_results[chunk_idx]
                result_ === nothing && continue
                merge_prefix_hash_scan_worker_stats!(stats, result_.stats)
            end
        end
        for strand in (:plus, :minus)
            for chunk_idx in chunk_range
                result_ = chunk_results[chunk_idx]
                result_ === nothing && continue
                hits = getfield(result_, strand)
                for hit in hits
                    commit_prefix_hash_scan_verified!(
                        out,
                        hit,
                        guides[hit.guide_idx],
                        chrom_name,
                        early_stopping,
                        es_acc,
                        is_es,
                        seen,
                        stats,
                        prelimited = prelimited,
                    )
                end
            end
        end
    end
    return nothing
end
