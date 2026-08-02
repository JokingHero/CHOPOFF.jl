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

@inline function prefix_hash_scan_raw_myers_distance(
    ::PrefixScanGeometry{:cas9},
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
        raw_idx = is_antisense ?
            candidate_start + 2 + ref_idx :
            candidate_start + pattern_length - ref_idx
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

@inline function prefix_hash_scan_raw_myers_distance(
    ::PrefixScanGeometry{:cas12a},
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
        raw_idx = is_antisense ?
            candidate_start + pattern_length - ref_idx :
            candidate_start + 3 + ref_idx
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

@inline prefix_hash_scan_raw_myers_distance(args...) =
    prefix_hash_scan_raw_myers_distance(
        CAS9_D3_PREFIX_SCAN_GEOMETRY, args...)

@inline function prefix_hash_scan_raw_myers_distance(
    geometry::PrefixScanGeometry{:generic},
    profile::PrefixHashScanMyersProfile,
    raw::AbstractVector{UInt8},
    candidate_start::Int,
    is_antisense::Bool,
    distance::Int)

    spec = prefix_scan_matcher_spec(geometry.matcher)
    offsets = is_antisense ? spec.rev_offsets : spec.fwd_offsets
    step = is_antisense ? spec.rev_step : spec.fwd_step
    pattern_length = Int(profile.length)
    reference_length = pattern_length + distance
    first_scored_prefix = pattern_length - distance
    pv = typemax(UInt64)
    mv = UInt64(0)
    score = pattern_length
    best = distance + 1
    @inbounds for ref_idx in 1:reference_length
        relative_idx = ref_idx <= pattern_length ? offsets[ref_idx] :
            (step > 0 ? spec.span + ref_idx - pattern_length - 1 :
                -(ref_idx - pattern_length))
        raw_idx = candidate_start + relative_idx
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

function materialize_normalized_candidate_cas9(
    chrom_seq::LongDNA{4},
    candidate_start::Int,
    dbi::DBInfo,
    is_antisense::Bool)

    distance = dbi.motif.distance
    geometry = CAS9_D3_PREFIX_SCAN_GEOMETRY
    candidate_last_offset = prefix_scan_candidate_last_offset(geometry)
    guide_last_offset = geometry.guide_bases - 1
    if is_antisense
        guide_start = candidate_start + geometry.pam_bases
        extension_end = candidate_start + candidate_last_offset + distance
        if extension_end <= length(chrom_seq)
            ot = complement(chrom_seq[guide_start:extension_end])
        else
            guide = chrom_seq[guide_start:(candidate_start + candidate_last_offset)]
            extension = getExt3(
                chrom_seq, length(chrom_seq),
                candidate_start + prefix_scan_candidate_bases(geometry), distance)
            ot = complement(guide * extension)
        end
        return ot, candidate_start
    end

    extension_start = candidate_start - distance
    if extension_start >= 1
        ot = reverse(chrom_seq[extension_start:(candidate_start + guide_last_offset)])
    else
        extension = getExt5(chrom_seq, candidate_start - 1, distance)
        guide = chrom_seq[candidate_start:(candidate_start + guide_last_offset)]
        ot = reverse(extension * guide)
    end
    return ot, candidate_start + candidate_last_offset
end

function materialize_normalized_candidate_cas12a(
    chrom_seq::LongDNA{4},
    candidate_start::Int,
    dbi::DBInfo,
    is_antisense::Bool)

    distance = dbi.motif.distance
    if is_antisense
        extension_start = candidate_start - distance
        if extension_start >= 1
            ot = reverse_complement(
                chrom_seq[extension_start:(candidate_start + 20)])
        else
            extension = getExt5(chrom_seq, candidate_start - 1, distance)
            guide = chrom_seq[candidate_start:(candidate_start + 20)]
            ot = reverse_complement(extension * guide)
        end
        return ot, candidate_start + 24
    end

    extension_end = candidate_start + 24 + distance
    if extension_end <= length(chrom_seq)
        ot = chrom_seq[(candidate_start + 4):extension_end]
    else
        guide = chrom_seq[(candidate_start + 4):(candidate_start + 24)]
        extension = getExt3(
            chrom_seq, length(chrom_seq), candidate_start + 25, distance)
        ot = guide * extension
    end
    return ot, candidate_start
end


function materialize_normalized_candidate_cas9(
    raw::AbstractVector{UInt8},
    candidate_start::Int,
    dbi::DBInfo,
    is_antisense::Bool)

    distance = dbi.motif.distance
    geometry = CAS9_D3_PREFIX_SCAN_GEOMETRY
    candidate_last_offset = prefix_scan_candidate_last_offset(geometry)
    guide_last_offset = geometry.guide_bases - 1
    n = length(raw)
    if is_antisense
        guide_start = candidate_start + geometry.pam_bases
        extension_end = candidate_start + candidate_last_offset + distance
        if extension_end <= n
            ot = complement(LongDNA{4}(@view raw[guide_start:extension_end]))
        else
            guide = LongDNA{4}(
                @view raw[guide_start:(candidate_start + candidate_last_offset)])
            available_start = candidate_start + prefix_scan_candidate_bases(geometry)
            available = available_start <= n ?
                LongDNA{4}(@view raw[available_start:n]) :
                LongDNA{4}()
            extension = available * LongDNA{4}(repeat("-", distance - length(available)))
            ot = complement(guide * extension)
        end
        return ot, candidate_start
    end

    extension_start = candidate_start - distance
    if extension_start >= 1
        ot = reverse(LongDNA{4}(
            @view raw[extension_start:(candidate_start + guide_last_offset)]))
    else
        available_end = candidate_start - 1
        available = available_end >= 1 ?
            LongDNA{4}(@view raw[1:available_end]) :
            LongDNA{4}()
        extension = LongDNA{4}(repeat("-", distance - length(available))) * available
        guide = LongDNA{4}(
            @view raw[candidate_start:(candidate_start + guide_last_offset)])
        ot = reverse(extension * guide)
    end
    return ot, candidate_start + candidate_last_offset
end

function materialize_normalized_candidate_cas12a(
    raw::AbstractVector{UInt8},
    candidate_start::Int,
    dbi::DBInfo,
    is_antisense::Bool)

    distance = dbi.motif.distance
    n = length(raw)
    if is_antisense
        extension_start = candidate_start - distance
        if extension_start >= 1
            ot = reverse_complement(LongDNA{4}(
                @view raw[extension_start:(candidate_start + 20)]))
        else
            available_end = candidate_start - 1
            available = available_end >= 1 ?
                LongDNA{4}(@view raw[1:available_end]) :
                LongDNA{4}()
            extension =
                LongDNA{4}(repeat("-", distance - length(available))) * available
            guide = LongDNA{4}(
                @view raw[candidate_start:(candidate_start + 20)])
            ot = reverse_complement(extension * guide)
        end
        return ot, candidate_start + 24
    end

    extension_end = candidate_start + 24 + distance
    if extension_end <= n
        ot = LongDNA{4}(@view raw[(candidate_start + 4):extension_end])
    else
        guide = LongDNA{4}(
            @view raw[(candidate_start + 4):(candidate_start + 24)])
        available_start = candidate_start + 25
        available = available_start <= n ?
            LongDNA{4}(@view raw[available_start:n]) :
            LongDNA{4}()
        extension = available * LongDNA{4}(repeat("-", distance - length(available)))
        ot = guide * extension
    end
    return ot, candidate_start
end

materialize_normalized_candidate_specialized(
    ::PrefixScanGeometry{:cas9}, args...) =
    materialize_normalized_candidate_cas9(args...)

materialize_normalized_candidate_specialized(
    ::PrefixScanGeometry{:cas12a}, args...) =
    materialize_normalized_candidate_cas12a(args...)

function materialize_normalized_candidate_specialized(
    geometry::PrefixScanGeometry{:generic},
    chrom_seq::LongDNA{4}, candidate_start::Int,
    dbi::DBInfo, is_antisense::Bool)

    candidate_range = candidate_start:(
        candidate_start + prefix_scan_candidate_last_offset(geometry))
    return materialize_normalized_candidate(
        chrom_seq, candidate_range, dbi, is_antisense)
end

function materialize_normalized_candidate_specialized(
    geometry::PrefixScanGeometry{:generic},
    raw::AbstractVector{UInt8}, candidate_start::Int,
    dbi::DBInfo, is_antisense::Bool)

    spec = prefix_scan_matcher_spec(geometry.matcher)
    offsets = is_antisense ? spec.rev_offsets : spec.fwd_offsets
    step = is_antisense ? spec.rev_step : spec.fwd_step
    bytes = Vector{UInt8}(undef, geometry.guide_bases + geometry.distance)
    @inbounds for ref_idx in eachindex(bytes)
        relative_idx = ref_idx <= geometry.guide_bases ? offsets[ref_idx] :
            (step > 0 ? spec.span + ref_idx - geometry.guide_bases - 1 :
                -(ref_idx - geometry.guide_bases))
        raw_idx = candidate_start + relative_idx
        bytes[ref_idx] = 1 <= raw_idx <= length(raw) ?
            raw[raw_idx] : UInt8('-')
    end
    ot = LongDNA{4}(String(bytes))
    is_antisense && (ot = complement(ot))
    pos_offset = is_antisense ? spec.rev_pos_offset : spec.fwd_pos_offset
    return ot, candidate_start + pos_offset
end

function verify_prefix_hash_scan_bitmask_candidate!(
    out,
    chrom_seq,
    candidate_range::UnitRange{Int},
    geometry::PrefixScanGeometry,
    dbi::DBInfo,
    is_antisense::Bool,
    candidate_mask::UInt64,
    guides::Vector{LongDNA{4}},
    guides_::Vector{LongDNA{4}},
    chrom_name::String,
    distance::Int,
    early_stopping::Vector{Int},
    es_acc,
    is_es,
    seen,
    stats::Union{Nothing, PrefixHashScanStats};
    distance_first::Bool = false,
    myers_profiles::Union{Nothing, Vector{PrefixHashScanMyersProfile}} = nothing)

    strand = is_antisense ? "-" : "+"
    ot = LongDNA{4}()
    pos = 0
    verify_start = time_ns()
    mask = candidate_mask
    while mask != 0
        guide_idx = trailing_zeros(mask) + 1
        mask &= mask - 1
        is_es[guide_idx] && continue

        if stats !== nothing
            stats.alignment_calls += 1
        end
        align_start = time_ns()
        if myers_profiles !== nothing
            if stats !== nothing
                stats.distance_calls += 1
            end
            dist = prefix_hash_scan_raw_myers_distance(
                geometry, myers_profiles[guide_idx], chrom_seq, first(candidate_range),
                is_antisense, distance)
            if dist > distance
                if stats !== nothing
                    stats.align_ns += time_ns() - align_start
                end
                continue
            end
        elseif distance_first
            if isempty(ot)
                materialize_start = time_ns()
                ot, pos = materialize_normalized_candidate_specialized(
                    geometry, chrom_seq, first(candidate_range), dbi, is_antisense)
                if stats !== nothing
                    stats.candidate_materialize_ns += time_ns() - materialize_start
                end
            end
            if stats !== nothing
                stats.distance_calls += 1
            end
            dist = levenshtein(guides_[guide_idx], ot, distance, iscompatible)
            if dist > distance
                if stats !== nothing
                    stats.align_ns += time_ns() - align_start
                end
                continue
            end
        end
        if isempty(ot)
            materialize_start = time_ns()
            ot, pos = materialize_normalized_candidate_specialized(
                geometry, chrom_seq, first(candidate_range), dbi, is_antisense)
            if stats !== nothing
                stats.candidate_materialize_ns += time_ns() - materialize_start
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

        key = (
            string(guides[guide_idx]),
            aln.dist,
            chrom_name,
            Int(pos),
            strand,
            aln_guide,
            aln_ref,
        )
        key in seen[guide_idx] && continue
        push!(seen[guide_idx], key)

        emit_start = time_ns()
        print(out, guides[guide_idx], ",", aln_guide, ",", aln_ref, ",",
            aln.dist, ",", chrom_name, ",", pos, ",", strand, "\n")
        if stats !== nothing
            stats.emit_ns += time_ns() - emit_start
            stats.emitted_rows += 1
        end
        es_acc[guide_idx, aln.dist + 1] += 1
        if es_acc[guide_idx, aln.dist + 1] >= early_stopping[aln.dist + 1]
            is_es[guide_idx] = true
        end
    end
    if stats !== nothing
        stats.verify_ns += time_ns() - verify_start
    end
    return nothing
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
    stats::S) where {S <: Union{Nothing, PrefixHashScanStats}}

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
    stats::S) where {S <: Union{Nothing, PrefixHashScanStats}}

    for hit in hits
        evaluate_prefix_hash_scan_candidate!(
            output, raw, geometry, hit.start, hit.mask, global_offset, dbi,
            is_antisense, guides_, myers_profiles, distance, stats)
    end
    return output
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
    stats::Union{Nothing, PrefixHashScanStats})

    guide_idx = hit.guide_idx
    is_es[guide_idx] && return nothing
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

    emit_start = prefix_hash_scan_timer(stats)
    print(out, guide, ",", hit.aln_guide, ",", hit.aln_ref, ",",
        hit.dist, ",", chrom_name, ",", hit.pos, ",", strand, "\n")
    if stats !== nothing
        stats.emit_ns += time_ns() - emit_start
        stats.emitted_rows += 1
    end
    es_acc[guide_idx, hit.dist + 1] += 1
    if es_acc[guide_idx, hit.dist + 1] >= early_stopping[hit.dist + 1]
        is_es[guide_idx] = true
    end
    return nothing
end
