struct ScanVariantWork
    fragments::UnitRange{Int}
    core_first::Int            # zero for a packed group of complete fragments
    core_last::Int
end

const ScanVariantRow = Tuple{String, String, String, Int, String, Int, String, Int, String, String, String, String}

function prefix_hash_scan_variant_work(data, geometry, chunk_bases)
    work = ScanVariantWork[]
    first = 1
    bases = 0
    separator = geometry.distance + 1
    for (i, fragment) in enumerate(data.fragments)
        if bases > 0 && bases + fragment.bases + separator > chunk_bases
            push!(work, ScanVariantWork(first:(i - 1), 0, 0))
            first, bases = i, 0
        end
        if fragment.bases > chunk_bases
            last = fragment.bases - prefix_scan_candidate_last_offset(geometry)
            for start in 1:chunk_bases:last
                push!(work, ScanVariantWork(i:i, start, min(last, start + chunk_bases - 1)))
            end
            first, bases = i + 1, 0
        else
            bases += fragment.bases + separator
        end
    end
    bases > 0 && push!(work, ScanVariantWork(first:length(data.fragments), 0, 0))
    return work
end

function reserve_variant_hit!(state, seen, row, guide_idx)
    lock(state.merge_lock) do
        row in seen && return false
        bit = UInt64(1) << (guide_idx - 1)
        state.active_mask[] & bit == 0 && return false
        bucket = row[4] + 1
        if state.counts[guide_idx, bucket] >= state.limits[bucket]
            retire_prefix_hash_scan_guide!(state, guide_idx)
            return false
        end
        push!(seen, row)
        state.counts[guide_idx, bucket] += 1
        return true
    end
end

function scan_prefix_hash_scan_variant_work(ctx, worker, scratch, data, work, guides, state, seen, stats)
    raw = worker.buffer
    empty!(raw)
    starts, stops, shifts = Int[], Int[], Int[]
    io_start = prefix_hash_scan_timer(stats)
    for fi in work.fragments
        fragment = data.fragments[fi]
        first = work.core_first == 0 ? 1 : max(1, work.core_first - ctx.distance)
        last = work.core_first == 0 ? fragment.bases : min(fragment.bases,
            work.core_last + prefix_scan_candidate_last_offset(ctx.geometry) + ctx.distance)
        start = length(raw) + 1
        push!(starts, start)
        push!(shifts, first - start)
        materialize_variant_fragment!(raw, scratch, worker.io, ctx.index, data, fragment, first, last)
        push!(stops, length(raw))
        fi == lastindex(data.fragments) || append!(raw, fill(UInt8('-'), ctx.distance + 1))
    end
    stats === nothing || (stats.variant_io_ns += time_ns() - io_start)
    bounds = prefix_scan_bounds(ctx.geometry, raw, ctx.dbi)
    bounds === nothing && return ScanVariantRow[]
    if work.core_first != 0
        core = (work.core_first - shifts[1]):(work.core_last - shifts[1])
        bounds = prefix_scan_bounds(intersect(bounds.plus, core), intersect(bounds.minus, core))
        bounds === nothing && return ScanVariantRow[]
    end
    lookup = worker.lookup
    scan_generic_prefix_hits_raw_range_bucketed!(
        worker.plus_hits, worker.minus_hits, lookup.plus_candidates, lookup.minus_candidates,
        lookup.plus_radix, lookup.minus_radix, lookup.radix_counts,
        raw, ctx.query, ctx.geometry, bounds, ctx.simd_backend)
    if ctx.dbi.motif.ambig_max > 0
        scan_ambiguous_prefix_hits_range!(worker.plus_hits, worker.minus_hits, raw,
            ctx.geometry, ctx.dbi, ctx.query, ctx.geometry.prefix_bases, bounds,
            Val(ctx.dbi.motif.ambig_max), stats, ctx.simd_backend)
    end
    rows = ScanVariantRow[]
    spec = prefix_scan_matcher_spec(ctx.geometry.matcher)
    for (antisense, hits) in ((false, worker.plus_hits), (true, worker.minus_hits)), hit in hits
        mask = hit.mask & prefix_hash_scan_active_mask(state)
        mask == 0 && continue
        stats === nothing || (stats.variant_candidates += 1)
        slot = searchsortedlast(starts, hit.start)
        slot == 0 && continue
        hit.start + spec.span - 1 <= stops[slot] || continue
        fi = first(work.fragments) + slot - 1
        fragment = data.fragments[fi]
        shift = shifts[slot]
        extended = (hit.start + prefix_scan_reference_offset(ctx.geometry.matcher, antisense, 1),
            hit.start + prefix_scan_reference_offset(ctx.geometry.matcher, antisense, ctx.geometry.guide_bases + ctx.distance))
        lo = min(hit.start, extended...) + shift
        hi = max(hit.start + spec.span - 1, extended...) + shift
        if !variant_anchor_overlaps(data, fragment, lo, hi)
            stats === nothing || (stats.anchor_rejected += 1)
            continue
        end
        while mask != 0
            gi = trailing_zeros(mask) + 1
            mask &= mask - 1
            dist = prefix_hash_scan_raw_myers_distance(ctx.geometry, ctx.myers_profiles[gi],
                raw, hit.start, antisense, ctx.distance)
            dist <= ctx.distance || continue
            ot, local_pos = materialize_normalized_candidate_specialized(
                ctx.geometry, raw, hit.start, ctx.dbi, antisense)
            aln = align(ctx.guides_[gi], ot, ctx.distance, iscompatible)
            aln.dist <= ctx.distance || continue
            lo, hi = prefix_hash_scan_site_span(ctx.geometry, ctx.dbi.motif, hit.start, antisense, aln.ref)
            lo += shift
            hi += shift
            alleles = variant_site_alleles(data, fragment, lo, hi)
            # Candidate context can change the best traceback even when that
            # traceback avoids the anchor. Keep its site alleles and deduplicate
            # the resulting record, rather than rejecting another valid alignment.
            isempty(alleles) && continue
            start, start_origin = lift_variant_position(data, fragment, local_pos + shift)
            other = (antisense ? spec.rev_pos_offset : spec.fwd_pos_offset) == 0 ? hi : lo
            stop, end_origin = lift_variant_position(data, fragment, other)
            aln_guide = ctx.dbi.motif.extends5 ? reverse(aln.guide) : aln.guide
            aln_ref = ctx.dbi.motif.extends5 ? reverse(aln.ref) : aln.ref
            annotations = join((data.variants[i].annotation for i in alleles), ';')
            row = (string(guides[gi]), aln_guide, aln_ref, aln.dist, ctx.dbi.gi.chrom[fragment.chrom],
                start, antisense ? "-" : "+", stop, annotations, "", start_origin, end_origin)
            reserve_variant_hit!(state, seen, row, gi) && push!(rows, row)
        end
    end
    return rows
end

function stream_prefix_hash_scan_variants!(out, data, geometry, query, dbi, guides,
    guides_, myers_profiles, distance, chunk_bases, scan_threads, limits, stats, simd_backend)

    work = prefix_hash_scan_variant_work(data, geometry, chunk_bases)
    isempty(work) && return nothing
    index = dbi.gi.is_fa ? FASTA.Index(dbi.gi.filepath * ".fai") : read_prefix_hash_scan_twobit_index(dbi.gi.filepath)
    ctx = PrefixHashScanContext(geometry, query, dbi, guides_, myers_profiles, distance, index, simd_backend)
    state = PrefixHashScanEarlyStopState(length(guides), limits)
    seen = Set{ScanVariantRow}()
    results = Vector{Vector{ScanVariantRow}}(undef, length(work))
    next = Threads.Atomic{Int}(1)
    tasks = map(1:min(scan_threads, length(work))) do _
        Threads.@spawn begin
            worker = PrefixHashScanWorker(open(dbi.gi.filepath, "r"))
            scratch = UInt8[]
            worker_stats = prefix_hash_scan_worker_stats(stats)
            try
                while true
                    i = Threads.atomic_add!(next, 1)
                    i > length(work) && break
                    if prefix_hash_scan_active_mask(state) == 0
                        results[i] = ScanVariantRow[]
                    else
                        results[i] = scan_prefix_hash_scan_variant_work(ctx, worker, scratch,
                            data, work[i], guides, state, seen, worker_stats)
                    end
                end
            finally
                close(worker.io)
            end
            worker_stats
        end
    end
    for task in tasks
        ws = fetch(task)
        if stats !== nothing
            stats.variant_candidates += ws.variant_candidates
            stats.anchor_rejected += ws.anchor_rejected
            stats.variant_io_ns += ws.variant_io_ns
        end
    end
    # Sorting accepted rows makes uncapped output independent of which worker
    # encountered an identical fragment first.
    rows = reduce(vcat, results; init = ScanVariantRow[])
    sort!(rows)
    for row in rows
        write_variant_csv_row(out, row)
    end
    if stats !== nothing
        stats.emitted_rows += length(rows)
        stats.variant_rows += length(rows)
        stats.variant_retired_guides += state.retired_guides[]
    end
    state.retired_guides[] > 0 && @warn "VCF search reached variant early-stopping limits; output is incomplete" guides=state.retired_guides[]
    return nothing
end
