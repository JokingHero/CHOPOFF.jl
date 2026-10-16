struct ScanVariantWork
    fragments::UnitRange{Int}
    core_first::Int            # zero for a packed group of complete fragments
    core_last::Int
end

const ScanVariantRow = Tuple{String, String, String, Int, String, Int, String, Int, String, String, String, String}

# Variant work items per scan thread; more than one balances uneven items.
const PREFIX_HASH_SCAN_VARIANT_ITEMS_PER_THREAD = 4
const PREFIX_HASH_SCAN_VARIANT_MIN_ITEM_BASES = 1 << 16

# Sparse VCFs pack into few chunk-sized items, which would idle most workers.
function prefix_hash_scan_variant_item_bases(blocks, geometry, chunk_bases, scan_threads)
    separator = geometry.distance + 1
    total = sum(b -> sum(f -> f.bases + separator, b.fragments; init = 0), blocks; init = 0)
    return min(chunk_bases, max(PREFIX_HASH_SCAN_VARIANT_MIN_ITEM_BASES,
        cld(total, PREFIX_HASH_SCAN_VARIANT_ITEMS_PER_THREAD * scan_threads)))
end

function prefix_hash_scan_variant_work(chrom, geometry, chunk_bases, target)
    work = ScanVariantWork[]
    first = 1
    bases = 0
    separator = geometry.distance + 1
    for (i, fragment) in enumerate(chrom.fragments)
        if bases > 0 && bases + fragment.bases + separator > target
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
    bases > 0 && push!(work, ScanVariantWork(first:length(chrom.fragments), 0, 0))
    return work
end

const ScanVariantGenotypeColumns = NTuple{6, String}

variant_row_seen(state, seen, row) = lock(() -> row in seen, state.merge_lock)

# Genotype rows carried by no selected copy are remembered but not counted.
function reserve_variant_hit!(state, seen, row, guide_idx, genotype_rows = nothing, columns = nothing)
    lock(state.merge_lock) do
        row in seen && return false
        if columns !== nothing && isempty(columns[1])
            push!(seen, row)
            return false
        end
        bit = UInt64(1) << (guide_idx - 1)
        state.active_mask[] & bit == 0 && return false
        bucket = row[4] + 1
        if state.counts[guide_idx, bucket] >= state.limits[bucket]
            retire_prefix_hash_scan_guide!(state, guide_idx)
            return false
        end
        push!(seen, row)
        columns === nothing || (genotype_rows[row] = columns)
        state.counts[guide_idx, bucket] += 1
        return true
    end
end

function scan_prefix_hash_scan_variant_work(ctx, worker, data, chrom, work, guides, state, seen, stats,
    genotype_rows = nothing)
    raw = worker.buffer
    empty!(raw)
    starts, stops, shifts = Int[], Int[], Int[]
    io_start = prefix_hash_scan_timer(stats)
    for fi in work.fragments
        fragment = chrom.fragments[fi]
        first = work.core_first == 0 ? 1 : max(1, work.core_first - ctx.distance)
        last = work.core_first == 0 ? fragment.bases : min(fragment.bases,
            work.core_last + prefix_scan_candidate_last_offset(ctx.geometry) + ctx.distance)
        start = length(raw) + 1
        push!(starts, start)
        push!(shifts, first - start)
        materialize_variant_fragment!(raw, data, chrom, fragment, first, last)
        push!(stops, length(raw))
        fi == lastindex(chrom.fragments) || append!(raw, fill(UInt8('-'), ctx.distance + 1))
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
        fragment = chrom.fragments[fi]
        shift = shifts[slot]
        extended = (hit.start + prefix_scan_reference_offset(ctx.geometry.matcher, antisense, 1),
            hit.start + prefix_scan_reference_offset(ctx.geometry.matcher, antisense, ctx.geometry.guide_bases + ctx.distance))
        lo = min(hit.start, extended...) + shift
        hi = max(hit.start + spec.span - 1, extended...) + shift
        if !variant_anchor_overlaps(chrom, fragment, lo, hi)
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
            alleles = variant_site_alleles(chrom, fragment, lo, hi)
            # Candidate context can change the best traceback even when that
            # traceback avoids the anchor. Keep its site alleles and deduplicate
            # the resulting record, rather than rejecting another valid alignment.
            isempty(alleles) && continue
            start, start_origin = lift_variant_position(data, chrom, fragment, local_pos + shift)
            other = (antisense ? spec.rev_pos_offset : spec.fwd_pos_offset) == 0 ? hi : lo
            stop, end_origin = lift_variant_position(data, chrom, fragment, other)
            data.intervals === nothing || prefix_hash_scan_overlaps(data.intervals[fragment.chrom],
                min(start, stop), max(start, stop)) || continue
            aln_guide = ctx.dbi.motif.extends5 ? reverse(aln.guide) : aln.guide
            aln_ref = ctx.dbi.motif.extends5 ? reverse(aln.ref) : aln.ref
            annotations = join((data.variants[i].annotation for i in alleles), ';')
            row = (string(guides[gi]), aln_guide, aln_ref, aln.dist, ctx.dbi.gi.chrom[fragment.chrom],
                start, antisense ? "-" : "+", stop, annotations, "", start_origin, end_origin)
            columns = nothing
            if data.genotypes !== nothing
                # Attribution depends only on the row, so duplicates skip it.
                variant_row_seen(state, seen, row) && continue
                taken, excluded = variant_site_requirements(data, chrom, fragment, lo, hi)
                carriers, certain, possible = haplotype_tokens(data.genotypes, taken, excluded)
                columns = (carriers, "", string(certain), string(possible), "", "")
            end
            reserve_variant_hit!(state, seen, row, gi, genotype_rows, columns) && push!(rows, row)
        end
    end
    return rows
end

function stream_prefix_hash_scan_variants!(out, data, geometry, query, dbi, guides,
    guides_, myers_profiles, distance, chunk_bases, scan_threads, limits, stats, simd_backend;
    progress = nothing, report = true)

    lengths = dbi.gi.is_fa ? FASTA.Index(dbi.gi.filepath * ".fai").lengths :
        read_prefix_hash_scan_twobit_index(dbi.gi.filepath).lengths
    chroms = [ci for ci in eachindex(data.anchors) if !isempty(data.anchors[ci])]
    interval = progress === nothing ? 0.0 : progress.interval_ns / 1e9
    enumerate_progress = progress === nothing ? nothing :
        PrefixHashScanProgress(progress.label * " VCF enumerate", interval)
    enumerate_progress === nothing || (enumerate_progress.total_bases = sum(ci -> length(data.anchors[ci]), chroms; init = 0))
    scan_progress = progress === nothing ? nothing : PrefixHashScanProgress(progress.label * " VCF", interval)
    index = dbi.gi.is_fa ? FASTA.Index(dbi.gi.filepath * ".fai") : read_prefix_hash_scan_twobit_index(dbi.gi.filepath)
    ctx = PrefixHashScanContext(geometry, query, dbi, guides_, myers_profiles, distance, index, simd_backend)
    state = PrefixHashScanEarlyStopState(length(guides), limits)
    seen = Set{ScanVariantRow}()
    genotype_rows = data.genotypes === nothing ? nothing : Dict{ScanVariantRow, ScanVariantGenotypeColumns}()
    rows = ScanVariantRow[]
    fragments, fragment_bases, max_anchor, chrom_bytes, enumeration_ns = 0, 0, 0, 0, UInt64(0)
    truncated = Tuple{Int, Int}[]
    dense = Pair{String, Int}[]
    dense_anchors = 0
    for ci in chroms
        prefix_hash_scan_active_mask(state) == 0 && break
        started = time_ns()
        blocks = enumerate_variant_chromosome(data, dbi, lengths, ci, scan_threads; progress = enumerate_progress)
        enumeration_ns += time_ns() - started
        fragments += sum(b -> length(b.fragments), blocks)
        fragment_bases += sum(b -> sum(f -> f.bases, b.fragments; init = 0), blocks)
        max_anchor = max(max_anchor, maximum(b -> maximum(b.per_anchor; init = 0), blocks))
        chrom_bytes = max(chrom_bytes, sum(variant_descriptor_bytes, blocks))
        per_anchor = Iterators.flatten(b.per_anchor for b in blocks)
        append!(truncated, ((ci, data.variants[i].pos) for b in blocks for i in b.truncated))
        for (i, n) in zip(data.anchors[ci], per_anchor)
            n > PREFIX_HASH_SCAN_VARIANT_FRAGMENT_WARN || continue
            dense_anchors += 1
            push!(dense, string(dbi.gi.chrom[ci], ":", data.variants[i].pos) => n)
            sort!(dense; by = x -> -x[2])
            length(dense) > 5 && pop!(dense)
        end
        append!(rows, scan_prefix_hash_scan_variant_chromosome(ctx, data, blocks, guides, state, seen,
            stats, genotype_rows, geometry, chunk_bases, scan_threads, scan_progress))
        # Without a collection here, eight enumerating threads left about twice
        # the live fragment memory as garbage (genome-wide GIAB: 8.6 vs 4.0 GB
        # peak RSS); a full collection cost 18 s, an incremental one nothing.
        blocks = nothing
        GC.gc(false)
    end
    if report
        if !isempty(truncated)
            @warn "VCF allele combinations truncated at vcf_max_alleles=$(data.max_alleles); output omits local haplotypes with more ALT alleles" anchors=length(truncated) examples=unique([string(dbi.gi.chrom[c], ":", p) for (c, p) in first(truncated, 5)])
        end
        dense_anchors > 0 && @warn "VCF allele combinations are dense; preparation time and memory grow with fragments" anchors=dense_anchors fragments=fragments worst=dense
        progress === nothing || @info "prefixHashScan variant fragments" fragments max_anchor_fragments=max_anchor truncated_anchors=length(truncated) fragment_bases largest_chromosome_bytes=chrom_bytes enumeration_s=enumeration_ns / 1e9
    end
    # Sorting accepted rows makes uncapped output independent of which worker
    # encountered an identical fragment first.
    sort!(rows)
    for row in rows
        write_variant_csv_row(out, genotype_rows === nothing ? row : (row..., genotype_rows[row]...))
    end
    if stats !== nothing
        stats.fragments += fragments
        stats.fragment_bases += fragment_bases
        stats.max_anchor_fragments = max(stats.max_anchor_fragments, max_anchor)
        stats.vcf_truncated_anchors += length(truncated)
        stats.descriptor_bytes += chrom_bytes
        stats.enumeration_ns += enumeration_ns
        stats.emitted_rows += length(rows)
        stats.variant_rows += length(rows)
        stats.variant_retired_guides += state.retired_guides[]
    end
    state.retired_guides[] > 0 && @warn "VCF search reached variant early-stopping limits; output is incomplete" guides=state.retired_guides[]
    return nothing
end

function scan_prefix_hash_scan_variant_chromosome(ctx, data, blocks, guides, state, seen, stats,
    genotype_rows, geometry, chunk_bases, scan_threads, progress)

    target = prefix_hash_scan_variant_item_bases(blocks, geometry, chunk_bases, scan_threads)
    work = [(block, w) for block in blocks for w in prefix_hash_scan_variant_work(block, geometry, chunk_bases, target)]
    isempty(work) && return ScanVariantRow[]
    work_bases((block, w)) = w.core_first == 0 ? sum(i -> block.fragments[i].bases, w.fragments) : w.core_last - w.core_first + 1
    progress === nothing || (progress.total_bases += sum(work_bases, work))
    results = Vector{Vector{ScanVariantRow}}(undef, length(work))
    next = Threads.Atomic{Int}(1)
    tasks = map(1:min(scan_threads, length(work))) do _
        Threads.@spawn begin
            worker = PrefixHashScanWorker(open(ctx.dbi.gi.filepath, "r"))
            worker_stats = prefix_hash_scan_worker_stats(stats)
            try
                while true
                    i = Threads.atomic_add!(next, 1)
                    i > length(work) && break
                    if prefix_hash_scan_active_mask(state) == 0
                        results[i] = ScanVariantRow[]
                    else
                        results[i] = scan_prefix_hash_scan_variant_work(ctx, worker,
                            data, work[i]..., guides, state, seen, worker_stats, genotype_rows)
                    end
                    report_prefix_hash_scan_progress!(progress, work_bases(work[i]))
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
    return reduce(vcat, results; init = ScanVariantRow[])
end
