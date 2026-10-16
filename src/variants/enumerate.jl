function prepare_prefix_hash_scan_variants(path, dbi, lengths, selected, annotate;
    verbose = false, max_alleles = PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES,
    progress_interval = PREFIX_HASH_SCAN_PROGRESS_INTERVAL_S, genotype = false, samples = nothing,
    intervals = nothing)

    max_alleles >= 1 || error("vcf_max_alleles must be at least 1")
    # Genotype walks follow carried haplotypes, so a large cap is the default there.
    !genotype && max_alleles > PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES && @warn(
        "vcf_max_alleles=$max_alleles exceeds the default $PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES. " *
        "Fragments per variant grow like C(n, k) in dense loci; " * PREFIX_HASH_SCAN_VCF_MAX_ALLELES_NOTE)
    started = time_ns()
    progress = verbose ? PrefixHashScanProgress("VCF read", progress_interval) : nothing
    variants, counts, genotypes = read_scan_variants(path, dbi, lengths, selected, annotate;
        progress, genotype, samples)
    by_chrom = UnitRange{Int}[]
    max_stop = Vector{Int}[]
    cursor = 1
    for ci in eachindex(lengths)
        first = cursor
        while cursor <= length(variants) && variants[cursor].chrom == ci
            cursor += 1
        end
        push!(by_chrom, first:(cursor - 1))
        push!(max_stop, accumulate(max, [v.stop for v in @view variants[first:(cursor - 1)]]))
    end
    genotypes === nothing || (counts[:genotype_conflicts] = check_variant_genotype_conflicts(variants, by_chrom, genotypes, dbi.gi.chrom))
    window = length(dbi.motif) + dbi.motif.distance
    anchors = [intervals === nothing ? collect(by_chrom[ci]) :
        variant_reachable_anchors(variants, by_chrom[ci], intervals[ci], window, lengths[ci])
        for ci in eachindex(lengths)]
    ignored = counts[:no_alt] + counts[:spanning_deletion] + counts[:non_ref]
    ignored > 0 && @warn "VCF non-sequence alleles ignored" no_alt=counts[:no_alt] spanning_deletion=counts[:spanning_deletion] non_ref=counts[:non_ref]
    data = PrefixHashScanVariants(variants, by_chrom, max_stop, anchors, intervals, counts,
        max_alleles, time_ns() - started, genotypes)
    if verbose
        @info "prefixHashScan variants" counts=counts variants=length(variants) anchors=sum(length, anchors; init = 0) max_alleles descriptor_bytes=variant_descriptor_bytes(data) read_s=data.read_ns / 1e9
    end
    return data
end

# Anchors whose walk can produce a site overlapping a BED interval. After the
# anchor a walk emits at most `window` bases; reference bases outside every
# variant footprint can only be passed by emitting them, so the walk ends before
# the (window + 1)-th such base. Deletions therefore extend the reach without a
# fixed margin, and no allele cap is assumed.
function variant_reachable_anchors(variants, range, intervals, window, chrom_length)
    isempty(intervals) && return Int[]
    starts, ends, covered = Int[], Int[], Int[0]
    for i in range
        v = variants[i]
        v.first < v.stop || continue
        if !isempty(ends) && v.first <= ends[end] + 1
            ends[end] = max(ends[end], v.stop - 1)
        else
            push!(starts, v.first)
            push!(ends, v.stop - 1)
        end
    end
    for j in eachindex(starts)
        push!(covered, covered[end] + ends[j] - starts[j] + 1)
    end
    # Reference bases in 1:x outside every footprint.
    function free(x)
        j = searchsortedlast(starts, x)
        return x - (j == 0 ? 0 : covered[j] + min(x, ends[j]) - starts[j] + 1)
    end
    kept = Int[]
    for i in range
        v = variants[i]
        lower = min(v.first - window, v.pos)
        target = free(v.stop - 1) + window + 1
        reach = if free(chrom_length) < target
            chrom_length
        else
            lo, hi = v.stop, chrom_length
            while lo < hi
                mid = (lo + hi) >>> 1
                free(mid) >= target ? (hi = mid) : (lo = mid + 1)
            end
            lo
        end
        prefix_hash_scan_overlaps(intervals, lower, reach) && push!(kept, i)
    end
    return kept
end

# Enumerates the walks of one chromosome's anchors in blocks of consecutive
# anchors on up to `threads` tasks. The variant scan reads the blocks directly;
# concatenating them would hold two copies of the chromosome's fragments.
function enumerate_variant_chromosome(data, dbi, lengths, ci, threads; progress = nothing)
    anchors = data.anchors[ci]
    window = length(dbi.motif) + dbi.motif.distance
    ploidy = data.genotypes === nothing ? 0 : variant_max_ploidy(data.genotypes)
    ranges = collect(Iterators.partition(anchors, max(1, cld(length(anchors), 4 * threads))))
    blocks = Vector{ScanVariantBlock}(undef, length(ranges))
    next = Threads.Atomic{Int}(1)
    tasks = map(1:min(threads, length(ranges))) do _
        Threads.@spawn while true
            b = Threads.atomic_add!(next, 1)
            b > length(ranges) && break
            blocks[b] = enumerate_variant_block(data, dbi, lengths, ranges[b], ci, window, ploidy, progress)
        end
    end
    foreach(wait, tasks)
    return blocks
end

function enumerate_variant_block(data, dbi, lengths, anchors, ci, window, ploidy, progress)
    variants, genotypes = data.variants, data.genotypes
    chrom_length = lengths[ci]
    segments = ScanVariantSegment[]
    fragments = ScanVariantFragment[]
    spans = Tuple{Int, Int}[]
    per_anchor = Int[]
    truncated = Int[]
    for anchor in anchors
        anchor_first = length(fragments) + 1
        anchor_segments = length(segments) + 1
        initial = variants[anchor]
        left = max(1, initial.first - window)
        path = ScanVariantSegment[]
        left < initial.first && push!(path, ScanVariantSegment(0, left, initial.first - 1))
        push!(path, ScanVariantSegment(anchor, 1, length(initial.sequence)))
        taken = Int[anchor]
        initial_length = initial.first - left + length(initial.sequence)
        live = genotypes === nothing ? nothing :
            [variant_carrier_slots(genotypes, initial.source, initial.alt_index, ploidy)]
        walk = (; variants, segments, fragments, path, taken, ci, anchor, max_alleles = data.max_alleles,
            last_variant = last(data.by_chrom[ci]), goal = initial_length + window, chrom_length,
            genotypes, live, ploidy)
        capped = walk_variant_anchor!(walk, variant_first_at(variants, anchor + 1, walk.last_variant, initial.stop),
            initial.stop, initial_length)
        # Merging per anchor keeps the span list near one entry per anchor.
        append!(spans, merge_variant_reference_spans!([(seg.first, seg.last)
            for seg in @view(segments[anchor_segments:end]) if seg.variant == 0 && seg.first <= seg.last]))
        progress === nothing || report_prefix_hash_scan_progress!(progress, 1)
        push!(per_anchor, length(fragments) - anchor_first + 1)
        capped && push!(truncated, anchor)
    end
    reference = read_variant_reference_spans(dbi, lengths, ci, merge_variant_reference_spans!(spans))
    return ScanVariantBlock(ci, segments, fragments, reference..., per_anchor, truncated)
end

function emit_variant_fragment!(w, cursor, emitted)
    remaining = min(w.goal - emitted, w.chrom_length - cursor + 1)
    remaining > 0 && push!(w.path, ScanVariantSegment(0, cursor, cursor + remaining - 1))
    first = length(w.segments) + 1
    append!(w.segments, w.path)
    push!(w.fragments, ScanVariantFragment(w.ci, first:length(w.segments), w.anchor, emitted + max(remaining, 0)))
    remaining > 0 && pop!(w.path)
    return nothing
end

# Top-level recursion: a recursive local closure is boxed and dispatches
# dynamically on every call. Returns true when vcf_max_alleles cut a combination.
function walk_variant_anchor!(w, next, cursor, emitted)
    if emitted >= w.goal || next > w.last_variant
        emit_variant_fragment!(w, cursor, emitted)
        return false
    end
    v = w.variants[next]
    if v.first - cursor >= w.goal - emitted
        emit_variant_fragment!(w, cursor, emitted)
        return false
    end
    # Skipping leaves reference emission to the next taken variant.
    capped = walk_variant_anchor!(w, next + 1, cursor, emitted)
    v.first < cursor && return capped
    any(i -> scan_variants_conflict(v, w.variants[i]), w.taken) && return capped
    # Genotype walks take only alleles some copy can carry with all taken ones.
    live = w.live === nothing ? nothing :
        variant_live_slots(w.live[end], variant_carrier_slots(w.genotypes, v.source, v.alt_index, w.ploidy))
    live === nothing || !isempty(live) || return capped
    length(w.taken) >= w.max_alleles && return true
    live === nothing || push!(w.live, live)
    n = length(w.path)
    cursor < v.first && push!(w.path, ScanVariantSegment(0, cursor, v.first - 1))
    # Sites overlapping the anchor end before the goal, so a long ALT (a tandem
    # copy) is cut there like the reference tail.
    used = min(length(v.sequence), w.goal - emitted - (v.first - cursor))
    push!(w.path, ScanVariantSegment(next, 1, used))
    push!(w.taken, next)
    # Variants starting inside a taken deletion are conflicts or behind the
    # cursor; skipping them one call at a time recursed once per variant.
    capped |= walk_variant_anchor!(w, variant_first_at(w.variants, next + 1, w.last_variant, v.stop),
        v.stop, emitted + v.first - cursor + used)
    pop!(w.taken)
    w.live === nothing || pop!(w.live)
    resize!(w.path, n)
    return capped
end

# First index in lo:hi whose variant starts at or after `stop` (hi + 1 if none);
# variants are sorted by `first` within a chromosome.
function variant_first_at(variants, lo, hi, stop)
    while lo <= hi
        mid = (lo + hi) >>> 1
        variants[mid].first < stop ? (lo = mid + 1) : (hi = mid - 1)
    end
    return lo
end

function merge_variant_reference_spans!(spans)
    sort!(spans)
    merged = Tuple{Int, Int}[]
    for (a, b) in spans
        if !isempty(merged) && a <= merged[end][2] + 1
            merged[end] = (merged[end][1], max(merged[end][2], b))
        else
            push!(merged, (a, b))
        end
    end
    return merged
end

# Reads every reference base a block's fragments can use once, so workers
# materialize fragments without file reads.
function read_variant_reference_spans(dbi, lengths, ci, merged)
    index = dbi.gi.is_fa ? FASTA.Index(dbi.gi.filepath * ".fai") :
        read_prefix_hash_scan_twobit_index(dbi.gi.filepath)
    offsets = Int[]
    bytes = UInt8[]
    sizehint!(bytes, sum(((a, b),) -> b - a + 1, merged; init = 0))
    open(dbi.gi.filepath, "r") do genome
        cursor = VariantReferenceCursor(genome, index, lengths)
        for (a, b) in merged
            push!(offsets, length(bytes))
            append!(bytes, variant_reference_view(cursor, ci, a, b))
        end
    end
    return first.(merged), offsets, bytes
end
