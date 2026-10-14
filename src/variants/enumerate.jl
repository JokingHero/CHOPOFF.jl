function prepare_prefix_hash_scan_variants(path, dbi, lengths, selected, annotate;
    verbose = false, max_alleles = PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES,
    progress_interval = PREFIX_HASH_SCAN_PROGRESS_INTERVAL_S)

    max_alleles >= 1 || error("vcf_max_alleles must be at least 1")
    max_alleles > PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES && @warn(
        "vcf_max_alleles=$max_alleles exceeds the default $PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES. " *
        "Fragments per variant grow like C(n, k) in dense loci; " * PREFIX_HASH_SCAN_VCF_MAX_ALLELES_NOTE)
    started = time_ns()
    progress = verbose ? PrefixHashScanProgress("VCF read", progress_interval) : nothing
    variants, counts = read_scan_variants(path, dbi, lengths, selected, annotate; progress)
    progress = verbose ? PrefixHashScanProgress("VCF enumerate", progress_interval) : nothing
    progress === nothing || (progress.total_bases = length(variants))
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
    segments = ScanVariantSegment[]
    fragments = ScanVariantFragment[]
    window = length(dbi.motif) + dbi.motif.distance
    per_anchor = Int[]
    truncated = Int[]
    reference_spans = [Tuple{Int, Int}[] for _ in eachindex(lengths)]
    for ci in eachindex(lengths), anchor in by_chrom[ci]
        anchor_first = length(fragments) + 1
        anchor_segments = length(segments) + 1
        initial = variants[anchor]
        left = max(1, initial.first - window)
        path = ScanVariantSegment[]
        left < initial.first && push!(path, ScanVariantSegment(0, left, initial.first - 1))
        push!(path, ScanVariantSegment(anchor, 1, length(initial.sequence)))
        taken = Int[anchor]
        initial_length = initial.first - left + length(initial.sequence)
        walk = (; variants, segments, fragments, path, taken, ci, anchor, max_alleles,
            last_variant = last(by_chrom[ci]), goal = initial_length + window, chrom_length = lengths[ci])
        capped = walk_variant_anchor!(walk, anchor + 1, initial.stop, initial_length)
        # Merging per anchor keeps the span list near one entry per anchor.
        append!(reference_spans[ci], merge_variant_reference_spans!([(seg.first, seg.last)
            for seg in @view(segments[anchor_segments:end]) if seg.variant == 0 && seg.first <= seg.last]))
        report_prefix_hash_scan_progress!(progress, 1)
        push!(per_anchor, length(fragments) - anchor_first + 1)
        capped && push!(truncated, anchor)
    end
    reference = read_variant_reference_spans(dbi, lengths, reference_spans)
    data = PrefixHashScanVariants(variants, by_chrom, max_stop, segments, fragments, reference..., counts,
        maximum(per_anchor; init = 0), max_alleles, length(truncated), time_ns() - started)
    if !isempty(truncated)
        locus(i) = string(dbi.gi.chrom[variants[i].chrom], ":", variants[i].pos)
        @warn "VCF allele combinations truncated at vcf_max_alleles=$max_alleles; output omits local haplotypes with more ALT alleles" anchors=length(truncated) examples=unique(locus.(first(truncated, 5)))
    end
    dense = findall(>(PREFIX_HASH_SCAN_VARIANT_FRAGMENT_WARN), per_anchor)
    if !isempty(dense)
        worst = first(sort(dense; by = i -> -per_anchor[i]), 5)
        @warn "VCF allele combinations are dense; preparation time and memory grow with fragments" anchors=length(dense) fragments=length(fragments) worst=[string(dbi.gi.chrom[variants[i].chrom], ":", variants[i].pos) => per_anchor[i] for i in worst]
    end
    ignored = counts[:no_alt] + counts[:spanning_deletion] + counts[:non_ref]
    ignored > 0 && @warn "VCF non-sequence alleles ignored" no_alt=counts[:no_alt] spanning_deletion=counts[:spanning_deletion] non_ref=counts[:non_ref]
    if verbose
        @info "prefixHashScan variants" counts=counts variants=length(variants) fragments=length(fragments) max_anchor_fragments=data.max_anchor_fragments max_alleles truncated_anchors=length(truncated) fragment_bases=sum(f -> f.bases, fragments; init=0) descriptor_bytes=variant_descriptor_bytes(data) enumeration_s=data.enumeration_ns / 1e9
    end
    return data
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
    length(w.taken) >= w.max_alleles && return true
    n = length(w.path)
    cursor < v.first && push!(w.path, ScanVariantSegment(0, cursor, v.first - 1))
    push!(w.path, ScanVariantSegment(next, 1, length(v.sequence)))
    push!(w.taken, next)
    capped |= walk_variant_anchor!(w, next + 1, v.stop, emitted + v.first - cursor + length(v.sequence))
    pop!(w.taken)
    resize!(w.path, n)
    return capped
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

# Reads every reference base a fragment can use once, so workers and guide
# batches materialize fragments without file reads.
function read_variant_reference_spans(dbi, lengths, reference_spans)
    index = dbi.gi.is_fa ? FASTA.Index(dbi.gi.filepath * ".fai") :
        read_prefix_hash_scan_twobit_index(dbi.gi.filepath)
    firsts = Vector{Int}[]
    offsets = Vector{Int}[]
    bytes = UInt8[]
    merged_spans = merge_variant_reference_spans!.(reference_spans)
    sizehint!(bytes, sum(spans -> sum(((a, b),) -> b - a + 1, spans; init = 0), merged_spans; init = 0))
    open(dbi.gi.filepath, "r") do genome
        cursor = VariantReferenceCursor(genome, index, lengths)
        for (ci, merged) in enumerate(merged_spans)
            push!(firsts, first.(merged))
            push!(offsets, Int[])
            for (a, b) in merged
                push!(offsets[end], length(bytes))
                append!(bytes, variant_reference_view(cursor, ci, a, b))
            end
        end
    end
    return firsts, offsets, bytes
end
