function prepare_prefix_hash_scan_variants(path, dbi, lengths, selected, annotate;
    verbose = false, max_alleles = PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES)

    max_alleles >= 1 || error("vcf_max_alleles must be at least 1")
    max_alleles > PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES && @warn(
        "vcf_max_alleles=$max_alleles exceeds the default $PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES. " *
        "Fragments per variant grow like C(n, k) in dense loci; " * PREFIX_HASH_SCAN_VCF_MAX_ALLELES_NOTE)
    started = time_ns()
    variants, counts = read_scan_variants(path, dbi, lengths, selected, annotate)
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
    for ci in eachindex(lengths), anchor in by_chrom[ci]
        anchor_first = length(fragments) + 1
        initial = variants[anchor]
        left = max(1, initial.first - window)
        path = ScanVariantSegment[]
        left < initial.first && push!(path, ScanVariantSegment(0, left, initial.first - 1))
        push!(path, ScanVariantSegment(anchor, 1, length(initial.sequence)))
        taken = Int[anchor]
        initial_length = initial.first - left + length(initial.sequence)
        goal = initial_length + window
        capped = false
        function emit(cursor, emitted)
            remaining = min(goal - emitted, lengths[ci] - cursor + 1)
            remaining > 0 && push!(path, ScanVariantSegment(0, cursor, cursor + remaining - 1))
            first = length(segments) + 1
            append!(segments, path)
            push!(fragments, ScanVariantFragment(ci, first:length(segments), anchor, emitted + max(remaining, 0)))
            remaining > 0 && pop!(path)
        end
        function walk(next, cursor, emitted)
            if emitted >= goal || next > last(by_chrom[ci])
                emit(cursor, emitted)
                return
            end
            v = variants[next]
            if v.first - cursor >= goal - emitted
                emit(cursor, emitted)
                return
            end
            # Skipping leaves reference emission to the next taken variant.
            walk(next + 1, cursor, emitted)
            v.first < cursor && return
            any(i -> scan_variants_conflict(v, variants[i]), taken) && return
            if length(taken) >= max_alleles
                capped = true
                return
            end
            n = length(path)
            cursor < v.first && push!(path, ScanVariantSegment(0, cursor, v.first - 1))
            push!(path, ScanVariantSegment(next, 1, length(v.sequence)))
            push!(taken, next)
            walk(next + 1, v.stop, emitted + v.first - cursor + length(v.sequence))
            pop!(taken)
            resize!(path, n)
        end
        walk(anchor + 1, initial.stop, initial_length)
        push!(per_anchor, length(fragments) - anchor_first + 1)
        capped && push!(truncated, anchor)
    end
    data = PrefixHashScanVariants(variants, by_chrom, max_stop, segments, fragments, counts,
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
        @info "prefixHashScan variants" counts=counts variants=length(variants) fragments=length(fragments) max_anchor_fragments=data.max_anchor_fragments max_alleles truncated_anchors=length(truncated) fragment_bases=sum(f -> f.bases, fragments; init=0) descriptor_bytes=Base.summarysize(data) enumeration_s=data.enumeration_ns / 1e9
    end
    return data
end
