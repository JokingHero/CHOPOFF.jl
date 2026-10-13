# Fragments per anchor double with each compatible variant in reach; warn above this.
const PREFIX_HASH_SCAN_VARIANT_FRAGMENT_WARN = 4096

# Maximum ALT alleles in one local haplotype (walk), counting the anchor. Exact
# enumeration is infeasible for cohort VCFs; see vcf.md for chr22 measurements.
const PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES = 3
const PREFIX_HASH_SCAN_VCF_MAX_ALLELES_NOTE =
    "on GRCh38 chr22, gnomAD AF>=1% needs 3.4M, 30M and ~252M fragments at k=2, 3 and 4 " *
    "(about 0.3 GB per million fragments; ~80 GB at k=4), and genome-wide runs hold all chromosomes. " *
    "Single-sample VCFs stay small at any k. See vcf.md before raising k for cohort VCFs."

struct ScanVariantSegment
    variant::Int              # zero for reference; nonzero segments retain deletions
    first::Int
    last::Int
end

struct ScanVariantFragment
    chrom::Int
    segments::UnitRange{Int}
    anchor::Int
    bases::Int
end

struct PrefixHashScanVariants
    variants::Vector{ScanVariant}
    by_chrom::Vector{UnitRange{Int}}
    max_stop::Vector{Vector{Int}}
    segments::Vector{ScanVariantSegment}
    fragments::Vector{ScanVariantFragment}
    reference_first::Vector{Vector{Int}}    # merged reference spans per chromosome
    reference_offset::Vector{Vector{Int}}   # zero-based span offsets in reference_bytes
    reference_bytes::Vector{UInt8}
    counts::Dict{Symbol, Int}
    max_anchor_fragments::Int
    max_alleles::Int
    truncated_anchors::Int
    enumeration_ns::UInt64
end

# Base.summarysize tracks object identity and took over three minutes on a
# genome-wide single-sample VCF. This linear count omits per-object headers,
# about 20% of the total for GIAB HG001.
function variant_descriptor_bytes(data::PrefixHashScanVariants)
    nested(v) = sum(sizeof, v; init = 0)
    variant(v) = sizeof(ScanVariant) + sizeof(v.ref) + sizeof(v.alt) + sizeof(v.sequence) +
        sizeof(v.allele_key) + sizeof(v.annotation)
    return sum(variant, data.variants; init = 0) + sizeof(data.by_chrom) + nested(data.max_stop) +
        sizeof(data.segments) + sizeof(data.fragments) + nested(data.reference_first) +
        nested(data.reference_offset) + sizeof(data.reference_bytes)
end

function materialize_variant_fragment!(raw, data, fragment, first, last)
    offset = 0
    for si in fragment.segments
        seg = data.segments[si]
        n = length(seg.first:seg.last)
        lo, hi = max(first, offset + 1), min(last, offset + n)
        if lo <= hi
            a, b = seg.first + lo - offset - 1, seg.first + hi - offset - 1
            if seg.variant == 0
                starts = data.reference_first[fragment.chrom]
                j = searchsortedlast(starts, a)
                o = data.reference_offset[fragment.chrom][j] + a - starts[j]
                append!(raw, @view data.reference_bytes[(o + 1):(o + b - a + 1)])
            else
                append!(raw, @view data.variants[seg.variant].sequence[a:b])
            end
        end
        offset += n
        offset >= last && break
    end
    return raw
end

function lift_variant_position(data, fragment, position)
    offset = 0
    for si in fragment.segments
        seg = data.segments[si]
        n = length(seg.first:seg.last)
        if offset < position <= offset + n
            within = seg.first + position - offset - 1
            seg.variant == 0 && return within, ""
            v = data.variants[seg.variant]
            length(v.sequence) == v.stop - v.first && return v.first + within - 1, ""
            return v.pos, v.allele_key * ":" * string(within)
        end
        offset += n
    end
    error("Variant position lies outside its fragment: $position")
end

function variant_site_alleles(data, fragment, lo, hi)
    alleles = Int[]
    offset = 0
    for si in fragment.segments
        seg = data.segments[si]
        n = length(seg.first:seg.last)
        if seg.variant != 0
            # A deletion participates only if the site crosses its junction.
            used = n == 0 ? lo <= offset < hi : max(lo, offset + 1) <= min(hi, offset + n)
            used && push!(alleles, seg.variant)
        end
        offset += n
    end
    return alleles
end

function variant_anchor_overlaps(data, fragment, lo, hi)
    offset = 0
    for si in fragment.segments
        seg = data.segments[si]
        n = length(seg.first:seg.last)
        if seg.variant == fragment.anchor
            return n == 0 ? lo <= offset < hi : max(lo, offset + 1) <= min(hi, offset + n)
        end
        offset += n
    end
    return false
end

function variant_reference_overlaps(data, chrom, lo, hi)
    range = data.by_chrom[chrom]
    isempty(range) && return ""
    # Prefix maxima include long deletions starting well before the site.
    local_first = searchsortedfirst(data.max_stop[chrom], lo + 1)
    annotations = String[]
    for i in (first(range) + local_first - 1):last(range)
        v = data.variants[i]
        v.first > hi && break
        overlaps = v.first == v.stop ? lo < v.first <= hi : v.first <= hi && v.stop > lo
        overlaps && push!(annotations, first(split(v.annotation, '|'; limit = 2)) * "|" * v.allele_key)
    end
    return join(unique(annotations), ';')
end

# align returns no endpoints. Its non-gap reference length determines the
# consumed guide-side bases; PAM positions must be added separately.
function prefix_hash_scan_site_span(geometry, motif, candidate_start, antisense, aln_ref)
    consumed = count(!=('-'), aln_ref)
    positions = [candidate_start + prefix_scan_reference_offset(geometry.matcher, antisense, i)
        for i in 1:consumed]
    pam = antisense ? motif.pam_loci_rve : motif.pam_loci_fwd
    append!(positions, (candidate_start + i - 1 for i in pam))
    return extrema(positions)
end

function write_variant_csv_row(out, fields)
    for (i, value) in enumerate(fields)
        i == 1 || print(out, ',')
        text = string(value)
        if any(c -> c in (',', '"', '\n', '\r'), text)
            print(out, '"', replace(text, "\"" => "\"\""), '"')
        else
            print(out, text)
        end
    end
    print(out, '\n')
end

const PREFIX_SCAN_VARIANT_HEADER = "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand,end,variants,overlapping_variants,start_origin,end_origin\n"

function write_variant_reference_hit(out, hit, guide, chrom_name, data, geometry, motif)
    spec = prefix_scan_matcher_spec(geometry.matcher)
    pos_offset = hit.is_antisense ? spec.rev_pos_offset : spec.fwd_pos_offset
    candidate = hit.pos - pos_offset
    lo, hi = prefix_hash_scan_site_span(geometry, motif, candidate, hit.is_antisense, hit.aln_ref)
    other = pos_offset == 0 ? hi : lo
    chrom = findfirst(==(chrom_name), data[2])
    overlaps = variant_reference_overlaps(data[1], chrom, lo, hi)
    write_variant_csv_row(out, (guide, hit.aln_guide, hit.aln_ref, hit.dist, chrom_name,
        hit.pos, hit.is_antisense ? "-" : "+", other, "", overlaps, "", ""))
end
