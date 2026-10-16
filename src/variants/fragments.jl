# Fragments per anchor double with each compatible variant in reach; warn above this.
const PREFIX_HASH_SCAN_VARIANT_FRAGMENT_WARN = 4096

# Maximum ALT alleles in one local haplotype (walk), counting the anchor. Exact
# enumeration is infeasible for cohort VCFs; see vcf.md for chr22 measurements.
const PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES = 3
const PREFIX_HASH_SCAN_VCF_MAX_ALLELES_NOTE =
    "on GRCh38 chr22, gnomAD AF>=1% needs 3.4M, 30M and ~252M fragments at k=2, 3 and 4 " *
    "(about 0.3 GB per million fragments; ~80 GB at k=4); memory follows the largest selected chromosome. " *
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

# Prepared once per search and reused across guide batches. Fragments are
# enumerated per chromosome during the variant scan (ScanVariantBlock), so
# peak memory follows the largest chromosome rather than the whole selection.
struct PrefixHashScanVariants
    variants::Vector{ScanVariant}
    by_chrom::Vector{UnitRange{Int}}
    max_stop::Vector{Vector{Int}}
    anchors::Vector{Vector{Int}}            # per chromosome: anchors whose walks are enumerated
    intervals::Union{Nothing, Vector{Vector{UnitRange{Int}}}}   # BED region, merged per chromosome
    counts::Dict{Symbol, Int}
    max_alleles::Int
    read_ns::UInt64
    genotypes::Union{Nothing, ScanGenotypes}
end

# Fragments of consecutive anchors on one chromosome.
struct ScanVariantBlock
    chrom::Int
    segments::Vector{ScanVariantSegment}
    fragments::Vector{ScanVariantFragment}
    reference_first::Vector{Int}            # merged reference spans
    reference_offset::Vector{Int}           # zero-based span offsets in reference_bytes
    reference_bytes::Vector{UInt8}
    per_anchor::Vector{Int}                 # fragments per enumerated anchor
    truncated::Vector{Int}                  # anchors cut by vcf_max_alleles
end

# Base.summarysize tracks object identity and took over three minutes on a
# genome-wide single-sample VCF. This linear count omits per-object headers,
# about 20% of the total for GIAB HG001.
function variant_descriptor_bytes(data::PrefixHashScanVariants)
    nested(v) = sum(sizeof, v; init = 0)
    variant(v) = sizeof(ScanVariant) + sizeof(v.ref) + sizeof(v.alt) + sizeof(v.sequence) +
        sizeof(v.allele_key) + sizeof(v.annotation)
    return sum(variant, data.variants; init = 0) + sizeof(data.by_chrom) + nested(data.max_stop) +
        nested(data.anchors) + variant_genotype_bytes(data.genotypes)
end

variant_descriptor_bytes(c::ScanVariantBlock) = sizeof(c.segments) + sizeof(c.fragments) +
    sizeof(c.reference_first) + sizeof(c.reference_offset) + sizeof(c.reference_bytes)

variant_genotype_bytes(::Nothing) = 0
variant_genotype_bytes(g::ScanGenotypes) = sum(sizeof, g.samples; init = 0) + sizeof(g.record_source) +
    sizeof(g.record_first) + sizeof(g.sample) + sizeof(g.group) + sizeof(g.allele_first) + sizeof(g.alleles)

function materialize_variant_fragment!(raw, data, chrom, fragment, first, last)
    offset = 0
    for si in fragment.segments
        seg = chrom.segments[si]
        n = length(seg.first:seg.last)
        lo, hi = max(first, offset + 1), min(last, offset + n)
        if lo <= hi
            a, b = seg.first + lo - offset - 1, seg.first + hi - offset - 1
            if seg.variant == 0
                starts = chrom.reference_first
                j = searchsortedlast(starts, a)
                o = chrom.reference_offset[j] + a - starts[j]
                append!(raw, @view chrom.reference_bytes[(o + 1):(o + b - a + 1)])
            else
                append!(raw, @view data.variants[seg.variant].sequence[a:b])
            end
        end
        offset += n
        offset >= last && break
    end
    return raw
end

function lift_variant_position(data, chrom, fragment, position)
    offset = 0
    for si in fragment.segments
        seg = chrom.segments[si]
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

function variant_site_alleles(chrom, fragment, lo, hi)
    alleles = Int[]
    offset = 0
    for si in fragment.segments
        seg = chrom.segments[si]
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

function variant_anchor_overlaps(chrom, fragment, lo, hi)
    offset = 0
    for si in fragment.segments
        seg = chrom.segments[si]
        n = length(seg.first:seg.last)
        if seg.variant == fragment.anchor
            return n == 0 ? lo <= offset < hi : max(lo, offset + 1) <= min(hi, offset + n)
        end
        offset += n
    end
    return false
end

# Indexes of input alleles intersecting reference bases lo:hi; an insertion
# counts when its gap lies inside the site.
function variant_reference_overlap_indices(data, chrom, lo, hi)
    range = data.by_chrom[chrom]
    indices = Int[]
    isempty(range) && return indices
    # Prefix maxima include long deletions starting well before the site.
    local_first = searchsortedfirst(data.max_stop[chrom], lo + 1)
    for i in (first(range) + local_first - 1):last(range)
        v = data.variants[i]
        v.first > hi && break
        overlaps = v.first == v.stop ? lo < v.first <= hi : v.first <= hi && v.stop > lo
        overlaps && push!(indices, i)
    end
    return indices
end

function variant_reference_overlaps(data, chrom, lo, hi)
    annotations = [first(split(data.variants[i].annotation, '|'; limit = 2)) * "|" * data.variants[i].allele_key
        for i in variant_reference_overlap_indices(data, chrom, lo, hi)]
    return join(unique(annotations), ';')
end

variant_genotype_allele(v::ScanVariant) = (v.source, v.alt_index)

# Genotype requirements of a variant site: the taken alleles, and the untaken
# alleles a carrying copy must lack. Those intersect a consumed reference base
# or insert at a gap the site path crosses (gap g lies before base g). Alleles
# hidden inside a taken deletion impose nothing.
function variant_site_requirements(data, chrom, fragment, lo, hi)
    ranges = Tuple{Int, Int}[]
    gaps = Int[]
    taken = Int[]
    offset = 0
    exit = 0
    for si in fragment.segments
        seg = chrom.segments[si]
        n = length(seg.first:seg.last)
        if seg.variant == 0
            a, b = max(lo, offset + 1), min(hi, offset + n)
            if a <= b
                ra, rb = seg.first + a - offset - 1, seg.first + b - offset - 1
                push!(ranges, (ra, rb))
                exit == 0 || push!(gaps, exit, ra)
                append!(gaps, (ra + 1):rb)
                exit = rb + 1
            end
        else
            used = n == 0 ? lo <= offset < hi : max(lo, offset + 1) <= min(hi, offset + n)
            if used
                v = data.variants[seg.variant]
                push!(taken, seg.variant)
                exit == 0 || push!(gaps, exit, v.first)
                exit = v.stop
            end
        end
        offset += n
    end
    excluded = Tuple{Int, Int}[]
    bounds = vcat(first.(ranges), gaps)
    if !isempty(bounds)
        low, high = extrema(vcat(bounds, last.(ranges)))
        range = data.by_chrom[fragment.chrom]
        for i in (first(range) + searchsortedfirst(data.max_stop[fragment.chrom], low) - 1):last(range)
            v = data.variants[i]
            v.first > high && break
            i in taken && continue
            hit = v.first == v.stop ? v.first in gaps :
                any(((a, b),) -> v.first <= b && v.stop > a, ranges)
            hit && push!(excluded, variant_genotype_allele(v))
        end
    end
    return [variant_genotype_allele(data.variants[i]) for i in taken], excluded
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
const PREFIX_SCAN_GENOTYPE_HEADER = chomp(PREFIX_SCAN_VARIANT_HEADER) *
    ",carriers,lost,n_carriers,n_possible,n_lost,n_lost_possible\n"

prefix_scan_variant_header(data) = data.genotypes === nothing ? PREFIX_SCAN_VARIANT_HEADER : PREFIX_SCAN_GENOTYPE_HEADER

function write_variant_reference_hit(out, hit, guide, chrom_name, data, geometry, motif)
    spec = prefix_scan_matcher_spec(geometry.matcher)
    pos_offset = hit.is_antisense ? spec.rev_pos_offset : spec.fwd_pos_offset
    candidate = hit.pos - pos_offset
    lo, hi = prefix_hash_scan_site_span(geometry, motif, candidate, hit.is_antisense, hit.aln_ref)
    other = pos_offset == 0 ? hi : lo
    chrom = findfirst(==(chrom_name), data[2])
    overlaps = variant_reference_overlaps(data[1], chrom, lo, hi)
    row = (guide, hit.aln_guide, hit.aln_ref, hit.dist, chrom_name,
        hit.pos, hit.is_antisense ? "-" : "+", other, "", overlaps, "", "")
    genotypes = data[1].genotypes
    if genotypes === nothing
        write_variant_csv_row(out, row)
    else
        excluded = [variant_genotype_allele(data[1].variants[i])
            for i in variant_reference_overlap_indices(data[1], chrom, lo, hi)]
        lost, certain, possible = isempty(excluded) ? ("", 0, 0) :
            haplotype_tokens(genotypes, Tuple{Int, Int}[], excluded; lost = true)
        write_variant_csv_row(out, (row..., "", lost, "", "", certain, possible))
    end
end
