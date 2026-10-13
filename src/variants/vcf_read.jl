struct ScanVariant
    chrom::Int
    pos::Int
    first::Int                 # half-open reference footprint, or insertion gap
    stop::Int
    source::Int                # original record; its ALT choices are exclusive
    ref::String
    alt::String
    sequence::Vector{UInt8}     # trimmed, concrete replacement
    allele_key::String         # shared by the concrete realizations of an ALT
    annotation::String
end

# Variant records and reference spans arrive sorted within a chromosome, so one
# forward block read replaces a seek per record.
const PREFIX_HASH_SCAN_VARIANT_REFERENCE_BLOCK = 1 << 20

mutable struct VariantReferenceCursor
    io::IOStream
    index::Union{FASTA.Index, PrefixHashScanTwoBitIndex}
    lengths::Vector{Int}
    chrom::Int
    first::Int
    bytes::Vector{UInt8}
end

VariantReferenceCursor(io, index, lengths) = VariantReferenceCursor(io, index, lengths, 0, 0, UInt8[])

function variant_reference_view(cursor::VariantReferenceCursor, chrom, a, b)
    if cursor.chrom != chrom || a < cursor.first || b >= cursor.first + length(cursor.bytes)
        stop = min(cursor.lengths[chrom], max(b, a + PREFIX_HASH_SCAN_VARIANT_REFERENCE_BLOCK - 1))
        cursor.bytes = read_prefix_hash_scan_range!(cursor.bytes, cursor.io, cursor.index, chrom, a, stop)
        cursor.chrom, cursor.first = chrom, a
    end
    return @view cursor.bytes[(a - cursor.first + 1):(b - cursor.first + 1)]
end

function prefix_hash_scan_chromosomes(region, names)
    region === nothing && return trues(length(names))
    requested = strip.(split(region, ','))
    all(name -> name in names, requested) || error(
        "region accepts reference chromosome names only (comma-separated); unknown chromosome or interval: $region")
    return BitVector(name in requested for name in names)
end

variant_token(s) = replace(string(s), "%" => "%25", ";" => "%3B", "|" => "%7C")

function foreach_concrete_variant(f, sequence::String)
    bytes = collect(codeunits(sequence))
    positions = findall(b -> prefix_hash_scan_raw_code(b) == 0xff, bytes)
    masks = prefix_hash_scan_iupac_mask.(bytes[positions])
    any(iszero, masks) && error("Unsupported VCF allele sequence: $sequence")
    function visit(i)
        if i > length(positions)
            f(copy(bytes))
            return
        end
        for (bit, base) in enumerate(codeunits("ACGT"))
            masks[i] & (UInt8(1) << (bit - 1)) == 0 && continue
            bytes[positions[i]] = base
            visit(i + 1)
        end
    end
    visit(1)
end

function variant_trim(ref, alt)
    prefix = 0
    while prefix < min(length(ref), length(alt)) &&
            ref[prefix + 1] == alt[prefix + 1] && ref[prefix + 1] in "ACGT"
        prefix += 1
    end
    suffix = 0
    while suffix < min(length(ref), length(alt)) - prefix &&
            ref[end - suffix] == alt[end - suffix] && ref[end - suffix] in "ACGT"
        suffix += 1
    end
    # A no-op allele is still an annotated choice in the input.
    prefix == length(ref) == length(alt) && return 0, 0
    return prefix, suffix
end

function variant_annotation_fields(header, requested)
    fields = Tuple{String, String, String}[]
    for field in requested
        parts = split(field, '/'; limit = 2)
        length(parts) == 2 && parts[1] in ("INFO", "FORMAT") ||
            error("vcf_annotate entries must be INFO/key or FORMAT/key: $field")
        kind, key = String.(parts)
        definitions = filter(m -> m["ID"] == key, findall(header, kind))
        isempty(definitions) && error("VCF header has no $field field")
        push!(fields, (kind, key, definitions[1]["Number"]))
    end
    return fields
end

function variant_allele_annotation(record, fields, samples, alt_idx)
    isempty(fields) && return ""
    infos = VCF.hasinfo(record) ? Dict(VCF.info(record)) : Dict{String, String}()
    formats = VCF.hasformat(record) ? VCF.format(record) : String[]
    parts = String[]
    function select_value(value, number)
        value == "." && return value
        values = split(value, ',')
        number == "A" && return String(values[alt_idx])
        number == "R" && return join((values[1], values[alt_idx + 1]), ',')
        return value
    end
    for (kind, key, number) in fields
        value = if kind == "INFO"
            select_value(get(infos, key, "."), number)
        else
            join((sample * "=" * (
                key in formats ? select_value(VCF.genotype(record, i, key), number) : ".")
                for (i, sample) in enumerate(samples)), ',')
        end
        push!(parts, "$kind/$key=" * variant_token(value))
    end
    return "|" * join(parts, '|')
end

function read_scan_variants(path, dbi, lengths, selected, annotate; progress = nothing)
    variants = ScanVariant[]
    counts = Dict(:records => 0, :outside_region => 0, :no_alt => 0,
        :spanning_deletion => 0, :non_ref => 0, :multiallelic => 0, :ambiguous => 0, :iupac_skipped => 0)
    iupac_examples = String[]
    index = dbi.gi.is_fa ? FASTA.Index(dbi.gi.filepath * ".fai") :
        read_prefix_hash_scan_twobit_index(dbi.gi.filepath)
    chrom_index = Dict(name => i for (i, name) in enumerate(dbi.gi.chrom))
    stream = open(path, "r")
    reader = nothing
    try
        endswith(lowercase(path), ".gz") && (stream = GzipDecompressorStream(stream))
        reader = VCF.Reader(stream)
        header = VariantCallFormat.header(reader)
        fields = variant_annotation_fields(header, annotate)
        seen_chrom = Set{String}()
        previous_chrom = ""
        previous_pos = 0
        # Progress counts reference bases of finished chromosomes plus POS,
        # since VCF chromosome blocks need not follow reference order.
        progress === nothing || (progress.total_bases = sum(lengths; init = 0))
        finished_bases = 0
        reported = 0
        open(dbi.gi.filepath, "r") do genome
            cursor = VariantReferenceCursor(genome, index, lengths)
            for record in reader
                counts[:records] += 1
                source = counts[:records]
                chrom = VCF.chrom(record)
                if !haskey(chrom_index, chrom)
                    renamed = startswith(chrom, "chr") ? chrom[4:end] : "chr" * chrom
                    error("VCF chromosome is absent from the reference: $chrom; " * (haskey(chrom_index, renamed) ?
                        "the reference names it $renamed, rename with `bcftools annotate --rename-chrs`" :
                        "drop records on contigs absent from the reference, e.g. `bcftools view -t ^$chrom`"))
                end
                ci = chrom_index[chrom]
                pos = VCF.pos(record)
                if chrom != previous_chrom
                    chrom in seen_chrom && error("VCF chromosome blocks are not contiguous: $chrom")
                    isempty(previous_chrom) || (finished_bases += lengths[chrom_index[previous_chrom]])
                    push!(seen_chrom, chrom)
                    previous_chrom, previous_pos = chrom, 0
                end
                pos >= previous_pos || error("VCF positions are not sorted at $chrom:$pos")
                previous_pos = pos
                report_prefix_hash_scan_progress!(progress, finished_bases + pos - reported)
                reported = finished_bases + pos
                if !selected[ci]
                    counts[:outside_region] += 1
                    continue
                end
                ref = uppercase(VCF.ref(record))
                1 <= pos && pos + length(ref) - 1 <= lengths[ci] ||
                    error("VCF REF is outside the reference at $chrom:$pos")
                reference = uppercase(String(variant_reference_view(cursor, ci, pos, pos + length(ref) - 1)))
                for (a, b) in zip(codeunits(ref), codeunits(reference))
                    ma, mb = prefix_hash_scan_iupac_mask(a), prefix_hash_scan_iupac_mask(b)
                    mb != 0 && ma & mb == mb || error("VCF REF mismatch at $chrom:$pos: $ref versus $reference; " *
                        "check the reference build and normalize with `bcftools norm -f reference.fa -c e`")
                end
                alts = VCF.hasalt(record) ? VCF.alt(record) : String[]
                isempty(alts) && (counts[:no_alt] += 1)
                length(alts) > 1 && (counts[:multiallelic] += 1)
                # The parser can expose cumulative ranges for multiple IDs.
                # Reading their complete field preserves the original ID list.
                id_field = VCF.hasid(record) ? String(record.data[
                    first(first(record.id)):last(last(record.id))]) : ""
                for (alt_idx, alt_raw) in enumerate(alts)
                    if alt_raw in (".", "*", "<NON_REF>", "<*>")
                        reason = alt_raw == "." ? :no_alt : alt_raw == "*" ? :spanning_deletion : :non_ref
                        counts[reason] += 1
                        continue
                    end
                    alt = uppercase(alt_raw)
                    all(b -> prefix_hash_scan_iupac_mask(b) != 0, codeunits(alt)) ||
                        error("Unsupported VCF ALT at $chrom:$pos: $alt_raw; only sequence-resolved alleles are supported")
                    ambiguous = any(b -> prefix_hash_scan_raw_code(b) == 0xff, codeunits(alt))
                    # Long runs of N encode unknown sequence; their expansion is 4^n.
                    if prod(count_ones(prefix_hash_scan_iupac_mask(b)) for b in codeunits(alt)) >
                            PREFIX_HASH_SCAN_MAX_GUIDE_EXPANSIONS
                        counts[:iupac_skipped] += 1
                        length(iupac_examples) < 5 && push!(iupac_examples, "$chrom:$pos")
                        continue
                    end
                    counts[:ambiguous] += ambiguous
                    prefix, suffix = variant_trim(reference, alt)
                    first = pos + prefix
                    stop = pos + length(ref) - suffix
                    allele_key = variant_token("$chrom:$pos:$ref>$alt")
                    id = isempty(id_field) ? "$chrom:$pos:$ref>$alt" : id_field
                    annotation = variant_token(id) * "|" * allele_key
                    extra = variant_allele_annotation(record, fields, header.sampleID, alt_idx)
                    foreach_concrete_variant(alt) do concrete
                        resolved = ambiguous ? "|resolved=" * String(copy(concrete)) : ""
                        push!(variants, ScanVariant(ci, pos, first, stop, source, ref, alt,
                            concrete[(prefix + 1):(end - suffix)], allele_key, annotation * resolved * extra))
                    end
                end
            end
        end
    finally
        reader === nothing ? close(stream) : close(reader)
    end
    counts[:iupac_skipped] > 0 && @warn "VCF ALT alleles with more than $PREFIX_HASH_SCAN_MAX_GUIDE_EXPANSIONS IUPAC realizations ignored" alleles=counts[:iupac_skipped] examples=iupac_examples
    sort!(variants; by = v -> (v.chrom, v.first, v.first == v.stop ? 0 : 1, v.source, v.annotation))
    return variants, counts
end

function scan_variants_conflict(a::ScanVariant, b::ScanVariant)
    a.source == b.source && return true
    a.chrom == b.chrom || return false
    agap, bgap = a.first == a.stop, b.first == b.stop
    agap && bgap && return a.first == b.first
    agap && return b.first < a.first < b.stop
    bgap && return a.first < b.first < a.stop
    return max(a.first, b.first) < min(a.stop, b.stop)
end
