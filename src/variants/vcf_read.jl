struct ScanVariant
    chrom::Int
    pos::Int
    first::Int                 # half-open reference footprint, or insertion gap
    stop::Int
    source::Int                # original record; its ALT choices are exclusive
    alt_index::Int             # ALT number in the original record, as in GT
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

# `region` is `nothing`, comma-separated chromosome names, or a BED file path.
# Returns the selected chromosomes and, for BED, 1-based inclusive intervals
# sorted and merged per chromosome (otherwise `nothing`).
function prefix_hash_scan_region(region, names, lengths)
    region === nothing && return trues(length(names)), nothing
    # A long chromosome list is not a valid path; stat then throws (ENAMETOOLONG).
    is_bed = try isfile(region) catch e; e isa Base.IOError || rethrow(); false end
    is_bed && return prefix_hash_scan_bed(region, names, lengths)
    requested = strip.(split(region, ','))
    all(name -> name in names, requested) || error(
        "region accepts reference chromosome names (comma-separated) or a BED file path; unknown chromosome or missing file: $region")
    return BitVector(name in requested for name in names), nothing
end

function prefix_hash_scan_bed(path, names, lengths)
    chrom_index = Dict(name => i for (i, name) in enumerate(names))
    intervals = [UnitRange{Int}[] for _ in names]
    for (line_number, line) in enumerate(eachline(path))
        text = strip(line)
        (isempty(text) || startswith(text, '#') || startswith(text, "track") || startswith(text, "browser")) && continue
        fields = split(text, '\t')
        length(fields) >= 3 || error("BED line $line_number needs chrom, start and end separated by tabs: $path")
        chrom = String(fields[1])
        ci = get(chrom_index, chrom, 0)
        if ci == 0
            renamed = startswith(chrom, "chr") ? chrom[4:end] : "chr" * chrom
            error("BED chromosome is absent from the reference: $chrom (line $line_number)" *
                (haskey(chrom_index, renamed) ? "; the reference names it $renamed" : ""))
        end
        a, b = tryparse(Int, fields[2]), tryparse(Int, fields[3])
        a !== nothing && b !== nothing && 0 <= a < b <= lengths[ci] ||
            error("BED line $line_number needs 0 <= start < end <= chromosome length: $text")
        push!(intervals[ci], (a + 1):b)
    end
    any(!isempty, intervals) || error("BED region has no intervals: $path")
    for list in intervals
        sort!(list; by = first)
        merged = UnitRange{Int}[]
        for r in list
            if !isempty(merged) && first(r) <= last(merged[end]) + 1
                merged[end] = first(merged[end]):max(last(merged[end]), last(r))
            else
                push!(merged, r)
            end
        end
        copy!(list, merged)
    end
    return BitVector(!isempty(list) for list in intervals), intervals
end

# True when lo:hi overlaps one of the sorted, disjoint `intervals`.
function prefix_hash_scan_overlaps(intervals, lo, hi)
    j = searchsortedfirst(intervals, lo:lo; by = last)
    return j <= length(intervals) && first(intervals[j]) <= hi
end

function variant_token(s)
    t = string(s)
    any(b -> b == UInt8('%') || b == UInt8(';') || b == UInt8('|'), codeunits(t)) || return t
    return replace(t, "%" => "%25", ";" => "%3B", "|" => "%7C")
end

function foreach_concrete_variant(f, sequence::String)
    all(b -> prefix_hash_scan_raw_code(b) != 0xff, codeunits(sequence)) &&
        return f(collect(codeunits(sequence)))
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

# Symbolic and breakend ALTs. Precise DEL, DUP and INV become literal alleles;
# the other kinds are skipped and counted. Returns `nothing` for sequence ALTs.
function variant_sv_kind(alt)
    any(c -> c in ('<', '[', ']', '.'), alt) || return nothing
    startswith(alt, '<') || return :bnd
    type = uppercase(alt)
    (type == "<DEL>" || startswith(type, "<DEL:")) && return :del
    type in ("<DUP>", "<DUP:TANDEM>") && return :dup
    startswith(type, "<DUP:") && return :dup_other
    type == "<INV>" && return :inv
    (type == "<INS>" || startswith(type, "<INS:")) && return :ins
    (type == "<CNV>" || startswith(type, "<CNV:")) && return :cnv
    return :other
end

function variant_info_column(text, tab5)
    tab6 = findnext('\t', text, tab5 + 1)
    tab7 = tab6 === nothing ? nothing : findnext('\t', text, tab6 + 1)
    tab7 === nothing && return SubString(text, 1, 0)
    tab8 = something(findnext('\t', text, tab7 + 1), ncodeunits(text) + 1)
    return SubString(text, tab7 + 1, tab8 - 1)
end

# Value of INFO `key`; "" for a flag, `nothing` when absent.
function variant_info_field(info, key)
    for item in eachsplit(info, ';')
        item == key && return ""
        startswith(item, key * "=") && return SubString(item, length(key) + 2)
    end
    return nothing
end

# Last reference base of a symbolic SV: INFO/END, else POS + |SVLEN| (per ALT).
function variant_sv_end(info, pos, alt_idx)
    value = variant_info_field(info, "END")
    value === nothing || return tryparse(Int, value)
    svlen = variant_info_field(info, "SVLEN")
    svlen === nothing && return nothing
    values = split(svlen, ',')
    n = tryparse(Int, values[length(values) == 1 ? 1 : min(alt_idx, length(values))])
    return n === nothing ? nothing : pos + abs(n)
end

const VARIANT_IUPAC_BY_MASK = codeunits("?ACMGRSVTWYHKDBN")

function variant_reverse_complement!(bytes)
    reverse!(bytes)
    for i in eachindex(bytes)
        mask = prefix_hash_scan_iupac_mask(bytes[i])
        mask == 0 || (bytes[i] = VARIANT_IUPAC_BY_MASK[prefix_hash_scan_complement_mask(mask) + 1])
    end
    return bytes
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

# Record lines per batch; about 13 MB for a 3,202-sample cohort line.
const PREFIX_HASH_SCAN_VCF_LINE_BATCH = 1024

function read_scan_variants(path, dbi, lengths, selected, annotate; progress = nothing,
    genotype = false, samples = nothing)
    variants = ScanVariant[]
    counts = Dict(:records => 0, :outside_region => 0, :no_alt => 0,
        :spanning_deletion => 0, :non_ref => 0, :multiallelic => 0, :ambiguous => 0, :iupac_skipped => 0)
    iupac_examples = String[]
    sv_examples = Dict{Symbol, Vector{String}}()
    index = dbi.gi.is_fa ? FASTA.Index(dbi.gi.filepath * ".fai") :
        read_prefix_hash_scan_twobit_index(dbi.gi.filepath)
    chrom_index = Dict(name => i for (i, name) in enumerate(dbi.gi.chrom))
    stream = open_variant_vcf(path)
    genome = open(dbi.gi.filepath, "r")
    store = nothing
    columns = Int[]
    try
        # Records are split by line so that unselected chromosomes skip the
        # full record parser; only CHROM and POS are read for them.
        lines = eachline(stream)
        header_lines = String[]
        line = iterate(lines)
        # Only INFO/FORMAT definitions and sample names are used. Other meta lines
        # are skipped, so malformed ones (`###INFO=` in the 1000 Genomes SV set)
        # cannot fail the header parser.
        while line !== nothing && startswith(line[1], '#')
            text = line[1]
            any(p -> startswith(text, p), ("##fileformat=", "##INFO=", "##FORMAT=", "#CHROM")) &&
                push!(header_lines, text)
            line = iterate(lines)
        end
        header = VariantCallFormat.header(VCF.Reader(IOBuffer(join(header_lines, '\n') * '\n')))
        fields = variant_annotation_fields(header, annotate)
        if genotype
            columns = variant_sample_columns(header.sampleID, samples)
            any(m -> m["ID"] == "GT", findall(header, "FORMAT")) ||
                error("VCF header has no FORMAT/GT; vcf_mode=:genotype needs GT, use vcf_mode=:union for sites-only VCFs")
            store = ScanGenotypes(header.sampleID)
            counts[:not_carried] = 0
        end
        seen_chrom = Set{String}()
        previous_chrom = ""
        previous_pos = 0
        # Progress counts reference bases of finished chromosomes plus POS,
        # since VCF chromosome blocks need not follow reference order.
        progress === nothing || (progress.total_bases = sum(lengths; init = 0))
        finished_bases = 0
        reported = 0
        cursor = VariantReferenceCursor(genome, index, lengths)
        # Batches let genotype columns, the bulk of cohort lines, parse on all threads.
        records = Iterators.flatten((line === nothing ? String[] : [line[1]], lines))
        for batch in Iterators.partition(records, PREFIX_HASH_SCAN_VCF_LINE_BATCH)
            sources = cumsum([!isempty(text) for text in batch]) .+ counts[:records]
            batch_carried = genotype ? push_variant_genotype_batch!(store, batch, sources,
                chrom_index, selected, columns) : nothing
            for (line_idx, text) in enumerate(batch)
                isempty(text) && continue
                counts[:records] += 1
                source = counts[:records]
                tab1 = findfirst('\t', text)
                tab2 = tab1 === nothing ? nothing : findnext('\t', text, tab1 + 1)
                pos = tab2 === nothing ? nothing : tryparse(Int, SubString(text, tab1 + 1, tab2 - 1))
                pos === nothing && error("Malformed VCF record $source: $(first(text, 80))")
                chrom = SubString(text, 1, tab1 - 1)
                if !haskey(chrom_index, chrom)
                    renamed = startswith(chrom, "chr") ? chrom[4:end] : "chr" * chrom
                    error("VCF chromosome is absent from the reference: $chrom; " * (haskey(chrom_index, renamed) ?
                        "the reference names it $renamed, rename with `bcftools annotate --rename-chrs`" :
                        "drop records on contigs absent from the reference, e.g. `bcftools view -t ^$chrom`"))
                end
                ci = chrom_index[chrom]
                if chrom != previous_chrom
                    chrom in seen_chrom && error("VCF chromosome blocks are not contiguous: $chrom")
                    isempty(previous_chrom) || (finished_bases += lengths[chrom_index[previous_chrom]])
                    push!(seen_chrom, chrom)
                    previous_chrom, previous_pos = String(chrom), 0
                end
                pos >= previous_pos || error("VCF positions are not sorted at $chrom:$pos")
                previous_pos = pos
                report_prefix_hash_scan_progress!(progress, finished_bases + pos - reported)
                reported = finished_bases + pos
                if !selected[ci]
                    counts[:outside_region] += 1
                    continue
                end
                # ID, REF and ALT come from the line; the full record parser runs
                # only when vcf_annotate needs INFO or FORMAT fields.
                tab3 = findnext('\t', text, tab2 + 1)
                tab4 = tab3 === nothing ? nothing : findnext('\t', text, tab3 + 1)
                tab4 === nothing && error("Malformed VCF record $source: $(first(text, 80))")
                tab5 = something(findnext('\t', text, tab4 + 1), ncodeunits(text) + 1)
                record = isempty(fields) ? nothing : VCF.Record(text)
                ref = uppercase(SubString(text, tab3 + 1, tab4 - 1))
                1 <= pos && pos + length(ref) - 1 <= lengths[ci] ||
                    error("VCF REF is outside the reference at $chrom:$pos")
                reference = uppercase(String(variant_reference_view(cursor, ci, pos, pos + length(ref) - 1)))
                for (a, b) in zip(codeunits(ref), codeunits(reference))
                    ma, mb = prefix_hash_scan_iupac_mask(a), prefix_hash_scan_iupac_mask(b)
                    mb != 0 && ma & mb == mb || error("VCF REF mismatch at $chrom:$pos: $ref versus $reference; " *
                        "check the reference build and normalize with `bcftools norm -f reference.fa -c e`")
                end
                alts = split(SubString(text, tab4 + 1, tab5 - 1), ',')
                carried = genotype ? batch_carried[line_idx] : UInt128(0)
                length(alts) > 1 && (counts[:multiallelic] += 1)
                id_field = SubString(text, tab2 + 1, tab3 - 1)
                id_field == "." && (id_field = SubString(text, 1, 0))
                for (alt_idx, alt_raw) in enumerate(alts)
                    if alt_raw in (".", "*", "<NON_REF>", "<*>")
                        reason = alt_raw == "." ? :no_alt : alt_raw == "*" ? :spanning_deletion : :non_ref
                        counts[reason] += 1
                        continue
                    end
                    if genotype && alt_idx < 127 && carried & (UInt128(1) << alt_idx) == 0
                        counts[:not_carried] += 1
                        continue
                    end
                    kind = variant_sv_kind(alt_raw)
                    if kind !== nothing
                        info = variant_info_column(text, tab5)
                        kind in (:del, :dup, :inv) && variant_info_field(info, "IMPRECISE") !== nothing &&
                            (kind = :imprecise)
                        if !(kind in (:del, :dup, :inv))
                            key = Symbol("sv_skipped_", kind)
                            counts[key] = get(counts, key, 0) + 1
                            examples = get!(sv_examples, kind, String[])
                            length(examples) < 5 && push!(examples, "$chrom:$pos")
                            continue
                        end
                        sv_end = variant_sv_end(info, pos, alt_idx)
                        sv_end === nothing && error("Symbolic VCF ALT at $chrom:$pos needs integer INFO/END or INFO/SVLEN: $alt_raw")
                        pos < sv_end <= lengths[ci] ||
                            error("Symbolic VCF ALT at $chrom:$pos ends outside (POS, chromosome end]: $alt_raw, END=$sv_end")
                        counts[Symbol("sv_", kind)] = get(counts, Symbol("sv_", kind), 0) + 1
                        uppercase(alt_raw) == "<DUP>" && (counts[:sv_dup_untyped] = get(counts, :sv_dup_untyped, 0) + 1)
                        # A tandem copy is inserted before the original, an inversion
                        # replaces it; reference ambiguity stays as in the reference.
                        sequence = kind == :del ? UInt8[] :
                            map(b -> UInt8(uppercase(Char(b))), variant_reference_view(cursor, ci, pos + 1, sv_end))
                        kind == :inv && variant_reverse_complement!(sequence)
                        allele_key = variant_token("$chrom:$pos:$ref>$alt_raw/$sv_end")
                        annotation = (isempty(id_field) ? allele_key : variant_token(id_field)) * "|" * allele_key
                        push!(variants, ScanVariant(ci, pos, pos + 1, kind == :dup ? pos + 1 : sv_end + 1, source, alt_idx,
                            ref, alt_raw, sequence, allele_key,
                            annotation * variant_allele_annotation(record, fields, header.sampleID, alt_idx)))
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
                    annotation = (isempty(id_field) ? allele_key : variant_token(id_field)) * "|" * allele_key
                    extra = variant_allele_annotation(record, fields, header.sampleID, alt_idx)
                    concrete_label = isempty(extra) ? annotation : annotation * extra
                    foreach_concrete_variant(alt) do concrete
                        label = ambiguous ? annotation * "|resolved=" * String(copy(concrete)) * extra : concrete_label
                        push!(variants, ScanVariant(ci, pos, first, stop, source, alt_idx, ref, alt,
                            concrete[(prefix + 1):(end - suffix)], allele_key, label))
                    end
                end
            end
        end
    finally
        close(stream)
        close(genome)
    end
    sv_skipped = Dict(k => n for (k, n) in counts if startswith(string(k), "sv_skipped_"))
    isempty(sv_skipped) || @warn "VCF structural variants skipped; only precise <DEL>, <DUP>, <DUP:TANDEM> and <INV> are searched" skipped=sv_skipped examples=sv_examples
    get(counts, :sv_dup_untyped, 0) > 0 && @warn "VCF <DUP> alleles searched as tandem duplications" alleles=counts[:sv_dup_untyped]
    counts[:iupac_skipped] > 0 && @warn "VCF ALT alleles with more than $PREFIX_HASH_SCAN_MAX_GUIDE_EXPANSIONS IUPAC realizations ignored" alleles=counts[:iupac_skipped] examples=iupac_examples
    sort!(variants; by = v -> (v.chrom, v.first, v.first == v.stop ? 0 : 1, v.source, v.annotation))
    return variants, counts, store
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
