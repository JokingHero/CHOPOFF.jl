# Sample genotypes for vcf_mode=:genotype. Entries are sparse: a sample has an
# entry at a record only when some allele is ALT or missing; hom-ref is implicit.
#
# Phase groups: 0 is the chromosome-wide set of `|` calls without PS, PS + 1 a
# PS block, and -source an unphased record, whose slots permute freely.
const VARIANT_GT_MISSING = 0xff

struct ScanGenotypes
    samples::Vector{String}
    record_source::Vector{Int}      # sorted sources of records with entries
    record_first::Vector{Int32}     # entry range of record i: record_first[i]:(record_first[i+1]-1)
    sample::Vector{Int32}
    group::Vector{Int64}
    allele_first::Vector{Int32}     # alleles of entry e: allele_first[e]:(allele_first[e+1]-1)
    alleles::Vector{UInt8}
end

ScanGenotypes(samples) = ScanGenotypes(samples, Int[], Int32[1], Int32[], Int64[], Int32[1], UInt8[])

# Parses one GT value. VCF 4.4 allows a leading phase indicator; without it the
# first allele is phased unless another indicator is `/`.
function parse_variant_gt!(alleles, gt::AbstractString)
    empty!(alleles)
    phased = true
    i = firstindex(gt)
    last = lastindex(gt)
    i <= last && (gt[i] == '|' || gt[i] == '/') && (phased &= gt[i] == '|'; i = nextind(gt, i))
    while i <= last
        j = i
        while j <= last && gt[j] != '|' && gt[j] != '/'
            j = nextind(gt, j)
        end
        token = SubString(gt, i, prevind(gt, j))
        if token == "."
            push!(alleles, VARIANT_GT_MISSING)
        else
            value = tryparse(Int, token)
            value !== nothing && 0 <= value < VARIANT_GT_MISSING || error("Malformed VCF GT: $gt")
            push!(alleles, UInt8(value))
        end
        j <= last || break
        phased &= gt[j] == '|'
        i = nextind(gt, j)
        i <= last || error("Malformed VCF GT: $gt")
    end
    isempty(alleles) && error("Malformed VCF GT: $gt")
    return phased
end

function variant_sample_columns(sample_names, requested)
    isempty(sample_names) && error("vcf_mode=:genotype needs sample columns with GT; use vcf_mode=:union for sites-only VCFs")
    requested === nothing && return collect(eachindex(sample_names))
    columns = Int[]
    for name in requested
        i = findfirst(==(name), sample_names)
        i === nothing && error("vcf_samples names a sample absent from the VCF: $name")
        push!(columns, i)
    end
    return sort!(unique!(columns))
end

# Appends the selected samples' genotypes of one record. `fields` is the text
# after ALT; returns the possibly carried alleles as a bit set (bit a for ALT a).
function push_variant_genotypes!(store, alleles, source, fields::AbstractString, columns, where)
    n = ncodeunits(fields)
    bytes = codeunits(fields)
    # A byte loop is ~20x faster than findnext('\t') on 3,200-sample lines.
    function tab(i)
        while i <= n && @inbounds(bytes[i]) != UInt8('\t')
            i += 1
        end
        return min(i, n + 1)
    end
    t = 0
    for _ in 1:3                    # QUAL, FILTER, INFO
        t = tab(t + 1)
    end
    t > ncodeunits(fields) && error("VCF record has no FORMAT column at $where")
    format_end = tab(t + 1)
    keys = split(SubString(fields, t + 1, format_end - 1), ':')
    gt_key = findfirst(==("GT"), keys)
    ps_key = findfirst(==("PS"), keys)
    gt_key === nothing && error("VCF record has no GT at $where; vcf_mode=:genotype needs GT")
    carried = UInt128(0)
    column = 1
    start = format_end + 1
    entry_count = 0
    for c in columns
        while column < c
            start = tab(start) + 1
            column += 1
        end
        start > ncodeunits(fields) && error("VCF record has fewer sample columns than the header at $where")
        stop = tab(start) - 1
        # Cohort calls are mostly diploid hom-ref; skip them without parsing.
        if gt_key == 1 && stop - start >= 2 && bytes[start] == UInt8('0') && bytes[start + 2] == UInt8('0') &&
                (bytes[start + 1] == UInt8('|') || bytes[start + 1] == UInt8('/')) &&
                (stop - start == 2 || bytes[start + 3] == UInt8(':'))
            continue
        end
        value = SubString(fields, start, stop)
        gt = variant_subfield(value, gt_key)
        gt === nothing && error("VCF sample value lacks GT at $where")
        phased = parse_variant_gt!(alleles, gt)
        all(iszero, alleles) && continue
        ps = ps_key === nothing ? nothing : variant_subfield(value, ps_key)
        ps === nothing && (ps = ".")
        group = !phased || length(alleles) == 1 ? -Int64(source) :
            ps == "." ? Int64(0) : something(tryparse(Int64, ps), -1) + 1
        group == 0 && ps != "." && error("Malformed VCF PS at $where: $ps")
        push!(store.sample, Int32(c))
        push!(store.group, group)
        append!(store.alleles, alleles)
        push!(store.allele_first, Int32(length(store.alleles) + 1))
        entry_count += 1
        for a in alleles
            carried |= a == VARIANT_GT_MISSING ? typemax(UInt128) : UInt128(1) << min(a, 127)
        end
    end
    if entry_count > 0
        push!(store.record_source, source)
        push!(store.record_first, Int32(length(store.sample) + 1))
    end
    return carried
end

# Parses the genotypes of a batch of record lines on all threads; returns each
# line's carried-ALT bits. Chunk stores are appended in line order, so entries
# match a sequential parse. Lines off `selected` chromosomes or too short to
# hold FORMAT are left to the caller, which skips or rejects them.
function push_variant_genotype_batch!(store, lines, sources, chrom_index, selected, columns)
    carried = zeros(UInt128, length(lines))
    tasks = map(Iterators.partition(eachindex(lines), cld(length(lines), Threads.nthreads()))) do part
        Threads.@spawn begin
            chunk = ScanGenotypes(store.samples)
            alleles = UInt8[]
            for i in part
                text = lines[i]
                bytes = codeunits(text)
                tabs = Int[]
                while length(tabs) < 5
                    t = findnext(==(UInt8('\t')), bytes, isempty(tabs) ? 1 : tabs[end] + 1)
                    t === nothing && break
                    push!(tabs, t)
                end
                length(tabs) >= 4 || continue
                ci = get(chrom_index, SubString(text, 1, tabs[1] - 1), 0)
                ci != 0 && selected[ci] || continue
                fields = SubString(text, length(tabs) == 5 ? tabs[5] + 1 : ncodeunits(text) + 1)
                where = SubString(text, 1, tabs[1] - 1) * ":" * SubString(text, tabs[1] + 1, tabs[2] - 1)
                carried[i] = push_variant_genotypes!(chunk, alleles, sources[i], fields, columns, where)
            end
            chunk
        end
    end
    chunks = try
        fetch.(tasks)
    catch e
        e isa TaskFailedException ? throw(e.task.exception) : rethrow()
    end
    for chunk in chunks
        entries, allele_count = Int32(length(store.sample)), Int32(length(store.alleles))
        append!(store.record_source, chunk.record_source)
        append!(store.record_first, @view(chunk.record_first[2:end]) .+ entries)
        append!(store.sample, chunk.sample)
        append!(store.group, chunk.group)
        append!(store.allele_first, @view(chunk.allele_first[2:end]) .+ allele_count)
        append!(store.alleles, chunk.alleles)
    end
    return carried
end

# Field k of a colon-separated sample value, or nothing when absent.
function variant_subfield(value::AbstractString, k)
    i = firstindex(value)
    for _ in 2:k
        c = findnext(':', value, i)
        c === nothing && return nothing
        i = c + 1
    end
    c = findnext(':', value, i)
    return SubString(value, i, c === nothing ? lastindex(value) : c - 1)
end

function variant_record_entries(store, source)
    r = searchsortedfirst(store.record_source, source)
    (r > length(store.record_source) || store.record_source[r] != source) && return 1:0
    return Int(store.record_first[r]):(Int(store.record_first[r + 1]) - 1)
end

variant_entry_alleles(store, e) = @view store.alleles[store.allele_first[e]:(store.allele_first[e + 1] - 1)]

# Requirement of one row on one record: `required` is the taken ALT index (0 if
# none) and `forbidden` the ALT indexes the copy must not carry.
struct VariantGenotypeRequirement
    source::Int
    required::UInt8
    forbidden::Vector{UInt8}
end

function variant_requirements(taken, excluded)
    by_source = Dict{Int, VariantGenotypeRequirement}()
    for (source, a) in taken
        by_source[source] = VariantGenotypeRequirement(source, UInt8(a), UInt8[])
    end
    for (source, a) in excluded
        r = get!(() -> VariantGenotypeRequirement(source, 0x00, UInt8[]), by_source, source)
        UInt8(a) == r.required || UInt8(a) in r.forbidden || push!(r.forbidden, UInt8(a))
    end
    return sort!(collect(values(by_source)); by = r -> r.source)
end

# Per phase group, slot masks of copies that certainly and possibly satisfy (or,
# for `lost`, certainly and possibly violate) the requirements.
struct VariantGroupMasks
    group::Int64
    ploidy::Int
    certain::UInt64
    possible::UInt64
end

full_slots(ploidy) = ploidy >= 64 ? typemax(UInt64) : (UInt64(1) << ploidy) - 1

function variant_group_masks!(groups, store, sample, entries, requirements, lost)
    empty!(groups)
    for (r, range) in zip(requirements, entries)
        e = searchsortedfirst(@view(store.sample[range]), sample)
        e = e <= length(range) && store.sample[range[e]] == sample ? range[e] : 0
        if e == 0                    # hom-ref
            r.required == 0 && continue
            lost && continue
            return false
        end
        alleles = variant_entry_alleles(store, e)
        length(alleles) <= 64 || error("VCF ploidy above 64 is unsupported")
        certain = UInt64(0)
        possible = UInt64(0)
        for (j, a) in enumerate(alleles)
            bit = UInt64(1) << (j - 1)
            if lost                  # copy carries a forbidden allele
                a != VARIANT_GT_MISSING && a in r.forbidden && (certain |= bit)
                (a == VARIANT_GT_MISSING || a in r.forbidden) && (possible |= bit)
            else
                ok = (r.required == 0 || a == r.required) && !(a in r.forbidden)
                ok && a != VARIANT_GT_MISSING && (certain |= bit)
                (ok || a == VARIANT_GT_MISSING) && (possible |= bit)
            end
        end
        group = store.group[e]
        i = findfirst(g -> g.group == group, groups)
        if i === nothing
            push!(groups, VariantGroupMasks(group, length(alleles), certain, possible))
        elseif lost
            g = groups[i]
            groups[i] = VariantGroupMasks(group, max(g.ploidy, length(alleles)), g.certain | certain, g.possible | possible)
        else
            g = groups[i]
            groups[i] = VariantGroupMasks(group, min(g.ploidy, length(alleles)), g.certain & certain, g.possible & possible)
        end
    end
    return true
end

# One `S:*` per copy that holds in every assignment without a known index.
variant_unknown_tokens!(io, name, mask) = foreach(_ -> print(io, ';', name, ":*"), 1:count_ones(mask))

function variant_slot_tokens!(io, name, mask, group)
    for j in 1:64
        mask & (UInt64(1) << (j - 1)) == 0 && continue
        print(io, ';', name, ':', j)
        group > 0 && print(io, '@', group - 1)
    end
end

# Tokens for the copies that carry a row (lost=false: every taken allele and no
# excluded allele) or lack a reference site (lost=true: some excluded allele).
# `S:j` holds in every assignment of permutable phase groups and missing
# alleles, `S:*` some copy holds but which is unknown, `S:?` some assignment
# holds. Returns (tokens, certain copies, samples with an `S:?` token). For
# lost=true, `S:?` may follow certain tokens when further copies may be lost.
function haplotype_tokens(store::ScanGenotypes, taken, excluded; lost = false)
    requirements = variant_requirements(taken, excluded)
    isempty(requirements) && return "", 0, 0
    entries = [variant_record_entries(store, r.source) for r in requirements]
    candidates = Int32[]
    if lost
        for range in entries
            append!(candidates, @view store.sample[range])
        end
        unique!(sort!(candidates))
    else
        required = [range for (r, range) in zip(requirements, entries) if r.required != 0]
        isempty(required) && error("haplotype_tokens needs a taken allele unless lost=true")
        candidates = store.sample[required[1]]
        for range in @view required[2:end]
            candidates = intersect(candidates, @view store.sample[range])
        end
    end
    io = IOBuffer()
    certain_count = 0
    possible_count = 0
    groups = VariantGroupMasks[]
    for sample in candidates
        variant_group_masks!(groups, store, sample, entries, requirements, lost) || continue
        name = variant_token(store.samples[sample])
        if lost
            any(g -> g.possible != 0, groups) || continue
            full = findfirst(g -> g.certain == full_slots(g.ploidy), groups)
            hit = findall(g -> g.certain != 0, groups)
            # A trailing `S:?` marks that more copies than the certain tokens
            # name may be lost.
            if full !== nothing
                variant_slot_tokens!(io, name, full_slots(groups[full].ploidy), 0)
                certain_count += groups[full].ploidy
                continue
            elseif length(hit) == 1
                g = groups[hit[1]]
                g.group >= 0 ? variant_slot_tokens!(io, name, g.certain, g.group) : variant_unknown_tokens!(io, name, g.certain)
                certain_count += count_ones(g.certain)
                more = count_ones(g.certain) < g.ploidy && (g.possible & ~g.certain != 0 ||
                    any(h -> h.group != g.group && h.possible != 0, groups))
            elseif !isempty(hit)
                print(io, ';', name, ":*")
                certain_count += 1
                more = true
            else
                more = true
            end
            if more
                print(io, ';', name, ":?")
                possible_count += 1
            end
        else
            all(g -> g.possible != 0, groups) || continue
            partial = findall(g -> g.certain != full_slots(g.ploidy), groups)
            if length(partial) == 0
                ploidy = minimum(g -> g.ploidy, groups)
                variant_slot_tokens!(io, name, full_slots(ploidy), 0)
                certain_count += ploidy
            elseif length(partial) == 1 && groups[partial[1]].certain != 0
                g = groups[partial[1]]
                if g.group >= 0
                    variant_slot_tokens!(io, name, g.certain, g.group)
                    certain_count += count_ones(g.certain)
                else
                    variant_unknown_tokens!(io, name, g.certain)
                    certain_count += count_ones(g.certain)
                end
            else
                print(io, ';', name, ":?")
                possible_count += 1
            end
        end
    end
    tokens = String(take!(io))
    return isempty(tokens) ? "" : tokens[2:end], certain_count, possible_count
end

variant_max_ploidy(store::ScanGenotypes) = maximum(diff(store.allele_first); init = 1)

# Slots (sample - 1) * ploidy + copy that may carry an ALT, for pruning walks.
# Only chromosome-wide phase keeps copy identity across records; any other
# carrier marks every copy of its sample.
function variant_carrier_slots(store::ScanGenotypes, source, alt_index, ploidy)
    slots = Int32[]
    for e in variant_record_entries(store, source)
        base = (Int(store.sample[e]) - 1) * ploidy
        alleles = variant_entry_alleles(store, e)
        if store.group[e] == 0
            for (j, a) in enumerate(alleles)
                (a == alt_index || a == VARIANT_GT_MISSING) && push!(slots, base + j)
            end
        elseif any(a -> a == alt_index || a == VARIANT_GT_MISSING, alleles)
            append!(slots, (base + 1):(base + ploidy))
        end
    end
    return slots
end

function variant_live_slots(a::Vector{Int32}, b::Vector{Int32})
    out = Int32[]
    i, j = 1, 1
    while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            push!(out, a[i])
            i += 1
            j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    return out
end

# Some call sets phase overlapping alleles onto one copy (or make both
# homozygous) instead of using `*`. Walks never take both, and the overlap rule
# credits such a copy only with rows of the allele that hides the other.
function check_variant_genotype_conflicts(variants, by_chrom, store, chrom_names)
    checked = Set{NTuple{4, Int}}()
    conflicts = 0
    examples = String[]
    for range in by_chrom, i in range
        a = variants[i]
        reach = max(a.stop, a.first + 1)
        for j in (i + 1):last(range)
            b = variants[j]
            b.first < reach || break
            a.source == b.source && continue
            scan_variants_conflict(a, b) || continue
            key = (a.source, a.alt_index, b.source, b.alt_index)
            key in checked && continue
            push!(checked, key)
            tokens, certain, _ = haplotype_tokens(store, [(a.source, a.alt_index), (b.source, b.alt_index)], Tuple{Int, Int}[])
            certain > 0 || continue
            conflicts += 1
            length(examples) < 5 && push!(examples, "$(chrom_names[a.chrom]):$(a.pos)+$(b.pos) ($(first(split(tokens, ';'))))")
        end
    end
    conflicts > 0 && @warn "VCF genotypes place overlapping alleles on one copy; a row using the hidden allele does not credit that copy" pairs=conflicts examples
    return conflicts
end

# Selected sample names (VCF order) and their largest GT ploidy per chromosome,
# over records on `chroms`. Hom-ref samples never appear in output tokens.
function variant_sample_ploidy(path, chroms, samples)
    stream = open_variant_vcf(path)
    names = String[]
    ploidy = Dict{Tuple{String, String}, Int}()
    alleles = UInt8[]
    columns = Int[]
    try
        for line in eachline(stream)
            if startswith(line, "##")
                continue
            elseif startswith(line, '#')
                header = split(line, '\t')
                names = String.(header[10:end])
                columns = variant_sample_columns(names, samples)
                continue
            end
            tab = findfirst('\t', line)
            chrom = SubString(line, 1, tab - 1)
            chrom in chroms || continue
            fields = split(line, '\t')
            keys = split(fields[9], ':')
            gt = findfirst(==("GT"), keys)
            gt === nothing && error("VCF record has no GT; per-sample summaries need GT")
            for c in columns
                parse_variant_gt!(alleles, variant_subfield(fields[9 + c], gt))
                key = (names[c], String(chrom))
                ploidy[key] = max(get(ploidy, key, 0), length(alleles))
            end
        end
    finally
        close(stream)
    end
    return names[columns], ploidy
end
