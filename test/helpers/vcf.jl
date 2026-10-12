using Test, CHOPOFF, BioSequences, FASTX, TwoBit, CSV, DataFrames, CodecZlib

function vtest_fasta(path, sequences; width = 60)
    offset = 0
    open(path, "w") do out
        open(path * ".fai", "w") do fai
            for (name, sequence) in sequences
                write(out, ">", name, "\n")
                for start in 1:width:length(sequence)
                    println(out, sequence[start:min(end, start + width - 1)])
                end
                offset += sizeof(name) + 2
                println(fai, "$name\t$(length(sequence))\t$offset\t$(min(width, length(sequence)))\t$(min(width, length(sequence)) + 1)")
                offset += sizeof(sequence) + cld(length(sequence), width)
            end
        end
    end
    return path
end

function vtest_vcf(path, records; header = String[], samples = String[])
    open(path, "w") do io
        println(io, "##fileformat=VCFv4.2")
        foreach(line -> println(io, line), header)
        println(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            isempty(samples) ? "" : "\tFORMAT\t" * join(samples, '\t'))
        foreach(row -> println(io, join(row, '\t')), records)
    end
    return path
end

# Test edits state their footprints explicitly. The oracle does not call the
# production parser, conflict predicate, enumerator, geometry, or lifter.
function vtest_edit(sequence, first, stop, alt, id; original = alt, group = id)
    if first == stop || isempty(alt)
        if first == 1
            pos = 1
            ref = sequence[1:max(stop, 1)]
            vcf_alt = original * sequence[max(stop, 1):max(stop, 1)]
        else
            pos = first - 1
            ref = sequence[pos:max(pos, stop - 1)]
            vcf_alt = sequence[pos:pos] * original
        end
    else
        pos, ref, vcf_alt = first, sequence[first:(stop - 1)], original
    end
    concrete = first == stop || isempty(alt) ? (first == 1 ? alt * sequence[max(stop, 1):max(stop, 1)] : sequence[pos:pos] * alt) : alt
    return (; first, stop, alt, id, group, pos, ref, vcf_alt, concrete)
end

vtest_record(v; chrom = "chr1") = (chrom, v.pos, v.id, v.ref, v.vcf_alt, ".", "PASS", ".")

function vtest_haplotype(sequence, selected; chrom = "chr1")
    edits = sort(selected; by = v -> (v.first, v.first == v.stop ? 0 : 1))
    for i in eachindex(edits), j in 1:(i - 1)
        a, b = edits[i], edits[j]
        a.group == b.group && return nothing
        if a.first == a.stop && b.first == b.stop
            a.first == b.first && return nothing
        elseif a.first == a.stop
            b.first < a.first < b.stop && return nothing
        elseif b.first == b.stop
            a.first < b.first < a.stop && return nothing
        else
            !isempty(intersect(a.first:(a.stop - 1), b.first:(b.stop - 1))) && return nothing
        end
    end
    bases = Char[]
    coordinates = Tuple{Int, String}[]
    events = Tuple{Int, Int, NamedTuple}[]
    cursor = 1
    for v in edits
        for p in cursor:(v.first - 1)
            push!(bases, sequence[p])
            push!(coordinates, (p, ""))
        end
        start = length(bases) + 1
        for (offset, base) in enumerate(v.alt)
            push!(bases, base)
            origin = vtest_key(v, chrom) * ":$offset"
            push!(coordinates, length(v.alt) == v.stop - v.first ? (v.first + offset - 1, "") : (v.pos, origin))
        end
        push!(events, (start, length(bases), v))
        cursor = v.stop
    end
    for p in cursor:length(sequence)
        push!(bases, sequence[p])
        push!(coordinates, (p, ""))
    end
    return String(bases), coordinates, events
end

# Encoding expectations are independent of the production token formatter.
vtest_escape(s) = join(get(Dict('%' => "%25", ';' => "%3B", '|' => "%7C"), c, string(c)) for c in s)
vtest_key(v, chrom) = vtest_escape("$chrom:$(v.pos):$(v.ref)>$(v.vcf_alt)")
function vtest_annotation(v, chrom; resolved = true)
    id = v.id == "." ? "$chrom:$(v.pos):$(v.ref)>$(v.vcf_alt)" : v.id
    result = vtest_escape(id) * "|" * vtest_key(v, chrom)
    if resolved && any(c -> !(c in "ACGT"), v.vcf_alt)
        result *= "|resolved=" * v.concrete
    end
    return result
end

# Preserve original-record exclusivity, including concrete IUPAC realizations.
function vtest_records(edits; chrom = "chr1")
    rows = Tuple[]
    for group in unique([v.group for v in edits])
        choices = filter(v -> v.group == group, edits)
        v = first(choices)
        @assert all(x -> (x.pos, x.ref, x.id) == (v.pos, v.ref, v.id), choices)
        push!(rows, (chrom, v.pos, v.id, v.ref, join(unique(x.vcf_alt for x in choices), ','), ".", "PASS", "."))
    end
    return sort!(rows; by = r -> r[2], alg = MergeSort)
end

vtest_text(value) = ismissing(value) ? "" : string(value)
vtest_tokens(value) = isempty(vtest_text(value)) ? () : Tuple(sort(split(vtest_text(value), ';')))
function vtest_signature(row)
    return (guide = string(row.guide), alignment_guide = string(row.alignment_guide),
        alignment_reference = string(row.alignment_reference), distance = Int(row.distance),
        chromosome = string(row.chromosome), start = Int(row.start), strand = string(row.strand),
        stop = Int(row.end), variants = vtest_tokens(row.variants),
        overlapping_variants = vtest_tokens(row.overlapping_variants),
        start_origin = vtest_text(row.start_origin), end_origin = vtest_text(row.end_origin))
end

# Walk the declared motif positions in alignment order, then extend outward.
# This also handles an internal PAM; no production geometry helpers are used.
function vtest_span(row, motif)
    width = length(motif)
    right = (row.strand == "+") == motif.extends5
    left = right ? row.start - width + 1 : row.start
    pam = row.strand == "+" ? motif.pam_loci_fwd : motif.pam_loci_rve
    positions = [left + i - 1 for i in 1:width if !(i in pam)]
    right && reverse!(positions)
    step = right ? -1 : 1
    consumed = count(!=('-'), row.alignment_reference)
    while length(positions) < consumed
        push!(positions, last(positions) + step)
    end
    resize!(positions, consumed)
    append!(positions, (left + i - 1 for i in pam))
    return extrema(positions)
end

function vtest_haplotypes(reference, edits; chrom = "chr1", max_alleles = typemax(Int))
    groups = [[i for i in eachindex(edits) if edits[i].group == group]
        for group in unique([v.group for v in edits])]
    haplotypes = []
    for choices in Iterators.product(([0; group] for group in groups)...)
        selected = [edits[i] for i in choices if i != 0]
        length(selected) <= max_alleles || continue
        hap = vtest_haplotype(reference, selected; chrom)
        hap === nothing || push!(haplotypes, (hap..., !isempty(selected)))
    end
    return haplotypes
end

function vtest_oracle(dir, reference, edits, guides, motif, distance; chrom = "chr1", max_alleles = typemax(Int))
    expected = Set{NamedTuple}()
    haplotypes = vtest_haplotypes(reference, edits; chrom, max_alleles)
    genome = vtest_fasta(joinpath(dir, "oracle.fa"), ["h$i" => h[1] for (i, h) in enumerate(haplotypes)])
    output = joinpath(dir, "oracle.csv")
    CHOPOFF.search_prefixHashScan(guides, genome, motif, output;
        distance, query_variant = :bruteforce, early_stopping = fill(typemax(Int), distance + 1))
    for row in eachrow(DataFrame(CSV.File(output)))
        _, coordinates, events, alternative = haplotypes[parse(Int, row.chromosome[2:end])]
        lo, hi = vtest_span(row, motif)
        alleles = String[]
        for (a, b, edit) in events
            used = a > b ? lo < a <= hi : max(a, lo) <= min(b, hi)
            used && push!(alleles, vtest_annotation(edit, chrom))
        end
        isempty(alleles) && alternative && continue
        overlaps = String[]
        if !alternative
            for v in edits
                footprint = v.first == v.stop ? lo < v.first <= hi : !isempty(intersect(lo:hi, v.first:(v.stop - 1)))
                footprint && push!(overlaps, vtest_annotation(v, chrom; resolved = false))
            end
        end
        start, start_origin = coordinates[row.start]
        stop, end_origin = coordinates[((row.strand == "+") == motif.extends5) ? lo : hi]
        push!(expected, (guide = string(row.guide), alignment_guide = string(row.alignment_guide),
            alignment_reference = string(row.alignment_reference), distance = Int(row.distance),
            chromosome = chrom, start, strand = string(row.strand), stop,
            variants = Tuple(sort(alleles)), overlapping_variants = Tuple(sort(unique(overlaps))),
            start_origin, end_origin))
    end
    return expected
end

function vtest_assert_rows(actual, expected, dir; case_id = "vcf")
    observed = Set(vtest_signature(row) for row in eachrow(actual))
    missing_rows, extra_rows = setdiff(expected, observed), setdiff(observed, expected)
    if !isempty(missing_rows) || !isempty(extra_rows) || nrow(actual) != length(observed)
        # Keep all inputs outside mktempdir cleanup for a reproducible failure.
        failure = mktempdir(; prefix = "chopoff_vcf_$(case_id)_", cleanup = false)
        for file in readdir(dir)
            isfile(joinpath(dir, file)) && cp(joinpath(dir, file), joinpath(failure, file))
        end
        write(joinpath(failure, "differences.txt"), "Missing:\n$(repr(missing_rows))\nExtra:\n$(repr(extra_rows))\n")
        @error "VCF parity failure" case_id failure missing=length(missing_rows) extra=length(extra_rows)
    end
    @test isempty(missing_rows)
    @test isempty(extra_rows)
    @test nrow(actual) == length(observed)
    return actual
end

# Exact by default: the cap admits every record of the case at once.
vtest_exact_cap(edits) = max(CHOPOFF.PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES, length(unique(v.group for v in edits)))

function vtest_compare(dir, reference, edits, guides, motif, distance;
    chunk = 2 * 1024 * 1024, threads = 1, backend = :auto, chrom = "chr1", case_id = "vcf",
    genome = nothing, vcf = nothing, max_alleles = vtest_exact_cap(edits))
    genome === nothing && (genome = vtest_fasta(joinpath(dir, "reference.fa"), [chrom => reference]))
    vcf === nothing && (vcf = vtest_vcf(joinpath(dir, "input.vcf"), vtest_records(edits; chrom)))
    output = joinpath(dir, "actual.csv")
    CHOPOFF.search_prefixHashScan(guides, genome, motif, output;
        distance, vcf, stream_chunk_bases = chunk, scan_threads = threads, simd_backend = backend,
        early_stopping = fill(typemax(Int), distance + 1), vcf_max_alleles = max_alleles)
    actual = DataFrame(CSV.File(output))
    expected = vtest_oracle(dir, reference, edits, guides, motif, distance; chrom, max_alleles)
    return vtest_assert_rows(actual, expected, dir; case_id)
end

# Chromosomes are independent, so the oracle runs once per chromosome and the
# expected sets are merged. Production searches one genome and one VCF.
function vtest_compare_multi(dir, chromosomes, guides, motif, distance;
    chunk = 2 * 1024 * 1024, threads = 1, backend = :auto, case_id = "vcf_multi")
    genome = vtest_fasta(joinpath(dir, "reference.fa"), [chrom => reference for (chrom, (reference, _)) in chromosomes])
    records = reduce(vcat, (vtest_records(edits; chrom) for (chrom, (_, edits)) in chromosomes if !isempty(edits)); init = Tuple[])
    vcf = vtest_vcf(joinpath(dir, "input.vcf"), records)
    output = joinpath(dir, "actual.csv")
    cap = maximum(vtest_exact_cap(edits) for (_, (_, edits)) in chromosomes)
    CHOPOFF.search_prefixHashScan(guides, genome, motif, output;
        distance, vcf, stream_chunk_bases = chunk, scan_threads = threads, simd_backend = backend,
        early_stopping = fill(typemax(Int), distance + 1), vcf_max_alleles = cap)
    expected = Set{NamedTuple}()
    for (chrom, (reference, edits)) in chromosomes
        sub = mkpath(joinpath(dir, "oracle_" * chrom))
        union!(expected, vtest_oracle(sub, reference, edits, guides, motif, distance; chrom))
    end
    return vtest_assert_rows(DataFrame(CSV.File(output)), expected, dir; case_id)
end

# An edit whose VCF record carries padding that the reader must trim.
vtest_padded(sequence, pos, ref_len, vcf_alt, first, stop, alt, id; group = id) =
    (; first, stop, alt, id, group, pos, ref = sequence[pos:(pos + ref_len - 1)], vcf_alt, concrete = vcf_alt)
