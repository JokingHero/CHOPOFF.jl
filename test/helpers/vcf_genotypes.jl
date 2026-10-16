# Genotype-mode oracle. A sample is (; name, ploidy, calls), where calls maps an
# edit group to (alleles, label): alleles are ALT numbers of the record (0 REF,
# -1 missing) and label is :chrom (`|`, no PS), :unphased or a PS integer.
# Copies are built from complete haplotypes; production code is not called.

vtest_record_alts(edits, group) = unique(e.vcf_alt for e in edits if e.group == group)

function vtest_gt_field(alleles, label)
    gt = join((a < 0 ? "." : string(a) for a in alleles), label === :unphased ? '/' : '|')
    return gt * ":" * (label isa Integer ? string(label) : ".")
end

function vtest_genotype_vcf(path, edits, samples; chrom = "chr1")
    rows = Tuple[]
    for group in unique(e.group for e in edits)
        v = first(e for e in edits if e.group == group)
        fields = [haskey(s.calls, group) ? vtest_gt_field(s.calls[group]...) :
            vtest_gt_field(zeros(Int, s.ploidy), :chrom) for s in samples]
        push!(rows, (chrom, v.pos, v.id, v.ref, join(vtest_record_alts(edits, group), ','), ".", "PASS", ".", "GT:PS", fields...))
    end
    sort!(rows; by = r -> r[2], alg = MergeSort)
    return vtest_vcf(path, rows; header = ["##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">"], samples = [s.name for s in samples])
end

vtest_permutations(n) = n == 1 ? [[1]] : [[p[1:(k - 1)]; n; p[k:end]] for p in vtest_permutations(n - 1) for k in 1:n]

# truth[frame] lists, per assignment, the frame copies for which `holds(edits)`.
function vtest_sample_truth(edits, sample, holds)
    records = [(; group, alleles, label) for (group, (alleles, label)) in sample.calls if any(!=(0), alleles)]
    label(r) = r.label === :unphased || sample.ploidy == 1 ? (:u, r.group) : r.label
    labels = unique(label.(records))
    isempty(labels) && return Dict(:none => [[j for j in 1:sample.ploidy if holds(NamedTuple[])]])
    gaps = [(i, j) for (i, r) in enumerate(records) for j in 1:sample.ploidy if r.alleles[j] < 0]
    choices = [0:length(vtest_record_alts(edits, records[i].group)) for (i, _) in gaps]
    truth = Dict{Any, Vector{Vector{Int}}}()
    for (fi, frame) in enumerate(labels)
        truth[frame] = Vector{Int}[]
        perms = [k == fi ? [collect(1:sample.ploidy)] : vtest_permutations(sample.ploidy) for k in eachindex(labels)]
        for perm in Iterators.product(perms...), fill in Iterators.product(choices...)
            copies = Int[]
            for j in 1:sample.ploidy
                taken = NamedTuple[]
                for (i, r) in enumerate(records)
                    slot = perm[findfirst(==(label(r)), labels)][j]
                    g = findfirst(==((i, slot)), gaps)
                    a = g === nothing ? r.alleles[slot] : fill[g]
                    a > 0 && push!(taken, only(e for e in edits if e.group == r.group &&
                        e.vcf_alt == vtest_record_alts(edits, r.group)[a]))
                end
                holds(taken) && push!(copies, j)
            end
            push!(truth[frame], copies)
        end
    end
    return truth
end

# Checks one sample's tokens against the assignment truth table. For `lost`,
# a trailing `?` after certain tokens means more copies may be lost.
function vtest_check_tokens(tokens, truth, ploidy; lost = false)
    held = any(!isempty, h for t in values(truth) for h in t)
    isempty(tokens) && return !held
    held || return false
    always = all(!isempty, first(values(truth)))
    question = "?" in tokens
    certain = filter(!=("?"), tokens)
    lost || !question || length(tokens) == 1 || return false
    if isempty(certain)
        return ploidy > 2 || !always
    elseif all(==("*"), certain)
        k = length(certain)
        at_least = all(h -> length(h) >= k, first(values(truth)))
        lost || return at_least && !question
        return at_least && question == any(h -> length(h) > k, first(values(truth)))
    end
    slots = [parse(Int, first(split(t, '@'))) for t in certain]
    ps = occursin('@', certain[1]) ? parse(Int, split(certain[1], '@')[2]) : nothing
    frames = ps !== nothing ? [ps] : length(slots) < ploidy ? [:chrom] : collect(keys(truth))
    return any(frames) do f
        haskey(truth, f) && all(h -> all(in(h), slots), truth[f]) &&
            (!lost || question == any(h -> !issubset(h, slots), truth[f]))
    end
end

function vtest_sample_tokens(text, name)
    tokens = String[]
    for token in split(vtest_text(text), ';'; keepempty = false)
        sample, slot = rsplit(token, ':'; limit = 2)
        sample == vtest_escape(name) && push!(tokens, slot)
    end
    return tokens
end

vtest_certain(tokens) = count(!=("?"), tokens)

# The copy contains the row's site when both lifted endpoints exist and the
# edits it touches between them are exactly the row's alleles.
function vtest_copy_has_site(reference, taken, row, chrom; hap = vtest_haplotype(reference, taken; chrom))
    hap === nothing && error("generator placed conflicting edits on one copy")
    _, coordinates, events = hap
    i = findfirst(c -> vtest_origin(c) == (Int(row.start), vtest_text(row.start_origin)), coordinates)
    k = findfirst(c -> vtest_origin(c) == (Int(row.end), vtest_text(row.end_origin)), coordinates)
    (i === nothing || k === nothing) && return false
    lo, hi = minmax(i, k)
    alleles = [vtest_annotation(edit, chrom) for (a, b, edit) in events if (a > b ? lo < a <= hi : max(a, lo) <= min(b, hi))]
    return Tuple(sort(alleles)) == vtest_tokens(row.variants)
end

vtest_overlaps(v, lo, hi) = v.first == v.stop ? lo < v.first <= hi : v.first <= hi && v.stop > lo

function vtest_compare_genotypes(dir, reference, edits, samples, guides, motif, distance;
    chrom = "chr1", chunk = 2 * 1024 * 1024, threads = 1, case_id = "genotype", vcf = nothing)
    genome = vtest_fasta(joinpath(dir, "reference.fa"), [chrom => reference])
    vcf === nothing && (vcf = vtest_genotype_vcf(joinpath(dir, "input.vcf"), edits, samples; chrom))
    output = joinpath(dir, "genotype.csv")
    stops = fill(typemax(Int), distance + 1)
    CHOPOFF.search_prefixHashScan(guides, genome, motif, output; distance, vcf, vcf_mode = :genotype,
        stream_chunk_bases = chunk, scan_threads = threads, early_stopping = stops)
    actual = DataFrame(CSV.File(output; types = Dict(:carriers => String, :lost => String,
        :variants => String, :overlapping_variants => String, :start_origin => String, :end_origin => String)))
    union_output = joinpath(dir, "union.csv")
    CHOPOFF.search_prefixHashScan(guides, genome, motif, union_output; distance, vcf,
        stream_chunk_bases = chunk, scan_threads = threads, early_stopping = stops, vcf_max_alleles = vtest_exact_cap(edits))
    union_rows = Set(vtest_signature(row) for row in eachrow(DataFrame(CSV.File(union_output))))
    failures = String[]
    fail(message) = (push!(failures, message); nothing)
    variant_rows = filter(row -> !ismissing(row.variants), actual)
    reference_rows = filter(row -> ismissing(row.variants), actual)
    signatures = Dict(vtest_signature(row) => row for row in eachrow(variant_rows))
    haps = Dict{Vector{NamedTuple}, Any}()        # copies recur across rows
    length(signatures) == nrow(variant_rows) || fail("duplicate variant rows")
    for (signature, row) in signatures
        signature in union_rows || fail("genotype row absent from union output: $signature")
        certain = 0
        possible = 0
        for s in samples
            tokens = vtest_sample_tokens(row.carriers, s.name)
            truth = vtest_sample_truth(edits, s, taken -> vtest_copy_has_site(reference, taken, row, chrom;
                hap = get!(() -> vtest_haplotype(reference, taken; chrom), haps, taken)))
            vtest_check_tokens(tokens, truth, s.ploidy) || fail("carrier tokens $(s.name)=$tokens for $signature")
            certain += vtest_certain(tokens)
            possible += tokens == ["?"]
        end
        certain == row.n_carriers && possible == row.n_possible || fail("carrier counts for $signature")
        any(!isempty(vtest_sample_tokens(row.carriers, s.name)) for s in samples) || fail("row without carriers: $signature")
    end
    # Reference rows equal a search without variants; `lost` follows overlaps.
    # ALTs no sample can carry are dropped, so they overlap nothing.
    carried = filter(edits) do e
        a = findfirst(==(e.vcf_alt), vtest_record_alts(edits, e.group))
        any(s -> haskey(s.calls, e.group) && any(x -> x == a || x < 0, s.calls[e.group][1]), samples)
    end
    expected_reference = vtest_oracle(mkpath(joinpath(dir, "reference_oracle")), reference, carried, guides, motif, distance;
        chrom, max_alleles = 0)
    Set(vtest_signature(row) for row in eachrow(reference_rows)) == expected_reference || fail("reference rows differ")
    for row in eachrow(reference_rows)
        lo, hi = vtest_span(row, motif)
        certain = 0
        possible = 0
        for s in samples
            tokens = vtest_sample_tokens(row.lost, s.name)
            truth = vtest_sample_truth(edits, s, taken -> any(v -> vtest_overlaps(v, lo, hi), taken))
            vtest_check_tokens(tokens, truth, s.ploidy; lost = true) || fail("lost tokens $(s.name)=$tokens at $(row.start)$(row.strand)")
            certain += vtest_certain(tokens)
            possible += "?" in tokens
        end
        certain == row.n_lost && possible == row.n_lost_possible || fail("lost counts at $(row.start)$(row.strand)")
    end
    # Completeness: every variant site found on a possible copy is reported for its sample.
    copies = Dict{Vector{String}, Tuple{Vector{NamedTuple}, Set{String}}}()
    for s in samples
        vtest_sample_truth(edits, s, taken -> begin
            key = sort([e.id * e.vcf_alt for e in taken])
            entry = get!(() -> (taken, Set{String}()), copies, key)
            push!(entry[2], s.name)
            false
        end)
    end
    for (taken, names) in values(copies)
        isempty(taken) && continue
        for signature in vtest_copy_rows(mkpath(joinpath(dir, "copy_" * string(hash(taken)))), reference, taken,
                guides, motif, distance, chrom)
            row = get(signatures, signature, nothing)
            row === nothing && (fail("copy site missing: $signature for $(collect(names))"); continue)
            for name in names
                isempty(vtest_sample_tokens(row.carriers, name)) && fail("copy site lacks carrier $name: $signature")
            end
        end
    end
    if !isempty(failures)
        failure = mktempdir(; prefix = "chopoff_vcf_$(case_id)_", cleanup = false)
        for file in readdir(dir)
            isfile(joinpath(dir, file)) && cp(joinpath(dir, file), joinpath(failure, file))
        end
        write(joinpath(failure, "differences.txt"), join(failures, '\n'))
        @error "VCF genotype failure" case_id failure failures = length(failures) first = first(failures)
    end
    @test isempty(failures)
    return actual
end

# Variant sites of one complete copy, as vtest_oracle reports them.
function vtest_copy_rows(dir, reference, taken, guides, motif, distance, chrom)
    sequence, coordinates, events = vtest_haplotype(reference, taken; chrom)
    genome = vtest_fasta(joinpath(dir, "copy.fa"), ["h1" => sequence])
    output = vtest_search(dir, genome, [sequence], guides, motif, distance, joinpath(dir, "copy.csv"))
    rows = NamedTuple[]
    for row in eachrow(DataFrame(CSV.File(output)))
        lo, hi = vtest_span(row, motif)
        alleles = [vtest_annotation(edit, chrom) for (a, b, edit) in events if (a > b ? lo < a <= hi : max(a, lo) <= min(b, hi))]
        isempty(alleles) && continue
        start, start_origin = vtest_origin(coordinates[row.start])
        stop, end_origin = vtest_origin(coordinates[((row.strand == "+") == motif.extends5) ? lo : hi])
        push!(rows, (guide = string(row.guide), distance = Int(row.distance),
            chromosome = chrom, start, strand = string(row.strand), stop,
            variants = Tuple(sort(alleles)), overlapping_variants = (), start_origin, end_origin))
    end
    return rows
end

# Random calls; edits that conflict never both appear (or go missing) in one
# sample, since a copy cannot carry both.
function vtest_random_samples(rng, reference, edits; count = 3)
    groups = unique(e.group for e in edits)
    conflict(g, h) = any(vtest_haplotype(reference, [a, b]) === nothing
        for a in edits if a.group == g for b in edits if b.group == h)
    samples = []
    for k in 1:count
        ploidy = rand(rng, (1, 2, 2, 2, 3))
        calls = Dict{String, Any}()
        for g in groups
            any(h -> conflict(g, h), keys(calls)) && continue
            alts = length(vtest_record_alts(edits, g))
            alleles = [rand(rng) < 0.1 ? -1 : rand(rng, 0:alts) for _ in 1:ploidy]
            all(==(0), alleles) && continue
            calls[g] = (alleles, rand(rng, (:chrom, :chrom, :unphased, 100, 200)))
        end
        push!(samples, (; name = "s$k", ploidy, calls))
    end
    return samples
end
