# Realistic-VCF correctness checks for prefixHashScan (ClinVar, gnomAD).
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl prep    VCF [CHROM]
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl search  VCF CHROM OUTDIR
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl windows VCF CHROM N SEED [MAXGROUPS]
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl gtwindows VCF CHROM N SEED [SAMPLES]
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl gtsearch  VCF CHROM OUTDIR [S1,S2,...]
# `gtwindows` compares windows of real records and genotypes (SAMPLES random
# samples, default 20) with the genotype oracle. `gtsearch` runs vcf_mode=:genotype,
# checks reference rows and parity, and validates carrier tokens against GT/PS.
# `search` checks reference-row parity, run-to-run/thread/chunk/SIMD byte parity,
# and validates every variant row against a haplotype rebuilt by this script.
# `windows` compares random real-variant windows with the independent oracle.
using CHOPOFF, BioSequences, FASTX, CSV, DataFrames, Random, Test

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "test", "helpers", "vcf.jl"))
include(joinpath(ROOT, "test", "helpers", "vcf_sweeps.jl"))
include(joinpath(ROOT, "test", "helpers", "vcf_genotypes.jl"))
include(joinpath(ROOT, "test", "helpers", "vcf_real.jl"))
const GENOME = get(ENV, "CHOPOFF_PROFILE_GENOME",
    "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
const GUIDES = LongDNA{4}.(filter(!isempty, strip.(readlines(joinpath(@__DIR__, "data", "guides_for_tests.txt")))))
const K = parse(Int, get(ENV, "VCF_K", string(CHOPOFF.PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES)))

function chromosome_sequence(chrom)
    reader = open(FASTA.Reader, GENOME; index = GENOME * ".fai")
    try
        seekrecord(reader, chrom)
        return uppercase(sequence(String, read(reader)))
    finally
        close(reader)
    end
end

# Plain-text records via bcftools; independent of the production reader.
function read_records(vcf, chrom)
    records = NamedTuple[]
    cmd = chrom === nothing ? `bcftools view -H $vcf` : `bcftools view -H -t $chrom $vcf`
    for line in eachline(cmd)
        f = split(line, '\t'; limit = 9)
        pos = parse(Int, f[2])
        alts = String[]
        for (i, a) in enumerate(split(f[5], ','))
            a in (".", "*", "<NON_REF>", "<*>") && continue
            sv = symbolic_alt(a, f[8], pos, i)
            # Production ignores ALTs with more than 64 IUPAC realizations.
            sv === nothing ? (prod(length(IUPAC[c]) for c in uppercase(a)) <= 64 && push!(alts, uppercase(a))) :
                (isempty(sv) || push!(alts, sv))
        end
        push!(records, (chrom = String(f[1]), pos, id = String(f[3]),
            ref = uppercase(String(f[4])), alts))
    end
    return records
end

function prep(vcf, chrom)
    motif = Motif("Cas9"; distance = 3)
    dbi, lengths = CHOPOFF.prefix_hash_scan_dbinfo(GENOME, motif)
    selected, intervals = CHOPOFF.prefix_hash_scan_region(chrom, dbi.gi.chrom, lengths)
    t = @timed CHOPOFF.prepare_prefix_hash_scan_variants(vcf, dbi, lengths, selected, String[];
        verbose = true, intervals = intervals)
    data = t.value
    # Fragments exist one chromosome at a time, as in the variant scan.
    e = @timed map(ci -> begin
            bs = CHOPOFF.enumerate_variant_chromosome(data, dbi, lengths, ci, Threads.nthreads())
            (; fragments = sum(b -> length(b.fragments), bs), bases = sum(b -> sum(f -> f.bases, b.fragments; init = 0), bs),
                max_anchor = maximum(b -> maximum(b.per_anchor; init = 0), bs), bytes = sum(CHOPOFF.variant_descriptor_bytes, bs))
        end, findall(!isempty, data.anchors))
    println("prep vcf=$vcf read_seconds=$(round(t.time; digits = 2)) enumerate_seconds=$(round(e.time; digits = 2)) ",
        "allocated=$(t.bytes + e.bytes) counts=$(data.counts) variants=$(length(data.variants)) ",
        "fragments=$(sum(x -> x.fragments, e.value; init = 0)) fragment_bases=$(sum(x -> x.bases, e.value; init = 0)) ",
        "max_anchor_fragments=$(maximum(x -> x.max_anchor, e.value; init = 0)) ",
        "largest_chromosome_bytes=$(maximum(x -> x.bytes, e.value; init = 0)) peak_rss=$(Sys.maxrss())")
end

# --- row validator ---------------------------------------------------------

function parse_token(token)
    parts = split(token, '|')
    key = replace(parts[2], "%3B" => ";", "%7C" => "|", "%25" => "%")
    m = match(r"^(.*?):(\d+):([^:>]+)>(.+)$", key)
    resolved = nothing
    for p in parts[3:end]
        startswith(p, "resolved=") && (resolved = String(p[10:end]))
    end
    chrom, pos, ref, alt = m.captures
    return (key = String(parts[2]), chrom = String(chrom), pos = parse(Int, pos), ref = String(ref),
        alt = String(alt), concrete = something(resolved, String(alt)))
end

function check_row(row, seq, guide_set)
    alleles = parse_token.(split(row.variants, ';'))
    edits = map(alleles) do a
        first, stop, alt = footprint(a.ref, a.concrete, a.pos, seq)
        (; first, stop, alt, key = a.key, pos = a.pos)
    end
    lo_ref = minimum(e -> e.first, edits) - 120
    hi_ref = maximum(e -> e.stop, edits) + 120
    lo_ref = max(lo_ref, 1); hi_ref = min(hi_ref, length(seq))
    sort!(edits; by = e -> (e.first, e.first == e.stop ? 0 : 1))
    # Inserted bases keep (POS, allele key, offset); building "key:offset" per base
    # was slow for megabase duplication copies.
    bases = Char[]; coords = Tuple{Int, String, Int}[]; spans = UnitRange{Int}[]; junction = Int[]
    cursor = lo_ref
    for e in edits
        cursor <= e.first || return "conflicting alleles"
        for p in cursor:(e.first - 1)
            push!(bases, seq[p]); push!(coords, (p, "", 0))
        end
        start = length(bases) + 1
        # length(::String) is O(n); a megabase copy made this loop quadratic.
        replacement = ncodeunits(e.alt) == e.stop - e.first
        for (k, b) in enumerate(e.alt)
            push!(bases, b)
            push!(coords, replacement ? (e.first + k - 1, "", 0) : (e.pos, e.key, k))
        end
        push!(spans, start:length(bases))
        cursor = e.stop
    end
    for p in cursor:hi_ref
        push!(bases, seq[p]); push!(coords, (p, "", 0))
    end
    hap = String(bases)
    origin(o, p) = ismissing(o) ? (Int(p), "", 0) : (parts = rsplit(o, ':'; limit = 2); (Int(p), String(parts[1]), parse(Int, parts[2])))
    i = findfirst(==(origin(row.start_origin, row.start)), coords)
    j = findfirst(==(origin(row.end_origin, row.end)), coords)
    (i === nothing || j === nothing) && return "endpoint not in haplotype"
    lo, hi = minmax(i, j)
    site = hap[lo:hi]
    site = row.strand == "+" ? site : vtest_rc(site)
    # Cas9: protospacer then NGG, 5'->3' on the reported strand.
    endswith(site, "GG") || return "no NGG PAM"
    protospacer = site[1:(end - 3)]
    replace(row.alignment_reference, "-" => "") == protospacer || return "alignment_reference differs from haplotype"
    LongDNA{4}(replace(row.alignment_guide, "-" => "")) in guide_set || return "unknown guide"
    d = count(zip(row.alignment_guide, row.alignment_reference)) do (a, b)
        a == '-' || b == '-' || !iscompatible(DNA(a), DNA(b))
    end
    d == row.distance || return "distance $d != $(row.distance)"
    # Every annotated allele lies inside the consumed site or crosses into it.
    for s in spans
        touched = isempty(s) ? lo < first(s) <= hi : max(first(s), lo) <= min(last(s), hi)
        touched || return "annotated allele does not touch site"
    end
    return ""
end

function search(vcf, chrom_arg, outdir)
    mkpath(outdir)
    chrom = chrom_arg == "all" ? nothing : chrom_arg
    run_one(name; kw...) = begin
        path = joinpath(outdir, name * ".csv")
        t = @timed CHOPOFF.search_prefixHashScan(GUIDES, GENOME, Motif("Cas9"; distance = 3), path;
            distance = 3, region = chrom,
            early_stopping = fill(typemax(Int), 4), kw...)
        println("$name seconds=$(round(t.time; digits = 2))"); flush(stdout)
        path
    end
    base = run_one("reference")
    main = run_one("vcf"; vcf, vcf_max_alleles = K)
    variants = DataFrame(CSV.File(main; types = Dict(:start_origin => String, :end_origin => String,
        :variants => String, :overlapping_variants => String)))
    reference = DataFrame(CSV.File(base))
    refrows = select(variants[ismissing.(variants.variants), :], names(reference))
    @test sort(refrows) == sort(reference)
    for (name, kw) in (("threads1", (scan_threads = 1,)), ("chunk65536", (stream_chunk_bases = 65536,)),
            ("portable", (simd_backend = :portable,)))
        alt = run_one("vcf_" * name; vcf, vcf_max_alleles = K, kw...)
        @test read(alt) == read(main)
    end
    seqs = Dict{String, String}()
    guide_set = Set(GUIDES)
    failures = Dict{String, Int}()
    vrows = variants[.!ismissing.(variants.variants), :]
    for row in eachrow(vrows)
        seq = get!(() -> chromosome_sequence(string(row.chromosome)), seqs, string(row.chromosome))
        msg = check_row(row, seq, guide_set)
        isempty(msg) || (failures[msg] = get(failures, msg, 0) + 1)
    end
    println("variant rows=$(nrow(vrows)) reference rows=$(nrow(refrows)) validator failures=$failures")
    println("rows with start_origin=$(count(!ismissing, variants.start_origin)) end_origin=$(count(!ismissing, variants.end_origin))")
    @test isempty(failures)
    # Reference-row overlaps against an independent interval scan of the VCF.
    records = read_records(vcf, chrom)
    bychrom = Dict{String, Vector{NamedTuple}}()
    for r in records
        push!(get!(bychrom, r.chrom, NamedTuple[]), r)
    end
    maxref = maximum(r -> max(length(r.ref), maximum((footprint(r.ref, a, r.pos)[2] - r.pos for a in r.alts if !isempty(a)); init = 0)), records)
    positions = Dict(c => [r.pos for r in rs] for (c, rs) in bychrom)
    bad = 0
    for row in eachrow(variants[ismissing.(variants.variants), :])
        lo, hi = minmax(row.start, row.end)
        expected = Set{String}()
        rs = get(bychrom, string(row.chromosome), NamedTuple[])
        ps = get(positions, string(row.chromosome), Int[])
        for r in @view(rs[searchsortedfirst(ps, lo - maxref):searchsortedlast(ps, hi)]), a in r.alts
            first, stop, _ = footprint(r.ref, a, r.pos)
            hit = first == stop ? lo < first <= hi : first <= hi && stop > lo
            id = r.id == "." ? "$(r.chrom):$(r.pos):$(r.ref)>$a" : r.id
            hit && push!(expected, vtest_escape(id) * "|" * vtest_escape("$(r.chrom):$(r.pos):$(r.ref)>$a"))
        end
        observed = Set(vtest_tokens(row.overlapping_variants))
        observed == expected || (bad += 1)
    end
    println("reference rows with wrong overlapping_variants=$bad")
    @test bad == 0
end

# --- windowed oracle -------------------------------------------------------

function windows(vcf, chrom, n, seed, maxgroups)
    rng = MersenneTwister(seed)
    seq = chromosome_sequence(chrom)
    records = read_records(vcf, chrom)
    println("windows: $(length(records)) records on $chrom, seed=$seed, vcf_max_alleles=$K")
    positions = [r.pos for r in records]
    maxref = maximum(r -> length(r.ref), records)
    done = 0; attempts = 0; failed = String[]
    while done < n && attempts < 100 * n
        attempts += 1
        seed_record = records[rand(rng, eachindex(records))]
        pad = maxgroups >= 10 ? (35:50) : (60:160)          # dense VCFs need tight windows
        ws = seed_record.pos - rand(rng, pad)
        we = seed_record.pos + length(seed_record.ref) + rand(rng, pad)
        ws < 1 && continue
        w = seq[ws:we]
        all(c -> c in "ACGT", w) || continue
        near = @view records[searchsortedfirst(positions, ws - maxref):searchsortedlast(positions, we)]
        inside = filter(r -> ws + 30 <= r.pos && r.pos + length(r.ref) - 1 <= we - 30, near)
        partial = any(r -> r.pos <= we && r.pos + length(r.ref) - 1 >= ws &&
            !(ws + 30 <= r.pos && r.pos + length(r.ref) - 1 <= we - 30), near)
        partial && continue
        (isempty(inside) || length(inside) > maxgroups) && continue
        any(r -> isempty(r.alts), inside) && continue
        haplotype_count([vcat(concrete_alleles.(r.alts)...) for r in inside]) <= (maxgroups >= 10 ? 1024 : 256) || continue
        edits = window_edits(seq, inside, ws, we)
        d = rand(rng, 1:3)
        guides = guides_for(rng, w, edits, d)
        isempty(guides) && continue
        done += 1
        dir = mktempdir(; cleanup = true)
        chunk = rand(rng, (64, 97, 2 * 1024 * 1024))
        try
            @testset "window $done $chrom:$ws-$we records=$(length(inside)) d=$d chunk=$chunk" begin
                vtest_compare(dir, w, edits, LongDNA{4}.(guides), Motif("Cas9"; distance = d), d;
                    chrom = "w", chunk, threads = min(4, Threads.nthreads()), case_id = "real_$(chrom)_$ws",
                max_alleles = K)
            end
        catch err
            err isa Test.TestSetException || rethrow()
            push!(failed, "$chrom:$ws-$we")
        end
    end
    println("windows done=$done attempts=$attempts failed=$(length(failed)) $failed")
    isempty(failed) || exit(1)
end


# --- genotype mode -------------------------------------------------------------

# Records with the selected samples' GT/PS: calls[sample] = (alleles, label) with
# ALT numbers of the original record, label :chrom, :unphased or a PS integer.
function read_gt_records(vcf, chrom, samples)
    names = samples === nothing ? split(readchomp(`bcftools query -l $vcf`), '\n') : samples
    list = join(names, ',')
    records = NamedTuple[]
    for line in eachline(`bcftools view -H -t $chrom -s $list $vcf`)
        f = split(line, '\t')
        keys = split(f[9], ':')
        gi, pi = findfirst(==("GT"), keys), findfirst(==("PS"), keys)
        calls = Dict{String, Any}()
        for (k, name) in enumerate(names)
            v = split(f[9 + k], ':')
            gt = v[gi]
            alleles = [a == "." ? -1 : parse(Int, a) for a in split(gt, r"[|/]") if !isempty(a)]
            ps = pi === nothing || length(v) < pi ? "." : v[pi]
            label = occursin('/', gt) && length(alleles) > 1 ? :unphased : ps == "." ? :chrom : parse(Int, ps)
            calls[String(name)] = (alleles, label)
        end
        pos = parse(Int, f[2])
        push!(records, (chrom = String(f[1]), pos, id = String(f[3]), ref = uppercase(String(f[4])),
            alts = [searched_alt(a, f[8], pos, i) for (i, a) in enumerate(split(f[5], ','))], calls))
    end
    return String.(names), records
end

function gtwindows(vcf, chrom, n, seed, nsamples)
    rng = MersenneTwister(seed)
    seq = chromosome_sequence(chrom)
    all_names = split(readchomp(`bcftools query -l $vcf`), '\n')
    names = String.(length(all_names) <= nsamples ? all_names : shuffle(rng, all_names)[1:nsamples])
    _, records = read_gt_records(vcf, chrom, names)
    println("gtwindows: $(length(records)) records on $chrom, samples=$(length(names)), seed=$seed")
    positions = [r.pos for r in records]
    maxref = maximum(r -> length(r.ref), records)
    done = 0; attempts = 0; failed = String[]
    while done < n && attempts < 200 * n
        attempts += 1
        seed_record = records[rand(rng, eachindex(records))]
        ws = seed_record.pos - rand(rng, 60:160)
        we = seed_record.pos + length(seed_record.ref) + rand(rng, 60:160)
        ws < 1 && continue
        w = seq[ws:we]
        all(c -> c in "ACGT", w) || continue
        near = @view records[searchsortedfirst(positions, ws - maxref):searchsortedlast(positions, we)]
        inside = filter(r -> ws + 30 <= r.pos && r.pos + length(r.ref) - 1 <= we - 30, near)
        any(r -> r.pos <= we && r.pos + length(r.ref) - 1 >= ws &&
            !(ws + 30 <= r.pos && r.pos + length(r.ref) - 1 <= we - 30), near) && continue
        (isempty(inside) || length(inside) > 12) && continue
        # The oracle needs literal ACGT ALTs, numbered as in the file.
        all(r -> all(a -> all(in("ACGT"), a), r.alts), inside) || continue
        edits = window_edits(seq, [(; r.chrom, r.pos, r.id, r.ref, r.alts) for r in inside], ws, we)
        samples = [(; name, ploidy = length(first(r.calls[name][1] for r in inside)),
            calls = Dict{String, Any}("g$k" => r.calls[name] for (k, r) in enumerate(inside)
                if any(!=(0), r.calls[name][1]))) for name in names]
        any(s -> !isempty(s.calls), samples) || continue
        any(s -> any(r -> length(r.calls[s.name][1]) != s.ploidy, inside), samples) && continue
        # A copy cannot carry conflicting edits; skip windows where one sample may.
        conflict(g, h) = any(vtest_haplotype(w, [a, b]) === nothing for a in edits if a.group == g for b in edits if b.group == h)
        any(s -> any(conflict(g, h) for g in keys(s.calls) for h in keys(s.calls) if g < h), samples) && continue
        haps = Dict{Vector{String}, Any}()
        for s in samples
            vtest_sample_truth(edits, s, taken -> (haps[sort([e.group * e.vcf_alt for e in taken])] = taken; false))
        end
        d = rand(rng, 1:3)
        guides = String[]
        for taken in shuffle(rng, collect(values(haps)))
            length(guides) >= 3 && break
            isempty(taken) && continue
            hap, _, events = vtest_haplotype(w, taken; chrom = "w")
            sites = [(p, rev) for p in 1:(length(hap) - 22) for rev in (false, true)
                if (rev ? hap[p] == 'C' && hap[p + 1] == 'C' : hap[p + 21] == 'G' && hap[p + 22] == 'G') &&
                    any(ev -> max(ev[1], p) <= min(ev[2], p + 22) || (ev[1] > ev[2] && p < ev[1] <= p + 22), events)]
            isempty(sites) && continue
            p, rev = rand(rng, sites)
            g = collect(rev ? vtest_rc(hap[p:(p + 22)])[1:20] : hap[p:(p + 19)])
            for _ in 1:rand(rng, 0:d)
                g[rand(rng, 1:20)] = rand(rng, collect("ACGT"))
            end
            all(in("ACGT"), g) && push!(guides, String(g))
        end
        isempty(guides) && continue
        done += 1
        dir = mktempdir(; cleanup = true)
        try
            @testset "gt window $done $chrom:$ws-$we records=$(length(inside)) d=$d" begin
                vtest_compare_genotypes(dir, w, edits, samples, LongDNA{4}.(unique(guides)), Motif("Cas9"; distance = d), d;
                    chrom = "w", chunk = rand(rng, (64, 97, 2 * 1024 * 1024)), threads = min(4, Threads.nthreads()),
                    case_id = "realgt_$(chrom)_$ws")
            end
        catch err
            err isa Test.TestSetException || rethrow()
            push!(failed, "$chrom:$ws-$we")
        end
    end
    println("gtwindows done=$done attempts=$attempts failed=$(length(failed)) $failed")
    isempty(failed) || exit(1)
end

vtest_unescape(x) = replace(x, "%3B" => ";", "%7C" => "|", "%25" => "%")

# Checks carrier tokens against the GT/PS calls of the alleles a row uses:
# `S:j[@PS]` needs allele on copy j of that phase set for every row allele and
# REF on copy j for other phased calls overlapping the site span; other tokens
# need every row allele in some copy or missing.
function check_carriers(row, by_pos, records, positions, maxref)
    tokens = [split(t, '|') for t in vtest_tokens(row.variants)]
    # GATK-SV writes several records at one POS, all with REF N.
    alleles = map(tokens) do t
        m = match(r"^(.*?):(\d+):([^:>]+)>(.+)$", vtest_unescape(t[2]))
        pos = parse(Int, m[2])
        id = vtest_unescape(t[1]) == vtest_unescape(t[2]) ? "." : vtest_unescape(t[1])
        r = only(r for r in by_pos[pos] if r.ref == m[3] && r.id == id && m[4] in r.alts)
        (r, findfirst(==(m[4]), r.alts))
    end
    lo, hi = minmax(row.start, row.end)
    bad = String[]
    for token in split(vtest_text(row.carriers), ';'; keepempty = false)
        name, slot = rsplit(token, ':'; limit = 2)
        name = vtest_unescape(name)
        if slot in ("?", "*")
            all(((r, a),) -> any(x -> x == a || x < 0, r.calls[name][1]), alleles) || push!(bad, "possible carrier lacks allele")
            continue
        end
        j, ps = occursin('@', slot) ? (parse(Int, split(slot, '@')[1]), parse(Int, split(slot, '@')[2])) :
            (parse(Int, slot), nothing)
        for (r, a) in alleles
            gt, label = r.calls[name]
            homo = all(==(a), gt)
            homo || (gt[j] == a && label == something(ps, :chrom)) || push!(bad, "copy lacks row allele")
        end
        for r in @view records[searchsortedfirst(positions, lo - maxref):searchsortedlast(positions, hi)]
            any(((u, _),) -> u === r, alleles) && continue
            gt, label = r.calls[name]
            # A bare "." call has one element and no copy j.
            label == something(ps, :chrom) && j <= length(gt) && gt[j] > 0 || continue
            isempty(r.alts[gt[j]]) && continue   # skipped SV kinds are not searched
            first_, stop_, _ = footprint(r.ref, r.alts[gt[j]], r.pos)
            inside_deletion = any(((u, b),) -> (f = footprint(u.ref, u.alts[b], u.pos); f[1] <= first_ && stop_ <= f[2]), alleles)
            hit = first_ == stop_ ? lo < first_ <= hi : first_ <= hi && stop_ > lo
            hit && !inside_deletion && push!(bad, "copy carries other overlapping allele")
        end
    end
    return bad
end

function gtsearch(vcf, chrom, outdir, samples)
    mkpath(outdir)
    run_one(name; kw...) = begin
        path = joinpath(outdir, name * ".csv")
        t = @timed CHOPOFF.search_prefixHashScan(GUIDES, GENOME, Motif("Cas9"; distance = 3), path;
            distance = 3, region = chrom, early_stopping = fill(typemax(Int), 4), kw...)
        println("$name seconds=$(round(t.time; digits = 2)) peak_rss=$(Sys.maxrss())"); flush(stdout)
        path
    end
    base = run_one("reference")
    main = run_one("genotype"; vcf, vcf_mode = :genotype, vcf_samples = samples)
    types = Dict(:start_origin => String, :end_origin => String, :variants => String,
        :overlapping_variants => String, :carriers => String, :lost => String)
    rows = DataFrame(CSV.File(main; types))
    reference = DataFrame(CSV.File(base))
    refrows = select(rows[ismissing.(rows.variants), :], names(reference))
    @test sort(refrows) == sort(reference)
    for (name, kw) in (("threads1", (scan_threads = 1,)), ("chunk65536", (stream_chunk_bases = 65536,)))
        alt = run_one("genotype_" * name; vcf, vcf_mode = :genotype, vcf_samples = samples, kw...)
        @test read(alt) == read(main)
    end
    _, records = read_gt_records(vcf, chrom, samples)
    by_pos = Dict{Int, Vector{NamedTuple}}()
    foreach(r -> push!(get!(by_pos, r.pos, NamedTuple[]), r), records)
    positions = [r.pos for r in records]
    maxref = maximum(r -> max(length(r.ref), maximum((footprint(r.ref, a, r.pos)[2] - r.pos for a in r.alts if !isempty(a)); init = 0)), records)
    vrows = rows[.!ismissing.(rows.variants), :]
    failures = Dict{String, Int}()
    seq = chromosome_sequence(chrom)
    guide_set = Set(GUIDES)
    for row in eachrow(vrows)
        for msg in [check_row(row, seq, guide_set); check_carriers(row, by_pos, records, positions, maxref)]
            isempty(msg) || (failures[msg] = get(failures, msg, 0) + 1)
        end
    end
    certain = sum(vrows.n_carriers; init = 0)
    println("variant rows=$(nrow(vrows)) certain carriers=$certain possible samples=$(sum(vrows.n_possible; init = 0)) " *
        "reference rows=$(nrow(refrows)) lost=$(sum(rows.n_lost[ismissing.(rows.variants)]; init = 0)) failures=$failures")
    @test isempty(failures)
end

mode = ARGS[1]
if mode == "prep"
    prep(ARGS[2], length(ARGS) >= 3 ? ARGS[3] : nothing)
elseif mode == "search"
    @testset "realistic search" begin
        search(ARGS[2], ARGS[3], ARGS[4])
    end
elseif mode == "windows"
    windows(ARGS[2], ARGS[3], parse(Int, ARGS[4]), parse(Int, ARGS[5]), length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 8)
elseif mode == "gtwindows"
    gtwindows(ARGS[2], ARGS[3], parse(Int, ARGS[4]), parse(Int, ARGS[5]), length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 20)
elseif mode == "gtsearch"
    @testset "realistic genotype search" begin
        gtsearch(ARGS[2], ARGS[3], ARGS[4], length(ARGS) >= 5 ? String.(split(ARGS[5], ',')) : nothing)
    end
else
    error("unknown mode $mode")
end
