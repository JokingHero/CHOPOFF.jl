# Realistic-VCF correctness checks for prefixHashScan (ClinVar, gnomAD).
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl prep    VCF [CHROM]
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl search  VCF CHROM OUTDIR
#   julia --project=. -t 8 test/local_human/vcf_realistic.jl windows VCF CHROM N SEED [MAXGROUPS]
# `search` checks reference-row parity, run-to-run/thread/chunk/SIMD byte parity,
# and validates every variant row against a haplotype rebuilt by this script.
# `windows` compares random real-variant windows with the independent oracle.
using CHOPOFF, BioSequences, FASTX, CSV, DataFrames, Random, Test

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "test", "helpers", "vcf.jl"))
include(joinpath(ROOT, "test", "helpers", "vcf_sweeps.jl"))
const GENOME = get(ENV, "CHOPOFF_PROFILE_GENOME",
    "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
const GUIDES = LongDNA{4}.(filter(!isempty, strip.(readlines(joinpath(@__DIR__, "data", "guides_for_tests.txt")))))
const K = parse(Int, get(ENV, "VCF_K", string(CHOPOFF.PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES)))
const IUPAC = Dict('A'=>"A",'C'=>"C",'G'=>"G",'T'=>"T",'R'=>"AG",'Y'=>"CT",'S'=>"CG",'W'=>"AT",'K'=>"GT",
    'M'=>"AC",'B'=>"CGT",'D'=>"AGT",'H'=>"ACT",'V'=>"ACG",'N'=>"ACGT")

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
        f = split(line, '\t'; limit = 6)
        alts = [String(a) for a in split(f[5], ',') if !(a in (".", "*", "<NON_REF>", "<*>"))]
        # Production ignores ALTs with more than 64 IUPAC realizations.
        filter!(a -> prod(length(IUPAC[c]) for c in uppercase(a)) <= 64, alts)
        push!(records, (chrom = String(f[1]), pos = parse(Int, f[2]), id = String(f[3]),
            ref = uppercase(String(f[4])), alts = uppercase.(alts)))
    end
    return records
end

function trim_footprint(ref, alt, pos)
    p = 0
    while p < min(length(ref), length(alt)) && ref[p + 1] == alt[p + 1] && ref[p + 1] in "ACGT"
        p += 1
    end
    s = 0
    while s < min(length(ref), length(alt)) - p && ref[end - s] == alt[end - s] && ref[end - s] in "ACGT"
        s += 1
    end
    p == length(ref) == length(alt) && ((p, s) = (0, 0))
    return pos + p, pos + length(ref) - s, alt[(p + 1):(end - s)]
end

function concrete_alleles(alt)
    out = [""]
    for c in alt
        out = [o * b for o in out for b in IUPAC[c]]
    end
    return out
end

function prep(vcf, chrom)
    motif = Motif("Cas9"; distance = 3)
    dbi, lengths = CHOPOFF.prefix_hash_scan_dbinfo(GENOME, motif)
    selected = CHOPOFF.prefix_hash_scan_chromosomes(chrom, dbi.gi.chrom)
    t = @timed CHOPOFF.prepare_prefix_hash_scan_variants(vcf, dbi, lengths, selected, String[]; verbose = true)
    data = t.value
    println("prep vcf=$vcf seconds=$(round(t.time; digits = 2)) allocated=$(t.bytes) counts=$(data.counts) ",
        "variants=$(length(data.variants)) fragments=$(length(data.fragments)) ",
        "fragment_bases=$(sum(f -> f.bases, data.fragments; init = 0)) max_anchor_fragments=$(data.max_anchor_fragments) ",
        "peak_rss=$(Sys.maxrss())")
end

# --- row validator ---------------------------------------------------------

function parse_token(token)
    parts = split(token, '|')
    key = replace(parts[2], "%3B" => ";", "%7C" => "|", "%25" => "%")
    m = match(r"^(.*):(\d+):([^:>]+)>([^:]+)$", key)
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
        first, stop, alt = trim_footprint(a.ref, a.concrete, a.pos)
        (; first, stop, alt, key = a.key, pos = a.pos)
    end
    lo_ref = minimum(e -> e.first, edits) - 120
    hi_ref = maximum(e -> e.stop, edits) + 120
    lo_ref = max(lo_ref, 1); hi_ref = min(hi_ref, length(seq))
    sort!(edits; by = e -> (e.first, e.first == e.stop ? 0 : 1))
    bases = Char[]; coords = Tuple{Int, String}[]; spans = UnitRange{Int}[]; junction = Int[]
    cursor = lo_ref
    for e in edits
        cursor <= e.first || return "conflicting alleles"
        for p in cursor:(e.first - 1)
            push!(bases, seq[p]); push!(coords, (p, ""))
        end
        start = length(bases) + 1
        for (k, b) in enumerate(e.alt)
            push!(bases, b)
            push!(coords, length(e.alt) == e.stop - e.first ? (e.first + k - 1, "") : (e.pos, e.key * ":$k"))
        end
        push!(spans, start:length(bases))
        cursor = e.stop
    end
    for p in cursor:hi_ref
        push!(bases, seq[p]); push!(coords, (p, ""))
    end
    hap = String(bases)
    so = ismissing(row.start_origin) ? "" : row.start_origin
    eo = ismissing(row.end_origin) ? "" : row.end_origin
    i = findfirst(==((Int(row.start), so)), coords)
    j = findfirst(==((Int(row.end), eo)), coords)
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
    maxref = maximum(r -> length(r.ref), records)
    positions = Dict(c => [r.pos for r in rs] for (c, rs) in bychrom)
    bad = 0
    for row in eachrow(variants[ismissing.(variants.variants), :])
        lo, hi = minmax(row.start, row.end)
        expected = Set{String}()
        rs = get(bychrom, string(row.chromosome), NamedTuple[])
        ps = get(positions, string(row.chromosome), Int[])
        for r in @view(rs[searchsortedfirst(ps, lo - maxref):searchsortedlast(ps, hi)]), a in r.alts
            first, stop, _ = trim_footprint(r.ref, a, r.pos)
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

function haplotype_count(groups)
    prod(1 + length(g) for g in groups; init = 1)
end

function window_edits(seq, records, ws, we)
    edits = NamedTuple[]
    for (k, r) in enumerate(records)
        for a in r.alts, c in concrete_alleles(a)
            first, stop, alt = trim_footprint(r.ref, c, r.pos)
            push!(edits, (; first = first - ws + 1, stop = stop - ws + 1, alt, id = r.id, group = "g$k",
                pos = r.pos - ws + 1, ref = r.ref, vcf_alt = a, concrete = c))
        end
    end
    return edits
end

function guides_for(rng, w, edits, d)
    haps = vtest_haplotypes(w, edits; chrom = "w")
    guides = String[]
    for (hap, _, events, alternative) in shuffle(rng, haps)
        length(guides) >= 3 && break
        sites = Tuple{Int, Bool}[]
        for p in 1:(length(hap) - 22)
            hap[p + 21] == 'G' && hap[p + 22] == 'G' && push!(sites, (p, false))
            hap[p] == 'C' && hap[p + 1] == 'C' && push!(sites, (p, true))
        end
        near = filter(sites) do (p, rev)
            any(ev -> max(ev[1], p) <= min(ev[2], p + 22) || (ev[1] > ev[2] && p < ev[1] <= p + 22), events) ||
                !alternative
        end
        isempty(near) && continue
        p, rev = rand(rng, near)
        g = collect(rev ? vtest_rc(hap[p:(p + 22)])[1:20] : hap[p:(p + 19)])
        for _ in 1:rand(rng, 0:d)
            i = rand(rng, 1:20); g[i] = rand(rng, collect("ACGT"))
        end
        any(c -> !(c in "ACGT"), g) && continue
        push!(guides, String(g))
    end
    return unique(guides)
end

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

mode = ARGS[1]
if mode == "prep"
    prep(ARGS[2], length(ARGS) >= 3 ? ARGS[3] : nothing)
elseif mode == "search"
    @testset "realistic search" begin
        search(ARGS[2], ARGS[3], ARGS[4])
    end
elseif mode == "windows"
    windows(ARGS[2], ARGS[3], parse(Int, ARGS[4]), parse(Int, ARGS[5]), length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 8)
else
    error("unknown mode $mode")
end
