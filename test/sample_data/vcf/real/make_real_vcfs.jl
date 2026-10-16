# Cuts small fixtures from the large real call sets (GRCh38) used to calibrate
# vcf.md. Each fixture contig is one GRCh38 window; records keep REF, ALT, ID and
# GT/PS and shift POS and END. Windows are chosen greedily so that every record
# and neighbourhood tag (variant kind, size, GT pattern, overlaps, clusters, SV
# breakpoints near small variants) appears at least twice, and each window stays
# small enough for the haplotype oracle. Needs bcftools and the local sources.
#   julia --project=. test/sample_data/vcf/real/make_real_vcfs.jl            # write fixtures
#   julia --project=. test/sample_data/vcf/real/make_real_vcfs.jl --check    # regenerate and compare bytes
# Expected outputs: CHOPOFF_VCF_REAL_UPDATE=1 julia --project=. -e 'include("test/src/variants.jl")'
using CHOPOFF, BioSequences, FASTX, CSV, DataFrames, Random, Test
const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
include(joinpath(ROOT, "test", "helpers", "vcf.jl"))
include(joinpath(ROOT, "test", "helpers", "vcf_sweeps.jl"))
include(joinpath(ROOT, "test", "helpers", "vcf_genotypes.jl"))
include(joinpath(ROOT, "test", "helpers", "vcf_real.jl"))

const GENOME = get(ENV, "CHOPOFF_PROFILE_GENOME",
    "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
const LOCAL = joinpath(ROOT, "test", "local_human", "data", "vcf")
const SV = normpath(joinpath(ROOT, "..", "sv_data"))
const TRIO = ["NA12878", "NA12891", "NA12892"]
# region: bcftools -t targets; streamed, so unindexed files work.
const SOURCES = [
    (set = "hg002", vcf = joinpath(SV, "HG002_GIAB_v5.0q_small_and_SV_literal_GRCh38.normalized.vcf.gz"),
        samples = ["HG002"], region = "22,X", cap = 0),
    (set = "hg001", vcf = joinpath(LOCAL, "giab_hg001_norm.vcf.gz"), samples = ["HG001"], region = "22", cap = 0),
    (set = "kg_panel", vcf = joinpath(LOCAL, "kg_chr22_norm.vcf.gz"), samples = TRIO, region = "22", cap = 0),
    (set = "kg_ps", vcf = joinpath(LOCAL, "hg001_kg_chr22_ps.vcf.gz"), samples = ["NA12878"], region = "22", cap = 0),
    (set = "kg_sv", vcf = joinpath(SV, "1kGP_3202_highcov_SV_GRCh38.renamed.vcf.gz"), samples = TRIO, region = nothing, cap = 0),
    (set = "clinvar", vcf = joinpath(LOCAL, "clinvar_norm.vcf.gz"), samples = String[], region = "22", cap = 0),
    (set = "gnomad", vcf = joinpath(LOCAL, "gnomad_chr22_pass.vcf.gz"), samples = String[],
        region = "22:30000000-33000000", cap = 3),
]
const KEEP_INFO = ("END", "SVLEN", "SVTYPE", "IMPRECISE", "AF", "CLNSIG")
const SEED = 20261015
const DISTANCE = 3
const PAD = 40          # flank around the records of a window
const GAP = 30          # records closer than this join one window
const MAX_SPAN = 12_000
const MAX_WINDOWS = 24
const MAX_HAPLOTYPES = 1024
const MAX_ORACLE_BASES = 100_000   # haplotypes x haplotype length; the brute-force oracle is slow
const MAX_ASSIGNMENTS = 64      # genotype oracle phase/missing-GT assignments per sample

function chromosome(reader, name)
    seekrecord(reader, name)
    return uppercase(sequence(String, read(reader)))
end

skipped(r, i) = isempty(searched_alt(r.alts[i], r.info, r.pos, i))

function info_value(info, key)
    for item in split(info, ';')
        startswith(item, key * "=") && return item[(length(key) + 2):end]
    end
    return nothing
end

# Last reference base touched by a searched allele; skipped kinds count as their POS.
function record_stop(pos, ref, alts, info)
    stop = pos + length(ref) - 1
    for (i, a) in enumerate(alts)
        startswith(a, '<') || continue
        sv = symbolic_alt(a, info, pos, i)
        isempty(something(sv, "")) || (stop = max(stop, parse(Int, split(sv, '/')[end])))
    end
    return stop
end

function kind_tags(ref, a, info)
    if startswith(a, '<')
        t = uppercase(a[2:(end - 1)])
        base = startswith(t, "DEL:ME") ? "DEL:ME" : startswith(t, "INS:ME") ? "INS:ME" : t
        tags = ["sv:" * base]
        occursin("IMPRECISE", info) && push!(tags, "sv_imprecise")
        info_value(info, "END") === nothing && push!(tags, "sv_svlen_only")
        return tags
    end
    a in ("*", ".") && return ["star"]
    any(c -> !(c in "ACGT"), a) && return [occursin('N', a) ? "alt_N" : "alt_iupac"]
    lr, la = length(ref), length(a)
    kind = lr == la == 1 ? "snv" : lr == la ? "mnp" :
        la > lr ? (lr == 1 && ref[1] == a[1] ? "ins" : "replace") : (la == 1 && ref[1] == a[1] ? "del" : "replace")
    size = abs(la - lr)
    tags = [kind]
    size >= 50 && push!(tags, kind * "50")
    size >= 1000 && push!(tags, kind * "1k")
    size >= 10_000 && push!(tags, kind * "10k")
    return tags
end

# One source record: shifted fields are rebuilt on output.
function parse_line(line, nsamples)
    f = split(line, '\t')
    pos = parse(Int, f[2])
    alts = String.(split(f[5], ','))
    info = join(filter(x -> first(split(x, '=')) in KEEP_INFO, split(f[8], ';')), ';')
    gts = String[]
    if nsamples > 0
        keys = split(f[9], ':')
        gi, pi = findfirst(==("GT"), keys), findfirst(==("PS"), keys)
        for v in f[10:end]
            x = split(v, ':')
            push!(gts, pi === nothing ? x[gi] : x[gi] * ":" * (length(x) >= pi ? x[pi] : "."))
        end
    end
    ref = uppercase(String(f[4]))
    tags = Set{String}()
    length(alts) > 1 && push!(tags, "multiallelic")
    kinds = Set{String}()
    for a in alts
        t = kind_tags(ref, a, String(f[8]))
        union!(tags, t); push!(kinds, first(t))
    end
    length(kinds) > 1 && push!(tags, "multiallelic_mixed")
    f[3] == "." && push!(tags, "missing_id")
    for g in gts
        gt = first(split(g, ':'))
        push!(tags, "gt:" * replace(gt, r"[3-9]" => "n"))
        occursin(':', g) && !endswith(g, ":.") && push!(tags, "ps")
    end
    # Literal insertions that repeat the following reference become <DUP:TANDEM> in the twin.
    return (; chrom = String(f[1]), pos, id = String(f[3]), ref, alts, info, filter = String(f[7]),
        gts, tags, stop = record_stop(pos, ref, alts, String(f[8])), symbolic = any(a -> startswith(a, '<'), alts),
        long = any(a -> !startswith(a, '<') && abs(length(a) - length(ref)) >= 50, alts))
end

function read_source(src)
    sample_args = isempty(src.samples) ? `` : `-s $(join(src.samples, ','))`
    region_args = src.region === nothing ? `` : `-t $(src.region)`
    header = String[]
    names = String[]
    bychrom = Dict{String, Vector{NamedTuple}}()
    for line in eachline(`bcftools view $region_args $sample_args $(src.vcf)`)
        if startswith(line, "##")
            push!(header, line)
        elseif startswith(line, '#')
            names = String.(split(line, '\t')[10:end])
        else
            r = parse_line(line, length(names))
            push!(get!(bychrom, r.chrom, NamedTuple[]), r)
        end
    end
    return header, names, bychrom
end

is_big(r) = r.symbolic || r.long

# Tags that need several records of one window.
function window_tags(rs, seq)
    tags = Set{String}()
    foreach(r -> union!(tags, r.tags), rs)
    for i in eachindex(rs), j in (i + 1):lastindex(rs)
        a, b = rs[i], rs[j]
        a.pos == b.pos && push!(tags, "same_pos")
        if b.pos <= a.stop
            push!(tags, is_big(a) && is_big(b) ? "sv_sv_overlap" : is_big(a) || is_big(b) ? "sv_overlap_small" : "overlap")
        elseif b.pos == a.stop + 1
            push!(tags, "adjacent")
        elseif b.pos - a.stop <= 23 && (is_big(a) || is_big(b))
            push!(tags, "sv_near_small")
        end
    end
    for i in eachindex(rs)
        count(r -> rs[i].pos <= r.pos <= rs[i].pos + 30, rs) >= 5 && push!(tags, "cluster5")
        r = rs[i]
        dels = any(a -> !startswith(a, '<') && length(a) < length(r.ref), r.alts) || any(a -> startswith(uppercase(a), "<DEL"), r.alts)
        dels && any(x -> x !== r && r.pos - 23 <= x.stop < r.pos, rs) && any(x -> x !== r && r.stop < x.pos <= r.stop + 23, rs) &&
            push!(tags, "del_bridge")
        # Tandem copies: a long insertion equal to the following reference.
        for a in r.alts
            n = length(a) - 1
            !startswith(a, '<') && length(r.ref) == 1 && n >= 50 && a[1] == r.ref[1] &&
                a[2:end] == seq[(r.pos + 1):min(end, r.pos + n)] && push!(tags, "ins50_tandem")
        end
    end
    if !isempty(rs) && !isempty(first(rs).gts)
        for k in eachindex(first(rs).gts)
            ps = unique(split(r.gts[k], ':')[end] for r in rs if occursin(':', r.gts[k]) && occursin('|', r.gts[k]))
            length(ps) > 1 && push!(tags, "ps_switch")
            phased = unique(occursin('/', r.gts[k]) for r in rs if count(c -> c in "|/", r.gts[k]) > 0)
            length(phased) > 1 && push!(tags, "phase_mix")
        end
        any(r -> length(unique(r.gts)) > 1, rs) && length(first(rs).gts) > 1 && push!(tags, "samples_differ")
    end
    return tags
end

# Literal records of a window, coordinates relative to `ws`.
function literal_records(rs, seq, ws)
    map(rs) do r
        pos = r.pos - ws + 1
        # END is absolute; the searched form "<TYPE>/END" is shifted with POS.
        alts = map(enumerate(r.alts)) do (i, a)
            x = searched_alt(a, r.info, r.pos, i)
            startswith(x, '<') ? x[1:findlast('/', x)] * string(parse(Int, x[(findlast('/', x) + 1):end]) - ws + 1) : x
        end
        if any(a -> startswith(a, '<'), alts)
            ref, alt = vreal_literal(seq, pos, only(a for a in alts if startswith(a, '<')))
            alts = [startswith(a, '<') ? alt : a for a in alts]
        else
            ref = r.ref
        end
        (; chrom = "w", pos, id = r.id, ref, alts)
    end
end

function hap_bound(sizes, cap)
    cap == 0 && return prod(1 + s for s in sizes; init = 1)
    e = zeros(Int, cap + 1); e[1] = 1             # elementary symmetric sums up to cap
    for s in sizes, j in (cap + 1):-1:2
        e[j] += e[j - 1] * s
    end
    return sum(e)
end

# Samples of a window for the genotype oracle, or `nothing` when a sample may
# carry conflicting alleles on one copy (the oracle cannot build that copy).
function window_samples(rs, edits, names, seq_w)
    isempty(names) && return nothing
    groups = Set(e.group for e in edits)
    samples = []
    for (k, name) in enumerate(names)
        calls = Dict{String, Any}()
        ploidy = 0
        for (i, r) in enumerate(rs)
            g = "g$i"
            alleles, label = vreal_call("GT:PS", occursin(':', r.gts[k]) ? r.gts[k] : r.gts[k] * ":.")
            ploidy = max(ploidy, length(alleles))
            any(!=(0), alleles) || continue
            g in groups || continue          # no searched ALT: skipped SV kinds
            # The oracle numbers searched ALTs only.
            any(i -> skipped(r, i), eachindex(r.alts)) && return nothing
            calls[g] = (alleles, label)
        end
        any(r -> length(vreal_call("GT:PS", occursin(':', r.gts[k]) ? r.gts[k] : r.gts[k] * ":.")[1]) != ploidy, rs) && return nothing
        conflict(g, h) = any(vtest_conflict(a, b) for a in edits if a.group == g for b in edits if b.group == h)
        any(conflict(g, h) for g in keys(calls) for h in keys(calls) if g < h) && return nothing
        # Mirrors vtest_sample_truth: one frame per label, copy orders of the other labels, missing-GT fills.
        labels = unique(l === :unphased || ploidy == 1 ? (:u, g) : l for (g, (_, l)) in calls)
        fills = prod(1 + length(vtest_record_alts(edits, g)) for (g, (a, _)) in calls for x in a if x < 0; init = 1)
        length(labels) * factorial(ploidy)^max(0, length(labels) - 1) * fills <= MAX_ASSIGNMENTS || return nothing
        push!(samples, (; name, ploidy, calls))
    end
    return samples
end

# Candidate window around records lo:hi of one chromosome, or `nothing`.
function candidate(rs, prefix_stop, lo, hi, seq, src)
    span_lo, span_hi = rs[lo].pos, maximum(r -> r.stop, @view rs[lo:hi])
    lo > 1 && prefix_stop[lo - 1] >= span_lo && return nothing
    ws = max(1, span_lo - PAD, lo > 1 ? prefix_stop[lo - 1] + 1 : 1)
    we = min(ncodeunits(seq), span_hi + PAD, hi < length(rs) ? rs[hi + 1].pos - 1 : ncodeunits(seq))
    (we - ws > MAX_SPAN || we - ws < 50 || we < span_hi) && return nothing
    w = seq[ws:we]
    all(c -> c in "ACGT", w) || return nothing
    # A window without a reference PAM site may yield no guides.
    any(p -> w[(p + 21):(p + 22)] == "GG" || w[p:(p + 1)] == "CC", 1:(length(w) - 22)) || return nothing
    window = rs[lo:hi]
    any(r -> r.symbolic && length(r.alts) > 1, window) && return nothing
    lit = literal_records(window, w, ws)
    # A literal INV with complementary ends would be trimmed to a smaller footprint.
    any(window) do r
        any(a -> uppercase(a) == "<INV>", r.alts) || return false
        p, stop = r.pos - ws + 1, r.stop - ws + 1
        w[p + 1] == only(vtest_rc(w[stop:stop]))
    end && return nothing
    # Concrete alleles of an IUPAC ALT, counted without expanding (saturates at 65).
    n_concrete(a) = foldl((n, c) -> min(n * length(IUPAC[c]), 65), a; init = 1)
    any(r -> any(a -> n_concrete(a) > 64, r.alts), lit) && return nothing
    haps = hap_bound([sum(a -> isempty(a) ? 0 : n_concrete(a), r.alts) for r in lit], src.cap)
    bases = length(w) + sum(r -> maximum(length, r.alts; init = 0), lit; init = 0)
    haps <= MAX_HAPLOTYPES && haps * bases <= MAX_ORACLE_BASES || return nothing
    return (; chrom = rs[lo].chrom, ws, we, records = window, tags = window_tags(window, seq))
end

function candidates(bychrom, reader, src)
    out = NamedTuple[]
    for chrom in sort(collect(keys(bychrom)))
        rs = sort(bychrom[chrom]; by = r -> r.pos, alg = MergeSort)
        seq = chromosome(reader, chrom)
        prefix_stop = accumulate(max, [r.stop for r in rs])
        seen = Set{Tuple{Int, Int}}()
        for i in eachindex(rs)
            lo = hi = i
            # Join neighbours within GAP while the window stays oracle-sized.
            while true
                grown = false
                if hi < length(rs) && rs[hi + 1].pos <= maximum(r -> r.stop, @view rs[lo:hi]) + GAP && hi - lo < 11
                    c = candidate(rs, prefix_stop, lo, hi + 1, seq, src)
                    c === nothing || (hi += 1; grown = true)
                end
                if !grown && lo > 1 && prefix_stop[lo - 1] >= rs[lo].pos - GAP && hi - lo < 11
                    c = candidate(rs, prefix_stop, lo - 1, hi, seq, src)
                    c === nothing || (lo -= 1; grown = true)
                end
                grown || break
            end
            (lo, hi) in seen && continue
            push!(seen, (lo, hi))
            c = candidate(rs, prefix_stop, lo, hi, seq, src)
            c === nothing || push!(out, c)
        end
    end
    return out
end

function select_windows(cands)
    covered = Dict{String, Int}()
    chosen = NamedTuple[]
    overlaps(a, b) = a.chrom == b.chrom && a.ws <= b.we && b.ws <= a.we
    score(c) = count(t -> get(covered, t, 0) < 2, c.tags)
    while length(chosen) < MAX_WINDOWS
        pool = filter(c -> !any(x -> overlaps(c, x), chosen), cands)
        isempty(pool) && break
        best = pool[argmin(map(c -> (-score(c), c.we - c.ws, c.chrom, c.ws), pool))]
        score(best) == 0 && break
        push!(chosen, best)
        foreach(t -> covered[t] = get(covered, t, 0) + 1, best.tags)
    end
    return sort(chosen; by = c -> (c.chrom, c.ws))
end

contig_name(set, c) = "$(set)_$(c.chrom)_$(c.ws)"

function shifted_line(r, name, ws, symbolic_twin, seq_w)
    pos = r.pos - ws + 1
    ref, alts, info = r.ref, copy(r.alts), r.info
    if symbolic_twin && length(alts) == 1 && !startswith(alts[1], '<') && ref[1] == alts[1][1] &&
            (length(ref) == 1 || length(alts[1]) == 1)
        n = length(alts[1]) - 1
        if length(alts[1]) == 1 && length(ref) - 1 >= 50
            ref, alts, info = ref[1:1], ["<DEL>"], "END=$(r.pos + length(r.ref) - 1)" * (isempty(info) ? "" : ";" * info)
        elseif length(ref) == 1 && n >= 50 && alts[1][2:end] == seq_w[(pos + 1):min(end, pos + n)]
            ref, alts, info = ref, ["<DUP:TANDEM>"], "END=$(r.pos + n)" * (isempty(info) ? "" : ";" * info)
        end
    end
    end_value = info_value(info, "END")
    end_value === nothing || (info = replace(info, "END=$end_value" => "END=$(parse(Int, end_value) - ws + 1)"))
    fields = [name, string(pos), r.id, ref, join(alts, ','), ".", r.filter, isempty(info) ? "." : info]
    isempty(r.gts) || append!(fields, [any(g -> occursin(':', g), r.gts) ? "GT:PS" : "GT"; r.gts])
    return join(fields, '\t')
end

function write_set_vcf(path, header, names, chosen, set, contigs, symbolic_twin, source)
    infos = Set(k for c in chosen for r in c.records for k in first.(split.(split(r.info, ';'), '=')) if !isempty(k))
    symbolic_twin && push!(infos, "END")
    has_ps = any(g -> occursin(':', g), (g for c in chosen for r in c.records for g in r.gts))
    open(path, "w") do io
        println(io, "##fileformat=VCFv4.2")
        println(io, "##source=make_real_vcfs.jl from $(basename(source)); contigs are GRCh38 windows, see windows.tsv")
        for c in chosen
            println(io, "##contig=<ID=$(contig_name(set, c)),length=$(c.we - c.ws + 1)>")
        end
        for key in sort(collect(infos))
            line = findfirst(l -> startswith(l, "##INFO=<ID=$key,"), header)
            println(io, line === nothing ? "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">" : header[line])
        end
        if !isempty(names)
            println(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
            has_ps && println(io, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">")
        end
        println(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", isempty(names) ? "" : "\tFORMAT\t" * join(names, '\t'))
        for c in chosen
            for r in c.records
                println(io, shifted_line(r, contig_name(set, c), c.ws, symbolic_twin, contigs[contig_name(set, c)]))
            end
        end
    end
end

function write_fixtures(destination)
    rng = MersenneTwister(SEED)
    reader = open(FASTA.Reader, GENOME; index = GENOME * ".fai")
    mkpath(joinpath(destination, "guides"))
    contigs = Pair{String, String}[]
    rows = String["set\tcontig\tsource\trecords\tgenotype_oracle\tguides\ttags"]
    readme = ["# Real-data VCF fixtures", "",
        "Generated by `make_real_vcfs.jl` from the large GRCh38 call sets used in `vcf.md`; do not edit by hand.",
        "Each contig of `reference.fa` is one GRCh38 window (`<set>_<chrom>_<start>`). Records keep REF, ALT, ID,",
        "GT and PS; POS and INFO/END are shifted to the window. `windows.tsv` lists the source coordinates, the",
        "tags each window covers, and its guides (Cas9, distance $DISTANCE).",
        "Tests: `test/variants/real.jl`. Expected outputs in `expected/` (gzipped CSV) are rewritten with",
        "`CHOPOFF_VCF_REAL_UPDATE=1` after the oracle checks pass.", "",
        "| Set | Source | Samples | Windows | Records | Bases |", "|---|---|---|---:|---:|---:|"]
    for src in SOURCES
        @info "Reading $(src.set)" src.vcf
        header, names, bychrom = read_source(src)
        cands = candidates(bychrom, reader, src)
        chosen = select_windows(cands)
        @info "$(src.set): $(length(cands)) candidates, $(length(chosen)) windows" tags = length(union(Set{String}(), (c.tags for c in chosen)...))
        flush(stderr)
        set_contigs = Dict{String, String}()
        guides = String[]
        for c in chosen
            name = contig_name(src.set, c)
            w = chromosome(reader, c.chrom)[c.ws:c.we]
            set_contigs[name] = w
            push!(contigs, name => w)
            lit = literal_records(c.records, w, c.ws)
            edits = window_edits(w, lit, 1, length(w))
            samples = window_samples(c.records, edits, names, w)
            haps = vtest_haplotypes(w, edits; chrom = "w", max_alleles = src.cap == 0 ? typemax(Int) : src.cap)
            g = guides_for(rng, w, edits, DISTANCE; n = 2, haps)
            if samples !== nothing
                # One more guide from a copy some sample carries.
                carried = []
                for s in samples
                    vtest_sample_truth(edits, s, taken -> (isempty(taken) || push!(carried, (vtest_haplotype(w, taken; chrom = "w")..., true)); false))
                end
                isempty(carried) || append!(g, guides_for(rng, w, edits, DISTANCE; n = 1, haps = carried))
            end
            g = unique(g)
            append!(guides, g)
            push!(rows, join([src.set, name, "$(c.chrom):$(c.ws)-$(c.we)", length(c.records),
                samples === nothing ? "no" : "yes", join(g, ','), join(sort(collect(c.tags)), ',')], '\t'))
        end
        write_set_vcf(joinpath(destination, "$(src.set).vcf"), header, names, chosen, src.set, set_contigs, false, src.vcf)
        if src.set == "hg002"
            write_set_vcf(joinpath(destination, "hg002_symbolic.vcf"), header, names, chosen, src.set, set_contigs, true, src.vcf)
        end
        write(joinpath(destination, "guides", "$(src.set).txt"), join(unique(guides), '\n') * "\n")
        push!(readme, "| $(src.set) | `$(basename(src.vcf))`$(src.region === nothing ? "" : " ($(src.region))") | " *
            "$(isempty(names) ? "none" : join(names, ", ")) | $(length(chosen)) | $(sum(c -> length(c.records), chosen; init = 0)) | " *
            "$(sum(c -> c.we - c.ws + 1, chosen; init = 0)) |")
    end
    close(reader)
    vtest_fasta(joinpath(destination, "reference.fa"), contigs)
    write(joinpath(destination, "windows.tsv"), join(rows, '\n') * "\n")
    append!(readme, ["", "`hg002_symbolic.vcf` is `hg002.vcf` with biallelic deletions of at least 50 bp written as `<DEL>` and",
        "insertions of at least 50 bp that repeat the following reference written as `<DUP:TANDEM>`.",
        "`gnomad` is searched with the default ALT-allele cap (3); the other sets with an exact cap."])
    write(joinpath(destination, "README.md"), join(readme, '\n') * "\n")
end

function check_fixtures()
    mktempdir() do temporary
        write_fixtures(temporary)
        for (root, _, files) in walkdir(temporary), file in files
            relative = relpath(joinpath(root, file), temporary)
            target = joinpath(@__DIR__, relative)
            isfile(target) && read(target) == read(joinpath(root, file)) || error("Stale real-data fixture: $relative")
        end
    end
    println("Real-data VCF fixtures match their sources.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? write_fixtures(@__DIR__) : ARGS == ["--check"] ? check_fixtures() : error("Expected no arguments or --check")
end
