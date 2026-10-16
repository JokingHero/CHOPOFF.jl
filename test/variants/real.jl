# Real-data fixtures cut from GRCh38 call sets (test/sample_data/vcf/real/README.md).
# Per contig, union rows equal the haplotype oracle and, where windows.tsv allows,
# genotype rows pass the per-copy oracle. Symbolic SVs are checked through a literal
# twin: the oracle runs on the twin and the symbolic rows must match it once allele
# keys are renamed. Per set, the output equals expected/ (gzipped: allele keys hold
# whole ALT sequences); CHOPOFF_VCF_REAL_UPDATE=1 rewrites expected/ instead.
isdefined(@__MODULE__, :vtest_rc) || include(joinpath(@__DIR__, "..", "helpers", "vcf_sweeps.jl"))
isdefined(@__MODULE__, :vtest_compare_genotypes) || include(joinpath(@__DIR__, "..", "helpers", "vcf_genotypes.jl"))
isdefined(@__MODULE__, :window_edits) || include(joinpath(@__DIR__, "..", "helpers", "vcf_real.jl"))

const VREAL_DIR = joinpath(@__DIR__, "..", "sample_data", "vcf", "real")
const VREAL_DISTANCE = 3
# cap 0: exact (every record of a contig may combine); gnomad uses the package default.
const VREAL_SETS = [
    (set = "hg002", windows = "hg002", genotype = true, cap = 0),
    (set = "hg002_symbolic", windows = "hg002", genotype = true, cap = 0),
    (set = "hg001", windows = "hg001", genotype = true, cap = 0),
    (set = "kg_panel", windows = "kg_panel", genotype = true, cap = 0),
    (set = "kg_ps", windows = "kg_ps", genotype = true, cap = 0),
    (set = "kg_sv", windows = "kg_sv", genotype = true, cap = 0),
    (set = "clinvar", windows = "clinvar", genotype = false, cap = 0),
    (set = "gnomad", windows = "gnomad", genotype = false, cap = CHOPOFF.PREFIX_HASH_SCAN_VCF_DEFAULT_MAX_ALLELES),
]

function vreal_golden(path, name)
    target = joinpath(VREAL_DIR, "expected", name * ".gz")
    if get(ENV, "CHOPOFF_VCF_REAL_UPDATE", "0") == "1"
        mkpath(dirname(target))
        write(target, transcode(GzipCompressor, read(path)))
    else
        @test isfile(target) && read(path) == transcode(GzipDecompressor, read(target))
    end
end

vreal_lines(path) = sort(readlines(path))

# Writes the contig's records, with symbolic SVs written literally when `literal`.
# Returns the literal records and the symbolic-to-literal allele key renames.
function vreal_contig(vcf, contig, seq, path; literal)
    header = filter(l -> startswith(l, "##") && !startswith(l, "##fileformat") &&
        (!startswith(l, "##contig") || occursin("ID=$contig,", l)), readlines(vcf))
    names, records = vreal_records(vcf; contig)
    renames = Pair{String, String}[]
    rows = Vector{String}[]
    lit = map(records) do r
        f = copy(r.fields)
        sv = findfirst(a -> startswith(a, '<'), r.alts)
        sv === nothing && (push!(rows, f); return r)
        ref, alt = vreal_literal(seq, r.pos, r.alts[sv])
        push!(renames, vtest_escape("$contig:$(r.pos):$(f[4])>$(split(f[5], ',')[sv])/$(split(r.alts[sv], '/')[end])") =>
            vtest_escape("$contig:$(r.pos):$ref>$alt"))
        literal && ((f[4], f[5]) = (ref, alt))
        push!(rows, f)
        (; r..., ref, alts = [i == sv ? alt : a for (i, a) in enumerate(r.alts)])
    end
    vtest_vcf(path, rows; header, samples = names)
    return names, lit, renames
end

# The symbolic search must give the literal rows once allele keys are renamed.
function vreal_symbolic_parity(dir, contig, seq, guides, motif, literal_output; kwargs...)
    vcf = joinpath(dir, "symbolic.vcf")
    _, _, renames = vreal_contig(joinpath(VREAL_DIR, contig.vcf), contig.name, seq, vcf; literal = false)
    output = joinpath(dir, "symbolic.csv")
    CHOPOFF.search_prefixHashScan(guides, joinpath(dir, "reference.fa"), motif, output; distance = VREAL_DISTANCE, vcf,
        early_stopping = fill(typemax(Int), VREAL_DISTANCE + 1), kwargs...)
    renamed = joinpath(dir, "symbolic_renamed.csv")
    write(renamed, replace(read(output, String), renames...))
    @test vreal_lines(renamed) == vreal_lines(literal_output)
end

function vreal_check_contig(spec, w, seq)
    motif = Motif("Cas9"; distance = VREAL_DISTANCE)
    guides = LongDNA{4}.(split(w.guides, ','))
    contig = (name = w.contig, vcf = spec.set * ".vcf")
    mktempdir() do dir
        vcf = joinpath(dir, "literal.vcf")
        names, lit, renames = vreal_contig(joinpath(VREAL_DIR, contig.vcf), w.contig, seq, vcf; literal = true)
        edits = window_edits(seq, lit, 1, length(seq))
        cap = spec.cap == 0 ? vtest_exact_cap(edits) : spec.cap
        vtest_compare(dir, seq, edits, guides, motif, VREAL_DISTANCE; chrom = w.contig, vcf,
            max_alleles = cap, case_id = "real_$(w.contig)")
        isempty(renames) || vreal_symbolic_parity(dir, contig, seq, guides, motif, joinpath(dir, "actual.csv");
            vcf_max_alleles = cap)
        (spec.genotype && w.genotype_oracle == "yes") || return
        samples = map(names) do name
            calls = Dict{String, Any}("g$i" => r.calls[name] for (i, r) in enumerate(lit)
                if any(!=(0), r.calls[name][1]) && any(!isempty, r.alts))
            (; name, ploidy = length(first(lit).calls[name][1]), calls)
        end
        gdir = mkpath(joinpath(dir, "genotype"))
        vtest_compare_genotypes(gdir, seq, edits, samples, guides, motif, VREAL_DISTANCE; chrom = w.contig, vcf,
            case_id = "realgt_$(w.contig)")
        isempty(renames) || vreal_symbolic_parity(gdir, contig, seq, guides, motif, joinpath(gdir, "genotype.csv");
            vcf_mode = :genotype)
    end
end

@testset "Real-data fixtures" begin
    windows = CSV.File(joinpath(VREAL_DIR, "windows.tsv"); delim = '\t', types = String)
    genome = joinpath(VREAL_DIR, "reference.fa")
    sequences = open(FASTA.Reader, genome) do reader
        Dict(identifier(r) => uppercase(sequence(String, r)) for r in reader)
    end
    motif = Motif("Cas9"; distance = VREAL_DISTANCE)
    for spec in VREAL_SETS
        @testset "$(spec.set)" begin
            ws = [w for w in windows if w.set == spec.windows]
            @test !isempty(ws)
            for w in ws
                @testset "$(w.contig)" begin
                    vreal_check_contig(spec, w, sequences[w.contig])
                end
            end
            guides = LongDNA{4}.(readlines(joinpath(VREAL_DIR, "guides", spec.windows * ".txt")))
            vcf = joinpath(VREAL_DIR, spec.set * ".vcf")
            region = join((w.contig for w in ws), ',')
            cap = spec.cap == 0 ? maximum(w -> parse(Int, w.records), ws) : spec.cap
            stops = fill(typemax(Int), VREAL_DISTANCE + 1)
            mktempdir() do dir
                output = joinpath(dir, "union.csv")
                stats = PrefixHashScanStats()
                CHOPOFF.search_prefixHashScan(guides, genome, motif, output; distance = VREAL_DISTANCE, vcf, region,
                    vcf_max_alleles = cap, early_stopping = stops, stats)
                rows = DataFrame(CSV.File(output))
                @test any(!ismissing, rows.variants)
                vreal_golden(output, "$(spec.set)_union.csv")
                spec.set == "kg_sv" && @test stats.vcf_sv_searched > 0 && stats.vcf_sv_skipped > 0
                spec.set == "gnomad" && @test stats.vcf_truncated_anchors > 0
                spec.genotype || return
                output = joinpath(dir, "genotype.csv")
                CHOPOFF.search_prefixHashScan(guides, genome, motif, output; distance = VREAL_DISTANCE, vcf, region,
                    vcf_mode = :genotype, early_stopping = stops)
                rows = DataFrame(CSV.File(output; types = Dict(:carriers => String, :lost => String)))
                @test any(!ismissing, rows.carriers)
                vreal_golden(output, "$(spec.set)_genotype.csv")
                summary = joinpath(dir, "per_sample.csv")
                CSV.write(summary, summarize_offtargets(rows; per_sample = true, vcf, distance = VREAL_DISTANCE))
                vreal_golden(summary, "$(spec.set)_per_sample.csv")
            end
        end
    end
end
