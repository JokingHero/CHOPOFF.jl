using Test
using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Random

const TRACEBACK_TEST_DIR = mktempdir()
const TRACEBACK_CORE_COLS = [:guide, :distance, :chromosome, :start, :strand]
const TRACEBACK_DNA_ALPHABET = ('A', 'C', 'G', 'T')
const TRACEBACK_PAD = repeat("ACGT", 24)
const TRACEBACK_CAS12A_GUIDES = LongDNA{4}.([
    "TCGATTGTTTGGCTCTCTAAA",
    "GCAGGGGGACGCAAGTACGAA",
    "GGGCCGAAACGCGACACCGCC",
])

function build_traceback_genome(genome_seq::String; chrom::String = "chr1", tag::String = "")
    label = tag == "" ? randstring(8) : tag
    genome_path = joinpath(TRACEBACK_TEST_DIR, "genome_$(label).fa")
    len = Base.length(genome_seq)
    open(genome_path, "w") do io
        write(io, ">$chrom\n$genome_seq\n")
    end
    header_len = Base.length(">$chrom\n")
    open(genome_path * ".fai", "w") do io
        write(io, "$chrom\t$len\t$header_len\t$len\t$(len + 1)\n")
    end
    return genome_path
end

function run_traceback_backend(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif;
    distance::Int,
    strict_pam::Bool,
    traceback_backend::Symbol,
)
    output_path = joinpath(
        TRACEBACK_TEST_DIR,
        "traceback_$(randstring(10))_$(traceback_backend)_$(distance)_$(strict_pam).csv",
    )
    CHOPOFF.search_sassy(
        guides,
        genome_path,
        motif,
        output_path;
        distance = distance,
        strict_pam = strict_pam,
        traceback_backend = traceback_backend,
    )
    return DataFrame(CSV.File(output_path))
end

function normalize_traceback_core(df::DataFrame)
    cols = Set(Symbol.(names(df)))
    for c in TRACEBACK_CORE_COLS
        c in cols || error("Missing required parity column '$c'.")
    end
    core = select(df, TRACEBACK_CORE_COLS)
    core.guide = String.(core.guide)
    core.distance = Int.(core.distance)
    core.chromosome = String.(core.chromosome)
    core.start = Int.(core.start)
    core.strand = String.(core.strand)
    sort!(core, TRACEBACK_CORE_COLS)
    return core
end

function assert_traceback_core_parity(
    align_df::DataFrame,
    custom_df::DataFrame,
    label::AbstractString,
)
    lhs = normalize_traceback_core(align_df)
    rhs = normalize_traceback_core(custom_df)
    if lhs != rhs
        lhs_only = antijoin(lhs, rhs, on = TRACEBACK_CORE_COLS)
        rhs_only = antijoin(rhs, lhs, on = TRACEBACK_CORE_COLS)
        println("Traceback core parity mismatch ($label)")
        println("align-only rows: ", nrow(lhs_only))
        if nrow(lhs_only) > 0
            println(first(lhs_only, min(nrow(lhs_only), 10)))
        end
        println("custom-only rows: ", nrow(rhs_only))
        if nrow(rhs_only) > 0
            println(first(rhs_only, min(nrow(rhs_only), 10)))
        end
    end
    @test lhs == rhs
end

@inline function random_dna(rng::AbstractRNG, len::Int)
    return String([rand(rng, TRACEBACK_DNA_ALPHABET) for _ in 1:len])
end

function mutate_substitutions(seq::String, edits::Int, rng::AbstractRNG)
    edits == 0 && return seq
    chars = collect(seq)
    positions = randperm(rng, Base.length(chars))[1:min(edits, Base.length(chars))]
    for pos in positions
        old = chars[pos]
        new_base = old
        while new_base == old
            new_base = rand(rng, TRACEBACK_DNA_ALPHABET)
        end
        chars[pos] = new_base
    end
    return String(chars)
end

@testset "Traceback Backend API Validation" begin
    guide = LongDNA{4}("ACGTACGTACGTACGTACGT")
    genome_path = build_traceback_genome(TRACEBACK_PAD * String(guide) * "AGG" * TRACEBACK_PAD; tag = "bad_backend")
    motif = Motif("Cas9"; distance = 1)
    output_path = joinpath(TRACEBACK_TEST_DIR, "bad_backend_out.csv")
    @test_throws ArgumentError CHOPOFF.search_sassy(
        [guide],
        genome_path,
        motif,
        output_path;
        distance = 1,
        traceback_backend = :not_a_backend,
    )
end

@testset "Traceback Backend Parity (Sample Fixtures, strict_pam=true)" begin
    data_root = joinpath(dirname(pathof(CHOPOFF)), "..", "test", "sample_data")
    cas9_genome = joinpath(data_root, "genome", "semirandom.fa")
    cas9_guides_all = LongDNA{4}.(readlines(joinpath(data_root, "guides.txt")))
    cas9_guides = cas9_guides_all[1:min(end, 4)]
    cases = [
        (
            label = "Cas9",
            motif_name = "Cas9",
            genome = cas9_genome,
            guides = cas9_guides,
        ),
        (
            label = "Cas12a",
            motif_name = "Cas12a",
            genome = joinpath(data_root, "genome", "semirandom.2bit"),
            guides = TRACEBACK_CAS12A_GUIDES,
        ),
    ]

    for case in cases
        for d in 1:3
            motif = Motif(case.motif_name; distance = d)
            align_df = run_traceback_backend(
                case.guides,
                case.genome,
                motif;
                distance = d,
                strict_pam = true,
                traceback_backend = :align,
            )
            custom_df = run_traceback_backend(
                case.guides,
                case.genome,
                motif;
                distance = d,
                strict_pam = true,
                traceback_backend = :custom,
            )
            assert_traceback_core_parity(
                align_df,
                custom_df,
                "$(case.label) d=$d strict=true",
            )
        end
    end
end

@testset "Traceback Backend Parity (Synthetic strict_pam=false)" begin
    rng = MersenneTwister(20260222)
    cases = [
        (label = "Cas9", motif_name = "Cas9", guide = "ACGTACGTACGTACGTACGT", pam = "AGG"),
        (label = "Cas12a", motif_name = "Cas12a", guide = "TGCATGCATGCATGCATGCAT", pam = "TTTA"),
    ]

    for case in cases
        for d in 1:3
            motif = Motif(case.motif_name; distance = d)
            mut_guide = mutate_substitutions(case.guide, min(d, 2), rng)
            site_exact = motif.extends5 ? (case.guide * case.pam) : (case.pam * case.guide)
            site_mut = motif.extends5 ? (mut_guide * case.pam) : (case.pam * mut_guide)
            genome_seq = TRACEBACK_PAD * site_exact * TRACEBACK_PAD * site_mut * TRACEBACK_PAD
            genome_path = build_traceback_genome(genome_seq; tag = "strict_false_$(case.label)_$d")
            guides = [LongDNA{4}(case.guide)]

            align_df = run_traceback_backend(
                guides,
                genome_path,
                motif;
                distance = d,
                strict_pam = false,
                traceback_backend = :align,
            )
            custom_df = run_traceback_backend(
                guides,
                genome_path,
                motif;
                distance = d,
                strict_pam = false,
                traceback_backend = :custom,
            )
            assert_traceback_core_parity(
                align_df,
                custom_df,
                "$(case.label) d=$d strict=false",
            )
        end
    end
end

@testset "Traceback Backend Randomized Core Parity" begin
    rng = MersenneTwister(12345)
    motif_cases = [("Cas9", 20, "AGG"), ("Cas12a", 21, "TTTA")]

    for (motif_name, guide_len, pam) in motif_cases
        for strict_pam in (true, false)
            for d in 1:3
                for trial in 1:2
                    guide = random_dna(rng, guide_len)
                    edits = rand(rng, 0:d)
                    mutated = mutate_substitutions(guide, edits, rng)
                    site = motif_name == "Cas9" ? (mutated * pam) : (pam * mutated)
                    genome_seq = random_dna(rng, 90) * site * random_dna(rng, 90)
                    genome_path = build_traceback_genome(
                        genome_seq;
                        tag = "rand_$(motif_name)_$(strict_pam)_$(d)_$(trial)",
                    )

                    motif = Motif(motif_name; distance = d)
                    guides = [LongDNA{4}(guide)]
                    align_df = run_traceback_backend(
                        guides,
                        genome_path,
                        motif;
                        distance = d,
                        strict_pam = strict_pam,
                        traceback_backend = :align,
                    )
                    custom_df = run_traceback_backend(
                        guides,
                        genome_path,
                        motif;
                        distance = d,
                        strict_pam = strict_pam,
                        traceback_backend = :custom,
                    )
                    assert_traceback_core_parity(
                        align_df,
                        custom_df,
                        "random motif=$motif_name strict=$(strict_pam) d=$d trial=$trial",
                    )
                end
            end
        end
    end
end

try
    rm(TRACEBACK_TEST_DIR; recursive = true, force = true)
catch
end
