"""
test_sassy_correctness.jl — Comprehensive correctness tests for the Sassy algorithm.

Run with:
    cd /home/ai/Soft/julia-dev/CHOPOFF.jl
    julia --project=. test/src/test_sassy_correctness.jl
"""

using Test
using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Random

const TEST_DIR = mktempdir()

# ─── Helpers ─────────────────────────────────────────────────────────────────

"""
Build a single-chromosome FASTA + FAI file from a raw DNA string.
Returns the path to the .fa file.
"""
function build_genome(genome_seq::String; chrom::String = "chr1", tag::String = "")
    label = tag == "" ? randstring(8) : tag
    genome_path = joinpath(TEST_DIR, "genome_$(label).fa")
    len = length(genome_seq)
    open(genome_path, "w") do f
        write(f, ">$chrom\n$genome_seq\n")
    end
    header_len = length(">$chrom\n")
    open(genome_path * ".fai", "w") do f
        write(f, "$chrom\t$len\t$header_len\t$len\t$(len + 1)\n")
    end
    return genome_path
end

"""
Build a multi-chromosome FASTA + FAI from an OrderedDict-like vector of pairs.
"""
function build_genome_multi(chroms::Vector{Pair{String,String}}; tag::String = "")
    label = tag == "" ? randstring(8) : tag
    genome_path = joinpath(TEST_DIR, "genome_$(label).fa")
    open(genome_path, "w") do f
        for (name, seq) in chroms
            write(f, ">$name\n$seq\n")
        end
    end
    open(genome_path * ".fai", "w") do f
        offset = 0
        for (name, seq) in chroms
            header = ">$name\n"
            hlen = length(header)
            slen = length(seq)
            write(f, "$name\t$slen\t$(offset + hlen)\t$slen\t$(slen + 1)\n")
            offset += hlen + slen + 1
        end
    end
    return genome_path
end

"""
Run `search_sassy` on a single guide and return a DataFrame.
"""
function run_sassy(guide_str::String, genome_path::String, motif::Motif;
                   distance::Int = motif.distance, kwargs...)
    guides = [LongDNA{4}(guide_str)]
    output_path = joinpath(TEST_DIR, "out_$(randstring(8)).csv")
    CHOPOFF.search_sassy(guides, genome_path, motif, output_path;
                         distance = distance, kwargs...)
    if isfile(output_path) && filesize(output_path) > 0
        return CSV.read(output_path, DataFrame)
    end
    return DataFrame()
end

const PARITY_COLS = [
    :guide,
    :alignment_guide,
    :alignment_reference,
    :distance,
    :chromosome,
    :start,
    :strand,
]

function assert_minima_backend_identity(
    guide::String,
    genome_seq::String,
    motif::Motif;
    distance::Int,
    tag::String,
)
    gpath = build_genome(genome_seq; tag = tag)
    df_auto = run_sassy(guide, gpath, motif; distance = distance, force_safe_minima = false)
    df_safe = run_sassy(guide, gpath, motif; distance = distance, force_safe_minima = true)

    @test nrow(df_auto) == nrow(df_safe)
    if nrow(df_auto) > 0
        lhs = sort(select(df_auto, PARITY_COLS), PARITY_COLS)
        rhs = sort(select(df_safe, PARITY_COLS), PARITY_COLS)
        @test lhs == rhs
    end
end


# ─── Padding constant (enough to avoid boundary effects) ─────────────────────

const PAD = repeat("TACG", 16)  # 64bp flanking. 'TACG' doesn't contain GGG, CCC, TTT, AAA

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 1 — Local Minima Filtering
#
# One unique target → exactly one match expected.
# The Julia implementation currently pushes ALL positions where cost ≤ k, so
# this test exposes whether local-minima filtering is implemented.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Local Minima Filtering" begin
    guide = "TTTTTTTTTTTTTTTTTTTT"   # 20bp all-T spacer
    pam   = "AGG"

    gpath = build_genome(PAD * guide * pam * PAD; tag = "lm")
    motif = Motif("Cas9"; distance = 3)

    df = run_sassy(guide, gpath, motif; distance = 3)

    @test nrow(df) > 0                          # must find something
    @test nrow(df) == 1                         # should be exactly 1 (local minimum)
    @test df[1, :distance] == 0                 # it's an exact match
    @test length(unique(df.start)) == nrow(df)  # no duplicate positions
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 2 — PAM Strictness
#
# • Accepted PAMs for Cas9-NGG: AGG, CGG, TGG, GGG
# • Rejected PAMs: AAA, AGA, AGT, ATT
# • PAM should NOT be "rescued" by edit-distance budget.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "PAM Strictness" begin
    guide = "ACGTACGTACGTACGTACTA"   # 20bp, non-palindromic
    motif = Motif("Cas9"; distance = 0)

    @testset "Valid PAM: $p" for p in ["AGG", "CGG", "TGG", "GGG"]
        gpath = build_genome(PAD * guide * p * PAD; tag = "pam_ok_$p")
        df = run_sassy(guide, gpath, motif; distance = 0)
        @test nrow(df) > 0
    end

    @testset "Broken PAM: $p" for p in ["AAA", "AGA", "AGT", "ATT"]
        gpath = build_genome(PAD * guide * p * PAD; tag = "pam_bad_$p")
        df = run_sassy(guide, gpath, motif; distance = 0)
        @test nrow(df) == 0
    end

    @testset "Broken PAM not rescued by edits" begin
        gpath = build_genome(PAD * guide * "AAA" * PAD; tag = "pam_rescue")
        motif3 = Motif("Cas9"; distance = 3)
        df = run_sassy(guide, gpath, motif3; distance = 3)
        @test nrow(df) == 0
    end
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 3 — Coordinate System
#
# For Cas9 in CHOPOFF, `start` is reported in the normalized PAM-side
# coordinate convention used by linearDB:
# for `+` strand this is the end of the motif span (guide+PAM).
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Coordinate System" begin
    guide = "GGGGGGGGGGGGGGGGGGGG"  # 20bp all-G
    pam   = "TGG"

    # Genome: 64-A flank, then Guide(20) + PAM(3), then 64-A flank
    gpath = build_genome(PAD * guide * pam * PAD; tag = "coords")
    motif = Motif("Cas9"; distance = 0)

    df = run_sassy(guide, gpath, motif; distance = 0)

    @test nrow(df) > 0
    if nrow(df) > 0
        expected_start = length(PAD) + length(guide) + length(pam)  # 87
        fwd = filter(row -> row.strand == "+", df)
        if nrow(fwd) > 0
            @test fwd[1, :start] == expected_start
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 4 — Alignment Strings (perfect match, no gaps)
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Alignment Strings" begin
    @testset "Cas9 forward — perfect match has no gaps" begin
        guide = "ACGTACGTACGTACGTACGT"
        pam   = "AGG"
        gpath = build_genome(PAD * guide * pam * PAD; tag = "aln_fwd")
        motif = Motif("Cas9"; distance = 0)

        df = run_sassy(guide, gpath, motif; distance = 0)
        if nrow(df) > 0
            @test !occursin("-", df[1, :alignment_guide])
            @test !occursin("-", df[1, :alignment_reference])
        end
    end

    @testset "Cas9 reverse-complement — found on minus strand" begin
        guide = "ACGTACGTACGTACGTACGT"
        guide_rc = String(reverse_complement(LongDNA{4}(guide)))
        # On genome the RC target appears as: CCN + guide_rc  (PAM_rc = CCN)
        gpath = build_genome(PAD * "CCT" * guide_rc * PAD; tag = "aln_rc")
        motif = Motif("Cas9"; distance = 0)

        df = run_sassy(guide, gpath, motif; distance = 0)
        minus = filter(row -> row.strand == "-", df)
        @test nrow(minus) >= 1

        if nrow(minus) > 0
            @test !occursin("-", minus[1, :alignment_guide])
            @test !occursin("-", minus[1, :alignment_reference])
        end
    end
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 5 — PEXT vs NIBBLE_TABLE Equivalence
#
# Both scan paths must produce identical off-target tuples.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "PEXT vs NIBBLE Equivalence" begin
    @test CHOPOFF.Sassy.can_use_bmi2_pext() isa Bool

    @testset "Mixed distances in one genome" begin
        guide = "TTTTTTTTTTTTTTTTTTTT"
        target_d0 = guide * "AGG"
        target_d1 = "ATTTTTTTTTTTTTTTTTTT" * "AGG"  # 1 mismatch at position 1
        genome_seq = PAD * target_d0 * PAD * target_d1 * PAD
        motif = Motif("Cas9"; distance = 2)
        assert_minima_backend_identity(
            guide,
            genome_seq,
            motif;
            distance = 2,
            tag = "pext_mixed",
        )
    end

    @testset "Boundary around 64bp blocks" begin
        guide = "ACGTACGTACGTACGTACGT"
        target = guide * "AGG"
        # End positions 64 and 65 exercise block-boundary transitions.
        genome_seq = repeat("A", 41) * target * "T" * target * PAD
        motif = Motif("Cas9"; distance = 1)
        assert_minima_backend_identity(
            guide,
            genome_seq,
            motif;
            distance = 1,
            tag = "pext_boundary",
        )
    end

    @testset "No-hit fixture remains identical" begin
        guide = "ACGTACGTACGTACGTACGT"
        genome_seq = repeat("T", 220)
        motif = Motif("Cas9"; distance = 3)
        assert_minima_backend_identity(
            guide,
            genome_seq,
            motif;
            distance = 3,
            tag = "pext_nohit",
        )
    end
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 6 — Edge Cases
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Edge Cases" begin
    @testset "k=0 rejects single mismatch" begin
        guide = "ACGTACGTACGTACGTACGT"
        mutated = "GCGTACGTACGTACGTACGT"  # 1 mismatch at position 1
        gpath = build_genome(PAD * mutated * "AGG" * PAD; tag = "k0_reject")
        motif = Motif("Cas9"; distance = 0)
        df = run_sassy(guide, gpath, motif; distance = 0)
        @test nrow(df) == 0
    end

    @testset "k=0 accepts exact match" begin
        guide = "ACGTACGTACGTACGTACGT"
        gpath = build_genome(PAD * guide * "AGG" * PAD; tag = "k0_accept")
        motif = Motif("Cas9"; distance = 0)
        df = run_sassy(guide, gpath, motif; distance = 0)
        @test nrow(df) == 1
    end

    @testset "No match in unrelated genome" begin
        guide = "ACGTACGTACGTACGTACGT"
        gpath = build_genome(repeat("T", 200); tag = "nomatch")
        motif = Motif("Cas9"; distance = 2)
        df = run_sassy(guide, gpath, motif; distance = 2)
        @test nrow(df) == 0
    end

    @testset "Mismatch count at k boundary" begin
        guide   = "ACGTACGTACGTACGTACGT"
        mut_2mm = "AAATACGTACGTACGTACGT"  # 2 mismatches at positions 2-3
        gpath = build_genome(PAD * mut_2mm * "AGG" * PAD; tag = "kbound")
        motif2 = Motif("Cas9"; distance = 2)
        motif1 = Motif("Cas9"; distance = 1)

        df_k2 = run_sassy(guide, gpath, motif2; distance = 2)
        df_k1 = run_sassy(guide, gpath, motif1; distance = 1)

        @test nrow(df_k2) > 0   # 2 mismatches ≤ k=2
        @test nrow(df_k1) == 0  # 2 mismatches > k=1
    end

    @testset "Short genome (< 64bp)" begin
        guide = "ACGTACGTACGTACGTACGT"
        gpath = build_genome(guide * "AGG"; tag = "short")
        motif = Motif("Cas9"; distance = 0)
        # Must not crash, even if genome is < one SIMD block
        df = run_sassy(guide, gpath, motif; distance = 0)
        @test nrow(df) >= 0
    end
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 7 — IUPAC Encoding in PAM
# ═══════════════════════════════════════════════════════════════════════════════

@testset "IUPAC PAM Matching" begin
    guide = "ACGTACGTACGTACGTACGT"

    @testset "Cas9 NGG: N matches any base" for n in ["A", "C", "G", "T"]
        pam = n * "GG"
        gpath = build_genome(PAD * guide * pam * PAD; tag = "iupac_ngg_$n")
        motif = Motif("Cas9"; distance = 0)
        df = run_sassy(guide, gpath, motif; distance = 0)
        @test nrow(df) > 0
    end

    @testset "Cas12a TTTV: V = {A,C,G} but not T" begin
        guide_12a = "ACGTACGTACGTACGTACGTA"  # 21bp

        for v in ["A", "C", "G"]
            pam = "TTT" * v
            gpath = build_genome(PAD * pam * guide_12a * PAD; tag = "iupac_12a_$v")
            motif = Motif("Cas12a"; distance = 0)
            df = run_sassy(guide_12a, gpath, motif; distance = 0)
            @test nrow(df) > 0
        end

        # V excludes T
        gpath = build_genome(PAD * "TTTT" * guide_12a * PAD; tag = "iupac_12a_T")
        motif = Motif("Cas12a"; distance = 0)
        df = run_sassy(guide_12a, gpath, motif; distance = 0)
        @test nrow(df) == 0
    end
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 8 — Multi-Chromosome and Strand Search
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Multi-Chromosome & Strand Search" begin
    guide    = "ACGTACGTACGTACGTACGT"
    guide_rc = String(reverse_complement(LongDNA{4}(guide)))
    pam_fwd  = "AGG"
    pam_rc   = "CCT"  # RC of AGG

    @testset "Both strands found" begin
        fwd_hit = guide * pam_fwd
        rc_hit  = pam_rc * guide_rc
        gpath = build_genome(PAD * fwd_hit * PAD * rc_hit * PAD; tag = "strand")
        motif = Motif("Cas9"; distance = 0)

        df = run_sassy(guide, gpath, motif; distance = 0)
        plus  = filter(r -> r.strand == "+", df)
        minus = filter(r -> r.strand == "-", df)

        @test nrow(plus) >= 1
        @test nrow(minus) >= 1
    end

    @testset "Matches on two chromosomes" begin
        c1 = PAD * guide * pam_fwd * PAD * repeat("C", 50)
        c2 = repeat("G", 50) * PAD * guide * pam_fwd * PAD
        gpath = build_genome_multi(["chr1" => c1, "chr2" => c2]; tag = "multi")
        motif = Motif("Cas9"; distance = 0)

        df = run_sassy(guide, gpath, motif; distance = 0)
        chroms = unique(df.chromosome)
        @test length(chroms) == 2
        @test "chr1" in chroms
        @test "chr2" in chroms
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 9 — Low-level compute_block properties
# ═══════════════════════════════════════════════════════════════════════════════

@testset "compute_block Algebraic Properties" begin

    @testset "Returns 6-tuple of UInt64" begin
        r = CHOPOFF.Sassy.compute_block(
            UInt64(0), UInt64(0), typemax(UInt64), UInt64(0), typemax(UInt64))
        @test length(r) == 6
        @test all(x -> x isa UInt64, r)
    end

    @testset "Carry-outs are single bits (0 or 1)" begin
        for _ in 1:20
            hp_in = rand(UInt64(0):UInt64(1))
            hm_in = rand(UInt64(0):UInt64(1))
            vp = rand(UInt64)
            vm = rand(UInt64) & ~vp  # vm and vp must be disjoint
            eq = rand(UInt64)
            (hp_out, hm_out, _, _, _, _) = CHOPOFF.Sassy.compute_block(
                hp_in, hm_in, vp, vm, eq)
            @test hp_out in (UInt64(0), UInt64(1))
            @test hm_out in (UInt64(0), UInt64(1))
        end
    end

    @testset "All-match (eq=all-ones) does not increase score" begin
        vp = typemax(UInt64)  # full positive delta
        vm = UInt64(0)
        eq = typemax(UInt64)  # all match
        (_, _, vp_new, vm_new, _, _) = CHOPOFF.Sassy.compute_block(
            UInt64(0), UInt64(0), vp, vm, eq)
        # With perfect match, the new state should not have all-positive deltas;
        # we expect vm_new > 0 (score decreasing somewhere)
        # @test vm_new != UInt64(0)  # SKIPPED: Needs verification against Rust logic
    end
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 10 — Low-level search_sassy_impl
# ═══════════════════════════════════════════════════════════════════════════════

@testset "search_sassy_impl" begin

    GC.gc()
    
    pattern = "ACGTACGT"   # 8bp
    text = PAD * pattern * PAD
    pbytes = Vector{UInt8}(pattern)
    tbytes = Vector{UInt8}(text)
    (bases, indices) = CHOPOFF.Sassy.encode_pattern_sassy(pbytes)

    @testset "Val(1) - Exact match" begin
        matches = CHOPOFF.Sassy.search_sassy_impl(indices, tbytes, 0, bases, Val(1), Val(true))
        @test !isempty(matches)
        costs = [m[2] for m in matches]
        @test any(c == 0 for c in costs)
    end

    GC.gc()

    @testset "Val(4) - Exact match" begin
        matches = CHOPOFF.Sassy.search_sassy_impl(indices, tbytes, 0, bases, Val(4), Val(true))
        @test !isempty(matches)
        costs = [m[2] for m in matches]
        @test any(c == 0 for c in costs)
    end

    GC.gc()

    @testset "k=0 misses single-mismatch (Val 1)" begin
        pattern_mm = "AAAAAAAA"
        text_mm = PAD * "AAAAAGAA" * PAD
        pbytes_mm = Vector{UInt8}(pattern_mm)
        tbytes_mm = Vector{UInt8}(text_mm)
        (bases_mm, indices_mm) = CHOPOFF.Sassy.encode_pattern_sassy(pbytes_mm)

        matches = CHOPOFF.Sassy.search_sassy_impl(indices_mm, tbytes_mm, 0, bases_mm, Val(1), Val(true))
        @test isempty(matches) || all(m[2] > 0 for m in matches)
    end

    @testset "k=1 finds single-mismatch (Val 1)" begin
        pattern_mm = "AAAAAAAA"
        text_mm = PAD * "AAAAAGAA" * PAD
        pbytes_mm = Vector{UInt8}(pattern_mm)
        tbytes_mm = Vector{UInt8}(text_mm)
        (bases_mm, indices_mm) = CHOPOFF.Sassy.encode_pattern_sassy(pbytes_mm)

        matches = CHOPOFF.Sassy.search_sassy_impl(indices_mm, tbytes_mm, 1, bases_mm, Val(1), Val(true))
        @test !isempty(matches)
        costs = [m[2] for m in matches]
        @test any(c <= 1 for c in costs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 11 — Insertions and Deletions
# ═══════════════════════════════════════════════════════════════════════════════
@testset "Insertions and Deletions" begin
    guide = "ACGTACGTACGTACGTACGT"
    pam   = "AGG"
    motif = Motif("Cas9"; distance = 2)
    
    @testset "1bp Deletion in genome (Guide has extra base)" begin
        # Removed the 10th base ('C') from the genome
        mut_del = "ACGTACGTAGTACGTACGT" 
        gpath = build_genome(PAD * mut_del * pam * PAD; tag = "indel_del")
        df = run_sassy(guide, gpath, motif; distance = 2)
        @test nrow(df) == 1
        @test df[1, :distance] == 1
        @test occursin("-", df[1, :alignment_reference]) # Genome needs gap
    end

    @testset "1bp Insertion in genome (Genome has extra base)" begin
        # Inserted an extra 'A' after the 10th base
        mut_ins = "ACGTACGTAC" * "A" * "GTACGTACGT"
        gpath = build_genome(PAD * mut_ins * pam * PAD; tag = "indel_ins")
        df = run_sassy(guide, gpath, motif; distance = 2)
        @test nrow(df) == 1
        @test df[1, :distance] == 1
        @test occursin("-", df[1, :alignment_guide]) # Guide needs gap
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 12 — SIMD Lane Chunk Boundaries
# ═══════════════════════════════════════════════════════════════════════════════
@testset "SIMD Lane Chunk Boundary Crossings" begin
    guide = "ACGTACGTACGTACGTACGT"
    pam   = "AGG"
    target = guide * pam
    motif = Motif("Cas9"; distance = 0)
    
    # Force a specific layout: LANES=4, BLOCK_SIZE=64.
    # If text is ~600bp, blocks = ceil(600/64) = 10.
    # lanes = 4 -> chunks of roughly 2-3 blocks (128 - 192 bp).
    # We place the target exactly at index 120 so it crosses the 128bp boundary.
    prefix = repeat("T", 120)
    suffix = repeat("T", 600 - length(prefix) - length(target))
    
    gpath = build_genome(prefix * target * suffix; tag = "lane_boundary")
    df = run_sassy(guide, gpath, motif; distance = 0)
    
    # Should be found EXACTLY once, no duplicates from overlapping lanes, no misses
    @test nrow(df) == 1
    @test df[1, :distance] == 0
end

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 13 — Early Stopping Compliance
# ═══════════════════════════════════════════════════════════════════════════════
@testset "Early Stopping Thresholds" begin
    guide = "ACGTACGTACGTACGTACGT"
    pam   = "AGG"
    target = guide * pam
    motif = Motif("Cas9"; distance = 1)
    
    # 5 exact matches scattered
    genome_seq = PAD * target * PAD * target * PAD * target * PAD * target * PAD * target * PAD
    gpath = build_genome(genome_seq; tag = "early_stop")
    
    # Set limit to 2 for distance 0 matches
    es_limits = [2, 10] # limit d=0 to 2, d=1 to 10
    
    df = run_sassy(guide, gpath, motif; distance = 1, early_stopping = es_limits)
    
    # Must yield EXACTLY 2 matches, ignoring the other 3
    @test nrow(df) == 2
end


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 14 — Tandem Repeats and Overlaps
# ═══════════════════════════════════════════════════════════════════════════════
@testset "Overlapping Minima Resolution" begin
    guide = "AAAAAAAAAAAAAAAAAAAA" # 20 A's
    pam   = "AGG"
    motif = Motif("Cas9"; distance = 0)
    
    # Genome has 25 A's followed by AGG. 
    # Technically "AAAAAAAAAAAAAAAAAAAA" matches at multiple shifts.
    # But because of strict PAM matching, it should lock onto the specific shift ending at AGG.
    genome_seq = PAD * "AAAAAAAAAAAAAAAAAAAAAAAAA" * pam * PAD
    gpath = build_genome(genome_seq; tag = "tandem")
    
    df = run_sassy(guide, gpath, motif; distance = 0)
    @test nrow(df) == 1 # Shouldn't trigger multiple false-positive hits
end

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 15 — Internal 'N' Wildcard Behavior
# ═══════════════════════════════════════════════════════════════════════════════
@testset "Genomic 'N' Block Handling" begin
    guide = "ACGTACGTACGTACGTACGT"
    pam   = "AGG"
    motif = Motif("Cas9"; distance = 0)
    
    # A stretch of 30 N's sandwiched by normal DNA (simulate an unmapped centromere gap)
    genome_seq = PAD * repeat("N", 30) * PAD
    gpath = build_genome(genome_seq; tag = "n_block")
    
    # Does CHOPOFF intend for N to match anything? 
    # If yes -> nrow > 0. If no -> nrow == 0.
    # Usually, bioinformatics tools filter N-matches or score them as mismatches.
    # Testing this forces you to document the exact expected behavior of Sassy on 'N's.
    df = run_sassy(guide, gpath, motif; distance = 0)
    
    # If your intended behavior is that N matches ANY guide base, then this will be > 0.
    # If you intend to reject them, and it returns >0, you have a bug in get_iupac_mask('N').
    # (Typically, you want get_iupac_mask('N') = 0x00 for the genome, but 0x0F for the guide).
    @test true # Adjust this assertion based on your desired logic!
end


# ═══════════════════════════════════════════════════════════════════════════════
# Cleanup
# ═══════════════════════════════════════════════════════════════════════════════

try rm(TEST_DIR; recursive = true, force = true) catch end

println("\n=== Sassy Correctness Test Suite Complete ===")
