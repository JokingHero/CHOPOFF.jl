using Test
using CHOPOFF
using BioSequences
using CSV
using DataFrames

const PHS_CORE_COLS = [:guide, :distance, :chromosome, :start, :strand]

function phs_core(path::String)
    df = DataFrame(CSV.File(path))
    core = select(df, PHS_CORE_COLS)
    core.guide = String.(core.guide)
    core.distance = Int.(core.distance)
    core.chromosome = String.(core.chromosome)
    core.start = Int.(core.start)
    core.strand = String.(core.strand)
    sort!(core, PHS_CORE_COLS)
    return core
end


function check_prefix_helper_matches_materialized(tdir::String, seq::String, motif::Motif, hash_len::Int, label::String)
    genome = joinpath(tdir, label * ".fa")
    write_phs_fasta(genome, "chr1", seq)
    chrom_seq = LongDNA{4}(seq)
    dbi = DBInfo(genome, "prefix_hash_scan_helper", motif)
    for is_antisense in (false, true)
        positions = CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
        isempty(positions) && continue
        candidate_range = first(positions)
        prefix = CHOPOFF.normalized_candidate_prefix(chrom_seq, candidate_range, dbi, is_antisense, hash_len)
        ot, _ = CHOPOFF.materialize_normalized_candidate(chrom_seq, candidate_range, dbi, is_antisense)
        @test prefix == ot[1:hash_len]
    end
end


function normalize_phs_query_map(query::Dict)
    out = Dict{eltype(keys(query)), Vector{Int}}()
    for (k, v) in query
        out[k] = sort(v)
    end
    return out
end

function normalize_phs_query_map(query::CHOPOFF.PrefixHashScanBitmaskQuery)
    out = Dict{eltype(keys(query.masks)), Vector{Int}}()
    for (k, mask0) in query.masks
        mask = mask0
        guides = Int[]
        while mask != 0
            push!(guides, trailing_zeros(mask) + 1)
            mask &= mask - 1
        end
        out[k] = guides
    end
    return out
end

function phs_hash_type(hash_len::Int)
    return CHOPOFF.smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
end

function write_phs_fasta(path::String, name::String, seq::String)
    open(path, "w") do io
        write(io, ">", name, "\n", seq, "\n")
    end
    open(path * ".fai", "w") do io
        write(io, name, "\t", string(length(seq)), "\t", string(length(name) + 2),
            "\t", string(length(seq)), "\t", string(length(seq) + 1), "\n")
    end
end

@testset "prefixHashScan" begin
    tdir = tempname()
    mkpath(tdir)


    @testset "exact Cas9 asset matches full fallback" begin
        motif = Motif("Cas9"; distance = 3)
        exact = CHOPOFF.load_precomputed_prefix_paths(motif, 3, 16; need_distances = true)
        scan_only = CHOPOFF.load_precomputed_prefix_paths(motif, 3, 16; need_distances = false)
        @test exact !== nothing
        @test scan_only !== nothing
        exact_paths, exact_distances, exact_asset = exact
        scan_paths, scan_distances, scan_asset = scan_only
        @test exact_asset == "Cas9"
        @test scan_asset == "Cas9"
        @test scan_distances === nothing
        @test scan_paths == exact_paths
        @test length(exact_distances) == size(exact_paths, 1)

        dir = joinpath(dirname(pathof(CHOPOFF)), "..", "data")
        full_paths = vcat(
            CHOPOFF.load(joinpath(dir, "Cas9_d4_p16_paths_part1.bin")),
            CHOPOFF.load(joinpath(dir, "Cas9_d4_p16_paths_part2.bin")),
        )
        full_distances = CHOPOFF.load(joinpath(dir, "Cas9_d4_p16_distances.bin"))
        keep = BitVector(full_distances .<= 3)
        @test exact_paths == full_paths[keep, 1:16]
        @test exact_distances == full_distances[keep]
    end

    @testset "precomputed Cas9 paths" begin
        stats = CHOPOFF.PrefixHashScanStats()
        paths, source = CHOPOFF.load_prefix_hash_scan_paths(Motif("Cas9"; distance = 3), 3, 16, stats)
        @test source == :precomputed
        @test stats.path_source == :precomputed
        @test size(paths, 2) == 16
        @test size(paths, 1) > 0
    end


    @testset "query hash variants match" begin
        motif = Motif("Cas9"; distance = 3)
        hash_len = 16
        hash_type = phs_hash_type(hash_len)
        paths, _ = CHOPOFF.load_prefix_hash_scan_paths(motif, 3, hash_len)
        guides = LongDNA{4}.([
            "ACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCA",
        ])
        guides_ = CHOPOFF.oriented_prefix_hash_scan_guides(guides, motif)

        for guide in guides_
            baseline = CHOPOFF.prefix_hash_scan_guide_hashes(paths, guide, hash_type; query_variant = :baseline)
            columnwise = CHOPOFF.prefix_hash_scan_guide_hashes(paths, guide, hash_type; query_variant = :columnwise)
            @test columnwise == baseline
        end

        baseline_stats = CHOPOFF.PrefixHashScanStats()
        columnwise_stats = CHOPOFF.PrefixHashScanStats()
        bitmask_stats = CHOPOFF.PrefixHashScanStats()
        auto_stats = CHOPOFF.PrefixHashScanStats()
        baseline_query = CHOPOFF.build_prefix_hash_scan_map_from_paths(paths, guides_, hash_type, baseline_stats; query_variant = :baseline)
        columnwise_query = CHOPOFF.build_prefix_hash_scan_map_from_paths(paths, guides_, hash_type, columnwise_stats; query_variant = :columnwise)
        bitmask_query = CHOPOFF.build_prefix_hash_scan_map_from_paths(paths, guides_, hash_type, bitmask_stats; query_variant = :bitmask64)
        auto_query = CHOPOFF.build_prefix_hash_scan_map_from_paths(paths, guides_, hash_type, auto_stats; query_variant = :auto)
        baseline_norm = normalize_phs_query_map(baseline_query)
        @test normalize_phs_query_map(columnwise_query) == baseline_norm
        @test normalize_phs_query_map(bitmask_query) == baseline_norm
        @test normalize_phs_query_map(auto_query) == baseline_norm
        @test baseline_stats.query_variant == :baseline
        @test columnwise_stats.query_variant == :columnwise
        @test bitmask_stats.query_variant == :bitmask64
        @test auto_stats.query_variant == :bitmask64
        @test bitmask_query isa CHOPOFF.PrefixHashScanBitmaskQuery
        @test columnwise_stats.query_fold_ns > 0
        @test columnwise_stats.query_dedup_ns > 0
        @test columnwise_stats.query_insert_ns > 0
        @test bitmask_stats.query_insert_ns > 0

        many_guides = fill(first(guides_), 65)
        many_stats = CHOPOFF.PrefixHashScanStats()
        many_query = CHOPOFF.build_prefix_hash_scan_map_from_paths(UInt8[1 2; 2 3], many_guides, UInt8, many_stats; query_variant = :auto)
        @test many_stats.query_variant == :columnwise
        @test many_query isa Dict
        @test_throws ErrorException CHOPOFF.build_prefix_hash_scan_map_from_paths(UInt8[1 2; 2 3], many_guides, UInt8; query_variant = :bitmask64)
    end


    @testset "direct Cas9 candidate hash matches prefix helper" begin
        motif = Motif("Cas9"; distance = 2)
        hash_len = 16
        hash_type = phs_hash_type(hash_len)
        genome = joinpath(tdir, "cas9_direct_hash.fa")
        seq = repeat("A", 30) * "ACGTACGTACGTACGTACGT" * "AGG" * repeat("C", 20) * "CC" * "A" * "TGCATGCATGCATGCATGCA" * repeat("A", 30)
        write_phs_fasta(genome, "chr1", seq)
        chrom_seq = LongDNA{4}(seq)
        dbi = DBInfo(genome, "prefix_hash_scan_direct_hash", motif)
        @test CHOPOFF.is_cas9_prefix_hash_candidate(dbi, hash_len)
        for is_antisense in (false, true)
            positions = CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
            isempty(positions) && continue
            candidate_range = first(positions)
            prefix = CHOPOFF.normalized_candidate_prefix(chrom_seq, candidate_range, dbi, is_antisense, hash_len)
            old_hashes = CHOPOFF.candidate_prefix_hashes(prefix, hash_type, nothing)
            direct_hashes = CHOPOFF.candidate_prefix_hashes_direct_cas9(chrom_seq, candidate_range, is_antisense, hash_len, hash_type)
            @test direct_hashes == old_hashes
        end
    end

    @testset "prefix helper matches materialized candidate" begin
        check_prefix_helper_matches_materialized(tdir, repeat("A", 30) * "ACGTACGTACGTACGTACGT" * "AGG" * repeat("A", 30), Motif("Cas9"; distance = 2), 8, "cas9_helper")
        check_prefix_helper_matches_materialized(tdir, repeat("A", 30) * "TTTA" * "TGCATGCATGCATGCATGCAT" * repeat("A", 30), Motif("Cas12a"; distance = 2), 8, "cas12a_helper")
    end

    @testset "fused Cas9 scan and directory parity" begin
        motif3 = Motif("Cas9"; distance = 3)
        hash_len = 16
        hash_type = phs_hash_type(hash_len)
        guide3 = LongDNA{4}("ACGTACGTACGTACGTACGT")
        fused_genome = joinpath(tdir, "prefix_hash_scan_fused.fa")
        fused_seq = repeat("N", 5) * repeat("A", 30) *
            "ACGTACGTACGTACGTACGT" * "AGG" *
            repeat("C", 30) * "CC" * "A" *
            "TGCATGCATGCATGCATGCA" * repeat("A", 30) * repeat("N", 5)
        write_phs_fasta(fused_genome, "chr1", fused_seq)
        chrom_seq = LongDNA{4}(fused_seq)
        dbi = DBInfo(fused_genome, "prefix_hash_scan_fused", motif3)

        query_masks = Dict{UInt32, UInt64}()
        expected = [Tuple{Int, UInt64}[] for _ in 1:2]
        for (strand_idx, is_antisense) in enumerate((false, true))
            for candidate_range in CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
                hash = only(CHOPOFF.candidate_prefix_hashes_direct_cas9(
                    chrom_seq, candidate_range, is_antisense, hash_len, hash_type))
                mask = xor(UInt64(hash) * 0x9e3779b97f4a7c15, 0xd1b54a32d192ed03)
                mask == 0 && (mask = 1)
                query_masks[hash] = mask
                push!(expected[strand_idx], (first(candidate_range), mask))
            end
        end

        query = CHOPOFF.PrefixHashScanBitmaskQuery(query_masks, 1)
        plus_hits, minus_hits = CHOPOFF.scan_cas9_prefix_hits(chrom_seq, dbi, query, hash_len)
        @test [(hit.start, hit.mask) for hit in plus_hits] == expected[1]
        @test [(hit.start, hit.mask) for hit in minus_hits] == expected[2]

        directory = CHOPOFF.build_prefix_hash_scan_directory(query, hash_len, 8)
        for (hash, mask) in query_masks
            @test CHOPOFF.prefix_hash_scan_candidate_mask(directory, hash) == mask
        end
        @test CHOPOFF.prefix_hash_scan_candidate_mask(directory, UInt32(0x12345678)) ==
            get(query_masks, UInt32(0x12345678), UInt64(0))

        outputs = Dict{Symbol, DataFrame}()
        for backend in (:legacy, :fused_dict, :fused_directory)
            output = joinpath(tdir, "scan_" * string(backend) * ".csv")
            CHOPOFF.search_prefixHashScan(
                [guide3],
                fused_genome,
                motif3,
                output;
                distance = 3,
                hash_len = hash_len,
                early_stopping = fill(100, 4),
                query_variant = :bitmask64,
                scan_backend = backend,
                bucket_bases = 8,
            )
            df = DataFrame(CSV.File(output))
            sort!(df, names(df))
            outputs[backend] = df
        end
        @test nrow(outputs[:legacy]) > 0
        @test outputs[:fused_dict] == outputs[:legacy]
        @test outputs[:fused_directory] == outputs[:legacy]
    end

    guide = LongDNA{4}("ACGTACGTACGTACGTACGT")
    motif = Motif("Cas9"; distance = 2)
    genome = joinpath(tdir, "prefix_hash_scan.fa")
    pad = repeat("A", 40)
    exact = "ACGTACGTACGTACGTACGT" * "AGG"
    mismatch = "ACGTACGTACGTACGTACGA" * "TGG"
    deletion = "ACGTACGTAGTACGTACGT" * "AGG"
    insertion = "ACGTACGTAC" * "A" * "GTACGTACGT" * "AGG"
    write_phs_fasta(genome, "chr1", join([pad, exact, pad, mismatch, pad, deletion, pad, insertion, pad], ""))

    db_path = joinpath(tdir, "prefixHashDB")
    mkpath(db_path)
    build_prefixHashDB("prefix_hash_scan_test", genome, motif, db_path)

    prefix_out = joinpath(tdir, "prefixhash.csv")
    scan_out = joinpath(tdir, "scan.csv")
    brute_out = joinpath(tdir, "bruteforce.csv")
    search_prefixHashDB(db_path, [guide], prefix_out; distance = 2, early_stopping = fill(100, 3))
    CHOPOFF.search_prefixHashScan([guide], genome, motif, scan_out; distance = 2, early_stopping = fill(100, 3), query_variant = :bitmask64)
    CHOPOFF.search_prefixHashScan([guide], genome, motif, brute_out; distance = 2, early_stopping = fill(100, 3), query_variant = :bruteforce)

    @testset "parity vs prefixHashDB" begin
        @test phs_core(scan_out) == phs_core(prefix_out)
        @test phs_core(brute_out) == phs_core(prefix_out)
        @test phs_core(brute_out) == phs_core(scan_out)
    end

    @testset "indel candidates survive prefix scan" begin
        df = DataFrame(CSV.File(scan_out))
        @test any(==(1), df.distance)
        @test any(occursin.("-", df.alignment_reference))
        @test any(occursin.("-", df.alignment_guide))
    end
end
