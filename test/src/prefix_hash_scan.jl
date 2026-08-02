using Test
using CHOPOFF
using BioSequences
using FASTX
using TwoBit
using CSV
using DataFrames
using Logging

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

function write_phs_twobit(
    path::String, name::String, seq::String;
    masks::Union{Nothing, Vector{UnitRange{Int}}} = nothing)

    writer = TwoBit.Writer(open(path, "w"), [name])
    write(writer, TwoBit.Record(name, LongDNA{4}(seq), masks))
    close(writer)
end

function write_phs_swapped_twobit(path::String, name::String, seq::String)
    starts = UInt32[]
    sizes = UInt32[]
    idx = 1
    while idx <= length(seq)
        if seq[idx] == 'N'
            first_idx = idx
            while idx <= length(seq) && seq[idx] == 'N'
                idx += 1
            end
            push!(starts, UInt32(first_idx - 1))
            push!(sizes, UInt32(idx - first_idx))
        else
            idx += 1
        end
    end
    codes = Dict('T' => UInt8(0), 'N' => UInt8(0), 'C' => UInt8(1),
        'A' => UInt8(2), 'G' => UInt8(3))
    packed = UInt8[]
    for first_idx in 1:4:length(seq)
        byte = UInt8(0)
        for offset in 0:3
            idx = first_idx + offset
            code = idx <= length(seq) ? codes[seq[idx]] : UInt8(0)
            byte |= code << (6 - 2 * offset)
        end
        push!(packed, byte)
    end
    write_u32(io, value) = write(io, bswap(UInt32(value)))
    open(path, "w") do io
        write_u32(io, CHOPOFF.PREFIX_HASH_SCAN_TWOBIT_SIGNATURE)
        write_u32(io, 0)
        write_u32(io, 1)
        write_u32(io, 0)
        write(io, UInt8(ncodeunits(name)))
        write(io, name)
        write_u32(io, 16 + 1 + ncodeunits(name) + 4)
        write_u32(io, length(seq))
        write_u32(io, length(starts))
        foreach(value -> write_u32(io, value), starts)
        foreach(value -> write_u32(io, value), sizes)
        write_u32(io, 0)
        write_u32(io, 0)
        write(io, packed)
    end
end

@testset "prefixHashScan" begin
    tdir = tempname()
    mkpath(tdir)

    @testset "2bit range decoding" begin
        root = normpath(joinpath(dirname(pathof(CHOPOFF)), ".."))
        genome = joinpath(root, "test", "sample_data", "genome", "semirandom.2bit")
        index = CHOPOFF.read_prefix_hash_scan_twobit_index(genome)
        @test index.names == ["semirandom$i" for i in 1:8]
        @test index.lengths == fill(108326, 8)

        io = open(genome, "r")
        reader = TwoBit.Reader(io)
        buffer = UInt8[]
        sequence = TwoBit.sequence(reader[index.names[1]])
        for range in (1:1, 2:101, 63:191, 1021:4099,
                (length(sequence) - 130):length(sequence))
            observed = CHOPOFF.read_prefix_hash_scan_twobit_range!(
                buffer, io, index, 1, first(range), last(range))
            @test String(observed) == string(sequence[range])
        end
        close(io)

        sequence_string = "ACGTACNNNTGCATGCATGC"
        native = joinpath(tdir, "range_native.2bit")
        swapped = joinpath(tdir, "range_swapped.2bit")
        write_phs_twobit(native, "chr1", sequence_string; masks = [2:5])
        write_phs_swapped_twobit(swapped, "chr1", sequence_string)
        for path in (native, swapped)
            path_index = CHOPOFF.read_prefix_hash_scan_twobit_index(path)
            @test path_index.names == ["chr1"]
            @test path_index.lengths == [length(sequence_string)]
            open(path, "r") do path_io
                @test String(CHOPOFF.read_prefix_hash_scan_twobit_range!(
                    UInt8[], path_io, path_index, 1, 2, 19)) == sequence_string[2:19]
            end
        end
    end

    @testset "Cas9 distance-3 geometry" begin
        geometry = CHOPOFF.CAS9_D3_PREFIX_SCAN_GEOMETRY
        @test geometry.guide_bases == 20
        @test geometry.pam_bases == 3
        @test geometry.prefix_bases == 16
        @test geometry.distance == 3
        @test CHOPOFF.prefix_scan_candidate_bases(geometry) == 23
        @test CHOPOFF.prefix_scan_candidate_last_offset(geometry) == 22

        work, ranges = CHOPOFF.prefix_hash_scan_chunk_work([22, 23, 24], 8)
        @test ranges == [1:0, 1:1, 2:2]
        @test [(item.chrom_idx, item.core_first, item.core_last) for item in work] ==
            [(2, 1, 1), (3, 1, 2)]
    end

    @testset "production Cas9 distance geometries" begin
        for distance in 0:4
            resolved = CHOPOFF.resolve_prefix_scan_geometry(
                Motif("Cas9"; distance = distance), distance, 16)
            @test CHOPOFF.prefix_scan_kind(resolved) == :cas9
            @test resolved.guide_bases == 20
            @test resolved.pam_bases == 3
            @test resolved.prefix_bases == 16
            @test resolved.distance == distance
        end
        @test CHOPOFF.resolve_prefix_scan_geometry(
            Motif("Cas9"; distance = 5), 5, 16) === nothing
    end

    @testset "generic motif geometries and backends" begin
        guide20 = "ACGTTGCAACGTAGCTTGCA"
        guide23 = guide20 * "TGC"
        pamless = Motif(
            "pamless", repeat("N", 20), repeat("X", 20),
            true, true, 1, true, 0)
        internal = Motif(
            "internal", "NNNNNNNNNNXXXNNNNNNNNNN",
            "XXXXXXXXXXAGGXXXXXXXXXX", true, true, 1, false, 0)
        forward_only = Motif(
            "forward_only", repeat("N", 20) * "XXX",
            repeat("X", 20) * "AGG", true, false, 1, true, 0)
        cases = [
            (Motif("Cas9_NGA"; distance = 1), guide20, guide20 * "AGA"),
            (Motif("Cas9_NNAGAA"; distance = 1), guide20, guide20 * "AAAGAA"),
            (Motif("CasX"; distance = 1), guide20, "TTCA" * guide20),
            (Motif("hfCas12Max"; distance = 1), guide23, "TAC" * guide23),
            (pamless, guide20, guide20),
            (internal, guide20, guide20[1:10] * "AGG" * guide20[11:20]),
            (forward_only, guide20, guide20 * "AGG"),
        ]
        for (case_idx, (motif, guide_string, forward_site)) in enumerate(cases)
            geometry = CHOPOFF.resolve_prefix_scan_geometry(motif, 1, 16)
            @test CHOPOFF.prefix_scan_kind(geometry) == :generic
            reverse_site = string(reverse_complement(LongDNA{4}(forward_site)))
            sequence = repeat("A", 70) * forward_site * repeat("C", 70) *
                reverse_site * repeat("G", 70)
            genome = joinpath(tdir, "generic_motif_$case_idx.fa")
            write_phs_fasta(genome, "chr1", sequence)
            guide = LongDNA{4}(guide_string)
            outputs = String[]
            for backend in (:legacy, :fused_directory, :streaming_fasta_simd)
                output = joinpath(tdir, "generic_motif_$(case_idx)_$(backend).csv")
                CHOPOFF.search_prefixHashScan(
                    [guide], genome, motif, output;
                    distance = 1, hash_len = 16,
                    early_stopping = fill(100, 2), scan_backend = backend,
                    stream_chunk_bases = 64)
                push!(outputs, read(output, String))
            end
            @test all(==(first(outputs)), outputs)

            public_output = joinpath(tdir, "generic_motif_$(case_idx)_public.csv")
            search_prefixHashScan(
                [guide], genome, public_output; motif = motif, distance = 1,
                early_stopping = fill(100, 2))
            @test read(public_output, String) == first(outputs)
            if case_idx == 1
                twobit = joinpath(tdir, "generic_motif.2bit")
                write_phs_twobit(twobit, "chr1", sequence)
                twobit_output = joinpath(tdir, "generic_motif_2bit.csv")
                CHOPOFF.search_prefixHashScan(
                    [guide], twobit, motif, twobit_output;
                    distance = 1, hash_len = 16,
                    early_stopping = fill(100, 2),
                    scan_backend = :streaming_2bit_simd,
                    stream_chunk_bases = 64)
                @test read(twobit_output, String) == first(outputs)
            end
            if motif.alias == "forward_only"
                @test all(==("+"), DataFrame(CSV.File(public_output)).strand)
            end
        end

        for name in ("Cas9_NGA", "Cas9_NNAGAA", "CasX")
            stats = CHOPOFF.PrefixHashScanStats()
            CHOPOFF.load_prefix_hash_scan_paths(
                Motif(name; distance = 1), 1, 16, stats)
            @test stats.path_source == :precomputed
        end

        ambiguous_motif = Motif(
            "Cas9_NGA"; distance = 1, ambig_max = 1)
        ambiguous_site = guide20[1:7] * "R" * guide20[9:20] * "AGA"
        ambiguous_genome = joinpath(tdir, "generic_motif_ambiguous.fa")
        write_phs_fasta(
            ambiguous_genome, "chr1", repeat("A", 70) * ambiguous_site * repeat("C", 70))
        ambiguous_outputs = String[]
        for backend in (:legacy, :streaming_fasta_simd)
            output = joinpath(tdir, "generic_motif_ambiguous_$(backend).csv")
            CHOPOFF.search_prefixHashScan(
                [LongDNA{4}(guide20)], ambiguous_genome, ambiguous_motif, output;
                distance = 1, hash_len = 16,
                early_stopping = fill(100, 2), scan_backend = backend,
                stream_chunk_bases = 64)
            push!(ambiguous_outputs, read(output, String))
        end
        @test only(unique(ambiguous_outputs)) == first(ambiguous_outputs)

        short_motif = Motif("test"; distance = 2)
        @test CHOPOFF.resolve_prefix_scan_geometry(short_motif, 2, 1) === nothing
    end

    @testset "supported Cas9 API" begin
        guide = LongDNA{4}("ACGTACGTACGTACGTACGT")
        pad = repeat("A", 40)
        genome = joinpath(tdir, "supported_api.fa")
        twobit_genome = joinpath(tdir, "supported_api.2bit")
        four_edit = "CCTTCCTTACGTACGTACGT" * "TGG"
        supported_sequence =
            pad * string(guide) * "AGG" * pad * four_edit * pad
        write_phs_fasta(genome, "chr1", supported_sequence)
        write_phs_twobit(twobit_genome, "chr1", supported_sequence)
        public_out = joinpath(tdir, "supported_api.csv")
        twobit_out = joinpath(tdir, "supported_api_2bit.csv")
        verbose_out = joinpath(tdir, "supported_api_verbose.csv")
        engine_out = joinpath(tdir, "supported_api_engine.csv")

        @test_logs min_level=Logging.Info search_prefixHashScan(
            [guide], genome, public_out; early_stopping = fill(100, 4))
        @test_logs (:info, r"prefixHashScan execution") search_prefixHashScan(
            [guide], genome, verbose_out;
            early_stopping = fill(100, 4), verbose = true)
        motif = setdist(Motif("Cas9"), 3)
        CHOPOFF.search_prefixHashScan(
            [guide], genome, motif, engine_out;
            distance = 3, early_stopping = fill(100, 4))
        twobit_stats = CHOPOFF.PrefixHashScanStats()
        CHOPOFF.search_prefixHashScan(
            [guide], twobit_genome, motif, twobit_out;
            distance = 3, early_stopping = fill(100, 4), stats = twobit_stats)
        @test read(public_out) == read(verbose_out) == read(engine_out) == read(twobit_out)
        expected_twobit_backend = CHOPOFF.can_use_prefix_hash_scan_simd() ?
            :streaming_2bit_simd : :fused_directory
        @test twobit_stats.scan_backend == expected_twobit_backend


        for distance in 0:4
            motif_d = Motif("Cas9"; distance = distance)
            limits = fill(100, distance + 1)
            distance_outputs = String[]
            for backend in (
                    :legacy, :fused_directory, :streaming_fasta_simd, :auto)
                output = joinpath(
                    tdir, "supported_api_d$(distance)_$(backend).csv")
                stats = CHOPOFF.PrefixHashScanStats()
                CHOPOFF.search_prefixHashScan(
                    [guide], genome, motif_d, output;
                    distance = distance,
                    hash_len = 16,
                    early_stopping = limits,
                    scan_backend = backend,
                    stream_chunk_bases = 64,
                    stats = stats,
                )
                push!(distance_outputs, read(output, String))
                @test stats.path_rows == (1, 129, 7873, 302337, 8196801)[distance + 1]
                @test stats.path_source == :precomputed
            end
            @test all(==(first(distance_outputs)), distance_outputs)
            public_distance_out =
                joinpath(tdir, "supported_public_d$(distance).csv")
            search_prefixHashScan(
                [guide], genome, public_distance_out;
                distance = distance, early_stopping = limits)
            @test read(public_distance_out, String) == first(distance_outputs)
            twobit_distance_out =
                joinpath(tdir, "supported_public_2bit_d$(distance).csv")
            search_prefixHashScan(
                [guide], twobit_genome, twobit_distance_out;
                distance = distance, early_stopping = limits)
            @test read(twobit_distance_out, String) == first(distance_outputs)
            if distance == 4
                @test 4 in DataFrame(CSV.File(public_distance_out)).distance
                brute_out = joinpath(tdir, "supported_bruteforce_d4.csv")
                CHOPOFF.search_prefixHashScan(
                    [guide], genome, motif_d, brute_out;
                    distance = 4, hash_len = 16, early_stopping = limits,
                    query_variant = :bruteforce, scan_backend = :legacy)
                @test read(brute_out, String) == first(distance_outputs)
            end
        end
        ambiguous_genome = joinpath(tdir, "supported_api_ambiguous.fa")
        ambiguous_guide = "ACGTACGTACNTACGTACGT"
        write_phs_fasta(
            ambiguous_genome, "chr1", pad * ambiguous_guide * "AGG" * pad)
        ambiguous_out = joinpath(tdir, "supported_api_ambiguous.csv")
        search_prefixHashScan(
            [guide], ambiguous_genome, ambiguous_out;
            early_stopping = fill(100, 4))
        @test nrow(DataFrame(CSV.File(ambiguous_out))) == 0

        large_out = joinpath(tdir, "supported_api_large.csv")
        batch64_out = joinpath(tdir, "supported_api_batch64.csv")
        batch1_out = joinpath(tdir, "supported_api_batch1.csv")
        large_guides = fill(guide, 65)
        search_prefixHashScan(
            large_guides, genome, large_out;
            early_stopping = fill(100, 4))
        CHOPOFF.search_prefixHashScan(
            large_guides[1:64], genome, motif, batch64_out;
            distance = 3, early_stopping = fill(100, 4))
        CHOPOFF.search_prefixHashScan(
            large_guides[65:65], genome, motif, batch1_out;
            distance = 3, early_stopping = fill(100, 4))
        large = DataFrame(CSV.File(large_out))
        manual_batches = vcat(
            DataFrame(CSV.File(batch64_out)),
            DataFrame(CSV.File(batch1_out)))
        sort!(large, names(large))
        sort!(manual_batches, names(manual_batches))
        @test large == manual_batches
        @test count(
            ==("guide,alignment_guide,alignment_reference,distance,chromosome,start,strand"),
            readlines(large_out)) == 1

        large129_out = joinpath(tdir, "supported_api_large129.csv")
        search_prefixHashScan(
            fill(guide, 129), genome, large129_out;
            early_stopping = fill(1, 4))
        @test nrow(DataFrame(CSV.File(large129_out))) == 129

        @test_throws ErrorException search_prefixHashScan(
            LongDNA{4}[], genome, public_out)
        @test_throws ErrorException search_prefixHashScan(
            [LongDNA{4}("ACGTACGTACGTACGTACG")], genome, public_out)
        @test_throws ErrorException search_prefixHashScan(
            [LongDNA{4}("ACGTACGTACGTACGTACGN")], genome, public_out)
        @test_throws ErrorException search_prefixHashScan(
            [guide], genome, public_out; early_stopping = fill(100, 3))
        @test_throws SystemError search_prefixHashScan(
            [guide], joinpath(tdir, "reference.2bit"), public_out)

        unindexed = joinpath(tdir, "unindexed.fa")
        open(unindexed, "w") do io
            write(io, ">chr1\n", pad, string(guide), "AGG", pad, "\n")
        end
        @test_throws ErrorException search_prefixHashScan(
            [guide], unindexed, public_out)
    end

    @testset "specialized Cas12a geometries and backends" begin
        geometry = CHOPOFF.CAS12A_D3_PREFIX_SCAN_GEOMETRY
        @test CHOPOFF.prefix_scan_kind(geometry) == :cas12a
        @test geometry.guide_bases == 21
        @test geometry.pam_bases == 4
        @test geometry.prefix_bases == 16
        @test geometry.distance == 3
        @test CHOPOFF.prefix_scan_candidate_bases(geometry) == 25
        @test CHOPOFF.prefix_scan_candidate_last_offset(geometry) == 24

        work, ranges = CHOPOFF.prefix_hash_scan_chunk_work(
            [24, 25, 26], 8, geometry)
        @test ranges == [1:0, 1:1, 2:2]
        @test [(item.chrom_idx, item.core_first, item.core_last) for item in work] ==
            [(2, 1, 1), (3, 1, 2)]

        motif = Motif("Cas12a"; distance = 3)
        stats = CHOPOFF.PrefixHashScanStats()
        paths, source =
            CHOPOFF.load_prefix_hash_scan_paths(motif, 3, 16, stats)
        @test source == :precomputed
        @test stats.path_source == :precomputed
        @test size(paths) == (302337, 16)

        guide = LongDNA{4}("TGCATGCATGCATGCATGCAT")
        forward_site = "TTTA" * string(guide)
        reverse_site = string(reverse_complement(LongDNA{4}(forward_site)))
        four_edit_site = "TTTA" * "ATCCTTCATGCATGCATGCAT"
        seq = reverse_site * repeat("ACGT", 45) * forward_site *
            repeat("ACGT", 20) * four_edit_site

        for distance in 0:4
            resolved = CHOPOFF.resolve_prefix_scan_geometry(
                Motif("Cas12a"; distance = distance), distance, 16)
            @test CHOPOFF.prefix_scan_kind(resolved) == :cas12a
            @test resolved.guide_bases == 21
            @test resolved.pam_bases == 4
            @test resolved.prefix_bases == 16
            @test resolved.distance == distance
        end
        @test CHOPOFF.resolve_prefix_scan_geometry(
            Motif("Cas12a"; distance = 5), 5, 16) === nothing
        genome = joinpath(tdir, "cas12a_specialized.fa")
        twobit_genome = joinpath(tdir, "cas12a_specialized.2bit")
        write_phs_fasta(genome, "chr1", seq)
        write_phs_twobit(twobit_genome, "chr1", seq)
        chrom_seq = LongDNA{4}(seq)
        dbi = DBInfo(genome, "cas12a_specialized", motif)

        for is_antisense in (false, true)
            for candidate_range in CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
                prefix = CHOPOFF.normalized_candidate_prefix(
                    chrom_seq, candidate_range, dbi, is_antisense, 16)
                expected_hash = CHOPOFF.candidate_prefix_hashes(
                    prefix, UInt32, nothing)
                direct_hash = CHOPOFF.candidate_prefix_hashes_direct_cas12a(
                    chrom_seq, candidate_range, is_antisense, 16, UInt32)
                @test direct_hash == expected_hash

                expected_materialized = CHOPOFF.materialize_normalized_candidate(
                    chrom_seq, candidate_range, dbi, is_antisense)
                observed = CHOPOFF.materialize_normalized_candidate_cas12a(
                    chrom_seq, first(candidate_range), dbi, is_antisense)
                raw_observed = CHOPOFF.materialize_normalized_candidate_cas12a(
                    collect(codeunits(seq)), first(candidate_range), dbi, is_antisense)
                @test observed == expected_materialized
                @test raw_observed == expected_materialized
            end
        end

        outputs = Dict{Symbol, String}()
        backend_stats = Dict{Symbol, CHOPOFF.PrefixHashScanStats}()
        for backend in (
                :legacy, :fused_dict, :fused_directory, :fused_fasta_simd,
                :streaming_fasta_simd, :streaming_fasta_simd_fused, :auto)
            output = joinpath(tdir, "cas12a_$(backend).csv")
            backend_stats[backend] = CHOPOFF.PrefixHashScanStats()
            CHOPOFF.search_prefixHashScan(
                [guide], genome, motif, output;
                distance = 3,
                hash_len = 16,
                early_stopping = fill(100, 4),
                scan_backend = backend,
                stats = backend_stats[backend],
            )
            outputs[backend] = read(output, String)
            @test backend_stats[backend].path_source == :precomputed
        end
        @test all(==(outputs[:legacy]), values(outputs))
        expected_auto = CHOPOFF.can_use_prefix_hash_scan_simd() ?
            :streaming_fasta_simd : :fused_directory
        @test backend_stats[:auto].scan_backend == expected_auto

        public_output = joinpath(tdir, "cas12a_public.csv")
        search_prefixHashScan(
            [guide], genome, public_output;
            motif = "Cas12a", early_stopping = fill(100, 4))
        @test read(public_output, String) == outputs[:legacy]
        twobit_output = joinpath(tdir, "cas12a_public_2bit.csv")
        search_prefixHashScan(
            [guide], twobit_genome, twobit_output;
            motif = "Cas12a", early_stopping = fill(100, 4))
        @test read(twobit_output, String) == outputs[:legacy]

        large_public_output = joinpath(tdir, "cas12a_public_large.csv")
        search_prefixHashScan(
            fill(guide, 65), genome, large_public_output;
            motif = "Cas12a", early_stopping = fill(100, 4))
        @test nrow(DataFrame(CSV.File(large_public_output))) ==
            65 * nrow(DataFrame(CSV.File(public_output)))
        @test count(
            ==("guide,alignment_guide,alignment_reference,distance,chromosome,start,strand"),
            readlines(large_public_output)) == 1

        invalid_genome = joinpath(tdir, "cas12a_invalid.fa")
        write_phs_fasta(
            invalid_genome, "chr1", "TTTT" * string(guide))

        for distance in (0, 1, 2, 4)
            motif_d = Motif("Cas12a"; distance = distance)
            limits = fill(100, distance + 1)
            distance_outputs = String[]
            for backend in (
                    :legacy, :fused_directory, :streaming_fasta_simd, :auto)
                output = joinpath(
                    tdir, "cas12a_d$(distance)_$(backend).csv")
                stats_d = CHOPOFF.PrefixHashScanStats()
                CHOPOFF.search_prefixHashScan(
                    [guide], genome, motif_d, output;
                    distance = distance,
                    hash_len = 16,
                    early_stopping = limits,
                    scan_backend = backend,
                    stream_chunk_bases = 64,
                    stats = stats_d,
                )
                push!(distance_outputs, read(output, String))
                @test stats_d.path_rows ==
                    (1, 129, 7873, 302337, 8196801)[distance + 1]
                @test stats_d.path_source == :precomputed
            end
            @test all(==(first(distance_outputs)), distance_outputs)
            public_distance_out =
                joinpath(tdir, "cas12a_public_d$(distance).csv")
            search_prefixHashScan(
                [guide], genome, public_distance_out;
                motif = "Cas12a", distance = distance,
                early_stopping = limits)
            @test read(public_distance_out, String) == first(distance_outputs)
            twobit_distance_out =
                joinpath(tdir, "cas12a_public_2bit_d$(distance).csv")
            search_prefixHashScan(
                [guide], twobit_genome, twobit_distance_out;
                motif = "Cas12a", distance = distance,
                early_stopping = limits)
            @test read(twobit_distance_out, String) == first(distance_outputs)
            if distance == 4
                @test 4 in DataFrame(CSV.File(public_distance_out)).distance
                brute_out = joinpath(tdir, "cas12a_bruteforce_d4.csv")
                CHOPOFF.search_prefixHashScan(
                    [guide], genome, motif_d, brute_out;
                    distance = 4, hash_len = 16, early_stopping = limits,
                    query_variant = :bruteforce, scan_backend = :legacy)
                @test read(brute_out, String) == first(distance_outputs)
            end
        end
        invalid_output = joinpath(tdir, "cas12a_invalid.csv")
        search_prefixHashScan(
            [guide], invalid_genome, invalid_output;
            motif = "Cas12a", early_stopping = fill(100, 4))
        @test nrow(DataFrame(CSV.File(invalid_output))) == 0

        ambiguous_genome = joinpath(tdir, "cas12a_ambiguous.fa")
        guide_string = string(guide)
        ambiguous_site =
            "TTTA" * guide_string[1:7] * "N" * guide_string[9:end]
        write_phs_fasta(
            ambiguous_genome, "chr1", ambiguous_site)
        ambiguous_output = joinpath(tdir, "cas12a_ambiguous.csv")
        search_prefixHashScan(
            [guide], ambiguous_genome, ambiguous_output;
            motif = "Cas12a", early_stopping = fill(100, 4))
        @test nrow(DataFrame(CSV.File(ambiguous_output))) == 0
    end

    @testset "reference ambiguity" begin
        function run_ambiguous_site(
            label::String, guide::LongDNA{4}, site::String, motif::Motif;
            backend::Symbol = :auto,
            twobit::Bool = false)

            genome = joinpath(tdir, label * (twobit ? ".2bit" : ".fa"))
            output = joinpath(tdir, label * "_" * string(backend) * ".csv")
            sequence = repeat("A", 40) * site * repeat("A", 40)
            if twobit
                write_phs_twobit(genome, "chr1", sequence)
            else
                write_phs_fasta(genome, "chr1", sequence)
            end
            stats = CHOPOFF.PrefixHashScanStats()
            CHOPOFF.search_prefixHashScan(
                [guide], genome, motif, output;
                distance = motif.distance,
                hash_len = 16,
                early_stopping = fill(100, motif.distance + 1),
                scan_backend = backend,
                stream_chunk_bases = 64,
                stats = stats,
            )
            return DataFrame(CSV.File(output)), read(output, String), stats
        end

        cas9_guide = LongDNA{4}("ACGTACGTACGTACGTACGT")
        cas9_string = string(cas9_guide)
        one_n = cas9_string[1:9] * "N" * cas9_string[11:end] * "AGG"
        motif0 = Motif("Cas9"; distance = 0, ambig_max = 0)
        motif1 = Motif("Cas9"; distance = 0, ambig_max = 1)
        rejected, _, _ = run_ambiguous_site(
            "cas9_ambig_rejected", cas9_guide, one_n, motif0)
        @test nrow(rejected) == 0
        rejected_twobit, _, _ = run_ambiguous_site(
            "cas9_ambig_rejected_2bit", cas9_guide, one_n, motif0;
            twobit = true)
        accepted_twobit, _, accepted_twobit_stats = run_ambiguous_site(
            "cas9_ambig_accepted_2bit", cas9_guide, one_n, motif1;
            twobit = true)
        @test nrow(rejected_twobit) == 0
        @test nrow(accepted_twobit) == 1
        @test accepted_twobit_stats.ambiguous_prefixes == 1

        backend_outputs = String[]
        for backend in (
                :legacy, :fused_dict, :fused_directory, :fused_fasta_simd,
                :streaming_fasta_simd, :streaming_fasta_simd_fused, :auto)
            accepted, output, stats = run_ambiguous_site(
                "cas9_ambig_backend", cas9_guide, one_n, motif1;
                backend = backend)
            @test nrow(accepted) == 1
            @test stats.ambiguous_prefixes == 1
            push!(backend_outputs, output)
        end
        @test all(==(first(backend_outputs)), backend_outputs)

        three_n = collect(cas9_string)
        three_n[[6, 10, 14]] .= 'N'
        accepted3, _, _ = run_ambiguous_site(
            "cas9_ambig_three", cas9_guide, String(three_n) * "AGG",
            Motif("Cas9"; distance = 0, ambig_max = 3))
        @test nrow(accepted3) == 1

        distance3_site = collect(cas9_string * "AGG")
        distance3_site[6] = 'Y'
        distance3_site[18] = 'N'
        distance3_site[22] = 'R'
        distance3_outputs = String[]
        for backend in (:legacy, :streaming_fasta_simd, :auto)
            result, output, stats = run_ambiguous_site(
                "cas9_d3_ambig_three", cas9_guide, String(distance3_site),
                Motif("Cas9"; distance = 3, ambig_max = 3);
                backend = backend)
            @test count(==(0), result.distance) == 1
            @test stats.ambiguous_prefixes >= 1
            push!(distance3_outputs, output)
        end
        @test all(==(first(distance3_outputs)), distance3_outputs)

        four_n = copy(three_n)
        four_n[18] = 'N'
        rejected4, _, _ = run_ambiguous_site(
            "cas9_ambig_four", cas9_guide, String(four_n) * "AGG",
            Motif("Cas9"; distance = 0, ambig_max = 3))
        @test nrow(rejected4) == 0

        compatible_pam, _, _ = run_ambiguous_site(
            "cas9_ambig_pam", cas9_guide, cas9_string * "ARG", motif1)
        incompatible_pam, _, _ = run_ambiguous_site(
            "cas9_ambig_bad_pam", cas9_guide, cas9_string * "AYG", motif1)
        @test nrow(compatible_pam) == 1
        @test nrow(incompatible_pam) == 0
        invalid_base, _, _ = run_ambiguous_site(
            "cas9_invalid_reference", cas9_guide,
            cas9_string[1:9] * "X" * cas9_string[11:end] * "AGG", motif1)
        @test nrow(invalid_base) == 0

        compatible_base = Dict(
            'R' => 'A', 'Y' => 'C', 'S' => 'C', 'W' => 'A', 'K' => 'G',
            'M' => 'A', 'B' => 'C', 'D' => 'A', 'H' => 'A', 'V' => 'A',
            'N' => 'A', 'n' => 'A')
        for (code, base) in compatible_base
            site = collect(cas9_string)
            site[findfirst(==(base), site)] = code
            result, _, _ = run_ambiguous_site(
                "cas9_iupac_" * string(code), cas9_guide,
                String(site) * "AGG", motif1)
            @test nrow(result) == 1
        end

        cas12a_guide = LongDNA{4}("TGCATGCATGCATGCATGCAT")
        cas12a_string = string(cas12a_guide)
        cas12a_site = "TTTR" * cas12a_string[1:7] * "N" * cas12a_string[9:end]
        cas12a_motif = Motif("Cas12a"; distance = 0, ambig_max = 2)
        cas12a_outputs = String[]
        for backend in (:legacy, :fused_directory, :streaming_fasta_simd, :auto)
            result, output, _ = run_ambiguous_site(
                "cas12a_ambig_backend", cas12a_guide, cas12a_site,
                cas12a_motif; backend = backend)
            @test nrow(result) == 2
            push!(cas12a_outputs, output)
        end
        @test all(==(first(cas12a_outputs)), cas12a_outputs)

        long_n, _, _ = run_ambiguous_site(
            "cas9_long_n", cas9_guide, repeat("N", 200),
            Motif("Cas9"; distance = 0, ambig_max = 3))
        @test nrow(long_n) == 0
        @test_throws ErrorException run_ambiguous_site(
            "cas9_ambig_too_high", cas9_guide, one_n,
            Motif("Cas9"; distance = 0, ambig_max = 4))
    end


    @testset "precomputed distance assets and Cas9 fallback" begin

    @testset "exact production distance assets" begin
        expected_rows = (1, 129, 7873, 302337, 8196801)
        for motif_name in ("Cas9", "Cas12a"), distance in 0:4
            motif = Motif(motif_name; distance = distance)
            stats = CHOPOFF.PrefixHashScanStats()
            paths, source = CHOPOFF.load_prefix_hash_scan_paths(
                motif, distance, 16, stats)
            @test source == :precomputed
            @test stats.path_source == :precomputed
            @test stats.path_rows == expected_rows[distance + 1]
            @test size(paths) == (expected_rows[distance + 1], 16)

            loaded = CHOPOFF.load_precomputed_prefix_paths(
                motif, distance, 16; need_distances = true)
            @test loaded !== nothing
            @test size(loaded[1]) == size(paths)
            @test length(loaded[2]) == expected_rows[distance + 1]
        end
    end
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

        compact_guides = [guides; first(guides)]
        serial_compact_stats = CHOPOFF.PrefixHashScanStats()
        parallel_compact_stats = CHOPOFF.PrefixHashScanStats()
        serial_compact, serial_guides = CHOPOFF.build_prefix_hash_scan_compact_query(
            compact_guides, motif, 3, hash_len, serial_compact_stats;
            bucket_bases = 11, prefilter_bits = 26,
            query_build_backend = :serial, query_threads = 4)
        parallel_compact, parallel_guides = CHOPOFF.build_prefix_hash_scan_compact_query(
            compact_guides, motif, 3, hash_len, parallel_compact_stats;
            bucket_bases = 11, prefilter_bits = 26,
            query_build_backend = :parallel, query_threads = 4)
        @test parallel_guides == serial_guides
        @test parallel_compact.presence == serial_compact.presence
        @test parallel_compact.prefix_bits == serial_compact.prefix_bits
        @test parallel_compact.directory.offsets == serial_compact.directory.offsets
        @test parallel_compact.directory.suffixes == serial_compact.directory.suffixes
        @test parallel_compact.directory.masks == serial_compact.directory.masks
        @test parallel_compact_stats.path_rows == serial_compact_stats.path_rows
        @test parallel_compact_stats.query_hashes == serial_compact_stats.query_hashes
        @test parallel_compact_stats.query_variant == :bitmask64
        @test serial_compact_stats.query_variant == :bitmask64
        one_serial, _ = CHOPOFF.build_prefix_hash_scan_compact_query(
            guides[1:1], motif, 3, hash_len, nothing; bucket_bases = 11,
            prefilter_bits = 26, query_build_backend = :serial, query_threads = 4)
        one_auto, _ = CHOPOFF.build_prefix_hash_scan_compact_query(
            guides[1:1], motif, 3, hash_len, nothing; bucket_bases = 11,
            prefilter_bits = 26, query_build_backend = :auto, query_threads = 4)
        @test one_auto.presence == one_serial.presence
        @test one_auto.directory.masks == one_serial.directory.masks

        @test_throws ErrorException CHOPOFF.build_prefix_hash_scan_compact_query(
            guides, motif, 3, hash_len, nothing; bucket_bases = 11,
            prefilter_bits = 26, query_build_backend = :invalid)
        @test_throws ErrorException CHOPOFF.build_prefix_hash_scan_compact_query(
            guides, motif, 3, hash_len, nothing; bucket_bases = 11,
            prefilter_bits = 26, query_threads = 0)

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

    @testset "direct Cas9 materialization matches generic boundaries" begin
        motif3 = Motif("Cas9"; distance = 3)
        guide3 = "ACGTACGTACGTACGTACGT"
        sequences = [
            guide3 * "AGG",
            "A" * guide3 * "AGG",
            "AA" * guide3 * "AGG",
            "CCA" * guide3,
            "A" * "CCA" * guide3,
            "AA" * "CCA" * guide3,
        ]
        for (seq_idx, seq) in enumerate(sequences)
            genome = joinpath(tdir, "cas9_materialize_boundary_$(seq_idx).fa")
            write_phs_fasta(genome, "chr1", seq)
            chrom_seq = LongDNA{4}(seq)
            dbi = DBInfo(genome, "cas9_materialize_boundary", motif3)
            for is_antisense in (false, true)
                for candidate_range in CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
                    expected = CHOPOFF.materialize_normalized_candidate(
                        chrom_seq, candidate_range, dbi, is_antisense)
                    observed = CHOPOFF.materialize_normalized_candidate_cas9(
                        chrom_seq, first(candidate_range), dbi, is_antisense)
                    raw_observed = CHOPOFF.materialize_normalized_candidate_cas9(
                        collect(codeunits(seq)), first(candidate_range), dbi, is_antisense)
                    @test observed == expected
                    @test raw_observed == expected
                end
            end
        end
    end

    @testset "raw SIMD block and ambiguity parity" begin
        motif3 = Motif("Cas9"; distance = 3)
        hash_len = 16
        target = collect(codeunits("ACGTACGTACGTACGTACGTAGG"))
        lower_target = collect(codeunits("acgtacgtacgtacgtacgtagg"))
        raw = fill(UInt8('A'), 340)
        for (idx, start) in enumerate((1, 31, 64, 127, 190, 250, 315))
            bytes = isodd(idx) ? target : lower_target
            copyto!(raw, start, bytes, 1, length(bytes))
        end
        raw[120] = UInt8('R')
        raw[180] = UInt8('N')
        seq = String(copy(raw))
        genome = joinpath(tdir, "prefix_hash_scan_simd_blocks.fa")
        write_phs_fasta(genome, "chr1", seq)
        chrom_seq = LongDNA{4}(seq)
        dbi = DBInfo(genome, "prefix_hash_scan_simd_blocks", motif3)

        masks = Dict{UInt32, UInt64}()
        for is_antisense in (false, true)
            for candidate_range in CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
                hash = only(CHOPOFF.candidate_prefix_hashes_direct_cas9(
                    chrom_seq, candidate_range, is_antisense, hash_len, UInt32))
                masks[hash] = get(masks, hash, UInt64(0)) | UInt64(1)
            end
        end
        query = CHOPOFF.PrefixHashScanBitmaskQuery(masks, 1)
        directory = CHOPOFF.build_prefix_hash_scan_directory(query, hash_len, 8)
        expected = CHOPOFF.scan_cas9_prefix_hits(
            chrom_seq, dbi, directory, hash_len; scan_threads = 1)
        for bits in (22, 24, 26)
            filtered = CHOPOFF.build_prefix_hash_scan_prefilter(
                directory, keys(masks), bits)
            for threads in (1, 4)
                observed = CHOPOFF.scan_cas9_prefix_hits_raw(
                    raw, dbi, filtered; scan_threads = threads)
                @test observed == expected
            end
        end
    end

    @testset "raw Myers boundary and IUPAC parity" begin
        motif3 = Motif("Cas9"; distance = 3)
        guide = LongDNA{4}("ACGTACGTACGTACGTACGT")
        guide_oriented = reverse(guide)
        profile = CHOPOFF.build_prefix_hash_scan_myers_profile(guide_oriented)
        fixture_genome = joinpath(tdir, "prefix_hash_scan_simd_blocks.fa")

        plus_raw = collect(codeunits(String(guide) * "AGG" * "RYN"))
        plus_dbi = DBInfo(fixture_genome, "raw_myers_plus", motif3)
        plus_ot, _ = CHOPOFF.materialize_normalized_candidate_cas9(
            plus_raw, 1, plus_dbi, false)
        @test CHOPOFF.prefix_hash_scan_raw_myers_distance(
            profile, plus_raw, 1, false, 3) ==
            CHOPOFF.levenshtein(guide_oriented, plus_ot, 3, iscompatible)

        minus_guide = String(reverse(complement(guide)))
        minus_raw = collect(codeunits("RYNCCN" * minus_guide))
        minus_start = 4
        minus_dbi = DBInfo(fixture_genome, "raw_myers_minus", motif3)
        minus_ot, _ = CHOPOFF.materialize_normalized_candidate_cas9(
            minus_raw, minus_start, minus_dbi, true)
        @test CHOPOFF.prefix_hash_scan_raw_myers_distance(
            profile, minus_raw, minus_start, true, 3) ==
            CHOPOFF.levenshtein(guide_oriented, minus_ot, 3, iscompatible)
    end

    @testset "FAI search metadata and streamed validation" begin
        motif3 = Motif("Cas9"; distance = 3)
        genome = joinpath(tdir, "prefix_hash_scan_fai.fa")
        seq = repeat("ACGT", 20)
        write_phs_fasta(genome, "chr1", seq)
        dbi, lengths = CHOPOFF.prefix_hash_scan_dbinfo(genome, motif3)
        @test dbi.gi.genomechecksum == 0
        @test dbi.gi.chrom == ["chr1"]
        @test lengths == [length(seq)]

        open(genome) do io
            records = CHOPOFF.PrefixHashScanFASTARecords(
                FASTA.Reader(io; copy = false), ["chr1"], [length(seq) + 1])
            @test_throws ErrorException first(records)
        end

        missing_fai = joinpath(tdir, "prefix_hash_scan_missing_fai.fa")
        open(missing_fai, "w") do io
            write(io, ">chr1
", seq, "
")
        end
        @test_throws ErrorException CHOPOFF.prefix_hash_scan_dbinfo(
            missing_fai, motif3)

        crlf_genome = joinpath(tdir, "prefix_hash_scan_crlf.fa")
        open(crlf_genome, "w") do io
            write(io, ">chr1\r\nACGT\r\nTGCA\r\nAC\r\n")
        end
        open(crlf_genome * ".fai", "w") do io
            write(io, "chr1\t10\t7\t4\t6\n")
        end
        crlf_index = FASTA.Index(crlf_genome * ".fai")
        open(crlf_genome) do io
            buffer = UInt8[]
            @test String(CHOPOFF.read_prefix_hash_scan_fasta_range!(
                buffer, io, crlf_index, 1, 1, 10)) == "ACGTTGCAAC"
            @test String(CHOPOFF.read_prefix_hash_scan_fasta_range!(
                buffer, io, crlf_index, 1, 3, 9)) == "GTTGCAA"
            @test String(CHOPOFF.read_prefix_hash_scan_fasta_range!(
                buffer, io, crlf_index, 1, 4, 6)) == "TTG"
        end
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
        raw = collect(codeunits(fused_seq))
        raw_plus, raw_minus = CHOPOFF.scan_cas9_prefix_hits_raw(
            raw, dbi, directory; scan_threads = 1)
        threaded_plus, threaded_minus = CHOPOFF.scan_cas9_prefix_hits_raw(
            raw, dbi, directory; scan_threads = 4)
        @test [(hit.start, hit.mask) for hit in raw_plus] == expected[1]
        @test [(hit.start, hit.mask) for hit in raw_minus] == expected[2]
        @test raw_plus == threaded_plus
        @test raw_minus == threaded_minus

        scratch_plus = CHOPOFF.PrefixHashScanHit[]
        scratch_minus = CHOPOFF.PrefixHashScanHit[]
        scratch_plus_id = objectid(scratch_plus)
        scratch_minus_id = objectid(scratch_minus)
        bounds = CHOPOFF.cas9_prefix_scan_bounds_raw(raw, dbi)
        motif_candidates = CHOPOFF.scan_cas9_prefix_hits_raw_range!(
            scratch_plus, scratch_minus, raw, directory, bounds...)
        @test scratch_plus == raw_plus
        @test scratch_minus == raw_minus
        @test motif_candidates >= length(scratch_plus) + length(scratch_minus)
        @test objectid(scratch_plus) == scratch_plus_id
        @test objectid(scratch_minus) == scratch_minus_id

        CHOPOFF.scan_cas9_prefix_hits_raw_range!(
            scratch_plus, scratch_minus, raw, directory,
            2, 1, 2, 1, 2, 1)
        @test isempty(scratch_plus)
        @test isempty(scratch_minus)
        @test objectid(scratch_plus) == scratch_plus_id
        @test objectid(scratch_minus) == scratch_minus_id

        bucket_directory = CHOPOFF.build_prefix_hash_scan_directory(
            query, hash_len, 11)
        bucket_query = CHOPOFF.build_prefix_hash_scan_prefilter(
            bucket_directory, collect(keys(query_masks)), 26)
        bucket_plus = CHOPOFF.PrefixHashScanHit[]
        bucket_minus = CHOPOFF.PrefixHashScanHit[]
        lookup_scratch = CHOPOFF.PrefixHashScanLookupScratch()
        bucket_motif_candidates =
            CHOPOFF.scan_cas9_prefix_hits_raw_range_bucketed!(
                bucket_plus, bucket_minus, lookup_scratch.plus_candidates,
                lookup_scratch.minus_candidates, lookup_scratch.plus_radix,
                lookup_scratch.minus_radix, lookup_scratch.radix_counts, raw,
                bucket_query, bounds...)
        @test bucket_plus == raw_plus
        @test bucket_minus == raw_minus
        @test bucket_motif_candidates == motif_candidates
        @test issorted(lookup_scratch.plus_radix)
        @test issorted(lookup_scratch.minus_radix)
        @test all(hash -> CHOPOFF.prefix_hash_scan_prefilter_contains(
            bucket_query, hash), keys(query_masks))

        duplicate_hash = first(keys(query_masks))
        duplicate_mask = query_masks[duplicate_hash]
        duplicate_candidates = [
            (UInt64(duplicate_hash) << 32) | UInt64(start)
            for start in (9, 3, 7)]
        duplicate_hits = CHOPOFF.PrefixHashScanHit[]
        duplicate_radix = UInt64[]
        CHOPOFF.resolve_prefix_hash_scan_bucketed_hits!(
            duplicate_hits, duplicate_candidates, duplicate_radix,
            lookup_scratch.radix_counts, bucket_query)
        @test [(hit.start, hit.mask) for hit in duplicate_hits] ==
            [(3, duplicate_mask), (7, duplicate_mask), (9, duplicate_mask)]

        guide_oriented = reverse(guide3)
        myers_profile = CHOPOFF.build_prefix_hash_scan_myers_profile(guide_oriented)
        for is_antisense in (false, true)
            for candidate_range in CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
                ot, _ = CHOPOFF.materialize_normalized_candidate_cas9(
                    raw, first(candidate_range), dbi, is_antisense)
                expected_distance = CHOPOFF.levenshtein(
                    guide_oriented, ot, 3, iscompatible)
                observed_distance = CHOPOFF.prefix_hash_scan_raw_myers_distance(
                    myers_profile, raw, first(candidate_range), is_antisense, 3)
                @test observed_distance == expected_distance
            end
        end
        for (hash, mask) in query_masks
            @test CHOPOFF.prefix_hash_scan_candidate_mask(directory, hash) == mask
        end
        @test CHOPOFF.prefix_hash_scan_candidate_mask(directory, UInt32(0x12345678)) ==
            get(query_masks, UInt32(0x12345678), UInt64(0))

        work, work_ranges = CHOPOFF.prefix_hash_scan_chunk_work(
            [0, 22, 23, 64, 65], 32)
        @test isempty(work_ranges[1])
        @test isempty(work_ranges[2])
        @test [(item.chrom_idx, item.core_first, item.core_last) for item in work] == [
            (3, 1, 1),
            (4, 1, 32),
            (4, 33, 42),
            (5, 1, 32),
            (5, 33, 43),
        ]

        stream_query = CHOPOFF.PrefixHashScanBitmaskQuery(
            Dict(hash => UInt64(1) for hash in keys(query_masks)), 1)
        stream_directory = CHOPOFF.build_prefix_hash_scan_directory(
            stream_query, hash_len, 8)
        reference_lengths = FASTA.Index(fused_genome * ".fai").lengths
        function stream_result_tuples(results, ranges)
            tuples = Tuple[]
            for chrom_idx in eachindex(ranges), strand in (:plus, :minus)
                for result_idx in ranges[chrom_idx]
                    for hit in getfield(something(results[result_idx]), strand)
                        push!(tuples, (
                            chrom_idx, strand, hit.guide_idx, hit.pos, hit.dist,
                            hit.is_antisense, hit.aln_guide, hit.aln_ref))
                    end
                end
            end
            return tuples
        end
        scheduler_stat_fields = (
            :motif_candidates, :prefix_hits, :guide_pairs,
            :alignment_calls, :distance_calls, :traceback_calls,
        )
        for mode in (:buffered_reuse, :fused), with_stats in (false, true)
            chunk_stats = with_stats ? CHOPOFF.PrefixHashScanStats() : nothing
            chrom_stats = with_stats ? CHOPOFF.PrefixHashScanStats() : nothing
            chunk_results, chunk_ranges = CHOPOFF.stream_prefix_hash_scan(
                fused_genome, reference_lengths, stream_directory, dbi,
                [guide_oriented], [myers_profile], 3, 64, 4, Val(mode),
                chunk_stats, Val(:chunk))
            chrom_results, chrom_ranges = CHOPOFF.stream_prefix_hash_scan(
                fused_genome, reference_lengths, stream_directory, dbi,
                [guide_oriented], [myers_profile], 3, 64, 4, Val(mode),
                chrom_stats, Val(:chromosome))
            @test stream_result_tuples(chunk_results, chunk_ranges) ==
                stream_result_tuples(chrom_results, chrom_ranges)
            if with_stats
                chunk_total = CHOPOFF.PrefixHashScanStats()
                chrom_total = CHOPOFF.PrefixHashScanStats()
                for result in chunk_results
                    CHOPOFF.merge_prefix_hash_scan_worker_stats!(
                        chunk_total, something(result).stats)
                end
                for result in chrom_results
                    CHOPOFF.merge_prefix_hash_scan_worker_stats!(
                        chrom_total, something(result).stats)
                end
                for field in scheduler_stat_fields
                    @test getfield(chunk_total, field) == getfield(chrom_total, field)
                end
            end
        end

        outputs = Dict{Symbol, DataFrame}()
        for backend in (
                :legacy, :fused_dict, :fused_directory, :fused_fasta_simd,
                :streaming_fasta_simd, :streaming_fasta_simd_fused, :auto)
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
                stream_chunk_bases = 64,
            )
            df = DataFrame(CSV.File(output))
            sort!(df, names(df))
            outputs[backend] = df
        end
        @test nrow(outputs[:legacy]) > 0
        @test outputs[:fused_dict] == outputs[:legacy]
        @test outputs[:fused_directory] == outputs[:legacy]
        @test outputs[:fused_fasta_simd] == outputs[:legacy]
        @test outputs[:streaming_fasta_simd] == outputs[:legacy]
        @test outputs[:streaming_fasta_simd_fused] == outputs[:legacy]
        @test read(joinpath(tdir, "scan_streaming_fasta_simd.csv")) ==
            read(joinpath(tdir, "scan_streaming_fasta_simd_fused.csv"))
        @test outputs[:auto] == outputs[:legacy]

        default_output = joinpath(tdir, "scan_default_tuning.csv")
        explicit_output = joinpath(tdir, "scan_explicit_tuning.csv")
        default_kwargs = (
            distance = 3,
            early_stopping = fill(100, 4),
            query_variant = :bitmask64,
            scan_backend = :auto,
        )
        CHOPOFF.search_prefixHashScan(
            [guide3], fused_genome, motif3, default_output; default_kwargs...)
        CHOPOFF.search_prefixHashScan(
            [guide3], fused_genome, motif3, explicit_output;
            default_kwargs...,
            bucket_bases = 11,
            stream_chunk_bases = 2 * 1024 * 1024,
            prefilter_bits = 26,
        )
        @test read(default_output) == read(explicit_output)
        lookup_outputs = Dict{Symbol, Vector{UInt8}}()
        lookup_stats = Dict{Symbol, CHOPOFF.PrefixHashScanStats}()
        for lookup_variant in (:inline, :bucketed, :auto)
            output = joinpath(
                tdir, "scan_lookup_" * string(lookup_variant) * ".csv")
            variant_stats = CHOPOFF.PrefixHashScanStats()
            CHOPOFF.search_prefixHashScan(
                [guide3], fused_genome, motif3, output;
                distance = 3, hash_len = hash_len,
                early_stopping = fill(100, 4), query_variant = :bitmask64,
                scan_backend = :streaming_fasta_simd, bucket_bases = 11,
                stream_chunk_bases = 64, prefilter_bits = 26,
                lookup_variant = lookup_variant, stats = variant_stats)
            lookup_outputs[lookup_variant] = read(output)
            lookup_stats[lookup_variant] = variant_stats
        end
        @test lookup_outputs[:bucketed] == lookup_outputs[:inline]
        @test lookup_outputs[:auto] == lookup_outputs[:inline]
        for field in (
                :motif_candidates, :prefix_hits, :guide_pairs,
                :alignment_calls, :distance_calls, :traceback_calls,
                :emitted_rows)
            @test getfield(lookup_stats[:bucketed], field) ==
                getfield(lookup_stats[:inline], field)
            @test getfield(lookup_stats[:auto], field) ==
                getfield(lookup_stats[:inline], field)
        end

        lookup_plain = joinpath(tdir, "scan_lookup_plain.csv")
        CHOPOFF.search_prefixHashScan(
            [guide3], fused_genome, motif3, lookup_plain;
            distance = 3, early_stopping = fill(100, 4),
            query_variant = :bitmask64,
            scan_backend = :streaming_fasta_simd, bucket_bases = 11,
            stream_chunk_bases = 64, prefilter_bits = 26,
            lookup_variant = :bucketed)
        @test read(lookup_plain) == lookup_outputs[:inline]

        query_build_outputs = Dict{Symbol, Vector{UInt8}}()
        query_build_stats = Dict{Symbol, CHOPOFF.PrefixHashScanStats}()
        query_build_guides = [
            guide3, LongDNA{4}("TGCATGCATGCATGCATGCA"), guide3]
        for query_build_backend in (:serial, :parallel, :auto)
            output = joinpath(
                tdir, "scan_query_" * string(query_build_backend) * ".csv")
            backend_stats = CHOPOFF.PrefixHashScanStats()
            CHOPOFF.search_prefixHashScan(
                query_build_guides, fused_genome, motif3, output;
                distance = 3, early_stopping = fill(100, 4),
                query_variant = :bitmask64, scan_backend = :auto,
                scan_threads = 4, stream_chunk_bases = 64,
                query_build_backend = query_build_backend, stats = backend_stats)
            query_build_outputs[query_build_backend] = read(output)
            query_build_stats[query_build_backend] = backend_stats
        end
        @test query_build_outputs[:parallel] == query_build_outputs[:serial]
        @test query_build_outputs[:auto] == query_build_outputs[:serial]
        for field in (
                :path_rows, :query_hashes, :motif_candidates, :prefix_hits,
                :guide_pairs, :alignment_calls, :distance_calls,
                :traceback_calls, :emitted_rows)
            @test getfield(query_build_stats[:parallel], field) ==
                getfield(query_build_stats[:serial], field)
            @test getfield(query_build_stats[:auto], field) ==
                getfield(query_build_stats[:serial], field)
        end

        plain_query_outputs = Dict{Symbol, Vector{UInt8}}()
        for query_build_backend in (:serial, :parallel)
            output = joinpath(
                tdir, "scan_query_plain_" * string(query_build_backend) * ".csv")
            CHOPOFF.search_prefixHashScan(
                query_build_guides, fused_genome, motif3, output;
                distance = 3, early_stopping = fill(100, 4),
                query_variant = :bitmask64, scan_backend = :auto,
                scan_threads = 4, stream_chunk_bases = 64,
                query_build_backend = query_build_backend)
            plain_query_outputs[query_build_backend] = read(output)
        end
        @test plain_query_outputs[:parallel] == plain_query_outputs[:serial]


        early_outputs = Dict{Symbol, Vector{UInt8}}()
        for backend in (
                :fused_fasta_simd, :streaming_fasta_simd,
                :streaming_fasta_simd_fused)
            output = joinpath(tdir, "scan_early_" * string(backend) * ".csv")
            CHOPOFF.search_prefixHashScan(
                [guide3],
                fused_genome,
                motif3,
                output;
                distance = 3,
                hash_len = hash_len,
                early_stopping = fill(1, 4),
                query_variant = :bitmask64,
                scan_backend = backend,
                bucket_bases = 8,
                stream_chunk_bases = 64,
            )
            early_outputs[backend] = read(output)
        end
        @test early_outputs[:streaming_fasta_simd] ==
            early_outputs[:fused_fasta_simd]
        @test early_outputs[:streaming_fasta_simd_fused] ==
            early_outputs[:streaming_fasta_simd]

        stream_stats = Dict{Symbol, CHOPOFF.PrefixHashScanStats}()
        stream_bytes = Dict{Symbol, Vector{UInt8}}()
        for backend in (:streaming_fasta_simd, :streaming_fasta_simd_fused)
            output = joinpath(tdir, "scan_stats_" * string(backend) * ".csv")
            backend_stats = CHOPOFF.PrefixHashScanStats()
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
                stream_chunk_bases = 64,
                stats = backend_stats,
            )
            stream_stats[backend] = backend_stats
            stream_bytes[backend] = read(output)
        end
        @test stream_bytes[:streaming_fasta_simd_fused] ==
            stream_bytes[:streaming_fasta_simd]
        for field in (
                :motif_candidates, :prefix_hits, :guide_pairs,
                :alignment_calls, :distance_calls, :traceback_calls,
                :emitted_rows)
            @test getfield(stream_stats[:streaming_fasta_simd_fused], field) ==
                getfield(stream_stats[:streaming_fasta_simd], field)
        end
        @test stream_stats[:streaming_fasta_simd].scan_backend ==
            :streaming_fasta_simd
        @test stream_stats[:streaming_fasta_simd_fused].scan_backend ==
            :streaming_fasta_simd_fused

        myers_output = joinpath(tdir, "scan_myers_raw.csv")
        myers_stats = CHOPOFF.PrefixHashScanStats()
        CHOPOFF.search_prefixHashScan(
            [guide3],
            fused_genome,
            motif3,
            myers_output;
            distance = 3,
            hash_len = hash_len,
            early_stopping = fill(100, 4),
            query_variant = :bitmask64,
            scan_backend = :fused_fasta_simd,
            bucket_bases = 8,
            verify_variant = :myers_raw,
            stats = myers_stats,
        )
        myers_df = DataFrame(CSV.File(myers_output))
        sort!(myers_df, names(myers_df))
        @test myers_df == outputs[:legacy]
        @test myers_stats.distance_calls > 0
        @test myers_stats.traceback_calls == myers_stats.emitted_rows
        @test_throws ErrorException CHOPOFF.search_prefixHashScan(
            [guide3], fused_genome, motif3, myers_output;
            distance = 3, early_stopping = fill(100, 4),
            lookup_variant = :invalid)
        @test_throws ErrorException CHOPOFF.search_prefixHashScan(
            [guide3], fused_genome, motif3, myers_output;
            distance = 3, early_stopping = fill(100, 4),
            scan_backend = :streaming_fasta_simd, bucket_bases = 11,
            prefilter_bits = 0, lookup_variant = :bucketed)
        @test_throws ErrorException CHOPOFF.search_prefixHashScan(
            [guide3], fused_genome, motif3, myers_output;
            distance = 3, early_stopping = fill(100, 4),
            scan_backend = :streaming_fasta_simd, bucket_bases = 8,
            lookup_variant = :bucketed)
        @test_throws ErrorException CHOPOFF.search_prefixHashScan(
            [guide3], fused_genome, motif3, myers_output;
            distance = 3, early_stopping = fill(100, 4),
            scan_backend = :streaming_fasta_simd_fused, bucket_bases = 11,
            lookup_variant = :bucketed)

        @test_throws ErrorException CHOPOFF.search_prefixHashScan(
            [guide3],
            fused_genome,
            motif3,
            myers_output;
            distance = 3,
            early_stopping = fill(100, 4),
            scan_backend = :fused_directory,
            verify_variant = :myers_raw,
        )
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
