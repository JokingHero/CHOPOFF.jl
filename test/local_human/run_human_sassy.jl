#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Dates
using FASTX
using Printf

const ROOT_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_GENOME = "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
const DEFAULT_GUIDES = joinpath(@__DIR__, "data", "guides_for_tests.txt")
const SAMPLE_GENOME = joinpath(ROOT_DIR, "test", "sample_data", "genome", "semirandom.fa")
const SAMPLE_GUIDES = joinpath(ROOT_DIR, "test", "sample_data", "guides.txt")
const CORE_PARITY_COLS = [:guide, :distance, :chromosome, :start, :strand]
const VALIDATED_PREFIX_REJECT_COLS = [
    :guide,
    :distance,
    :chromosome,
    :start,
    :strand,
    :reason,
]

function parse_bool_env(name::String, default::Bool)
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    raw in ("1", "true", "yes", "y") && return true
    raw in ("0", "false", "no", "n") && return false
    error("Invalid boolean in environment variable $name: '$raw'")
end

function parse_int_env(name::String, default::Int)
    raw = strip(get(ENV, name, string(default)))
    try
        return parse(Int, raw)
    catch
        error("Invalid integer in environment variable $name: '$raw'")
    end
end

function parse_early_stopping(distance::Int)
    raw = strip(get(ENV, "CHOPOFF_HUMAN_EARLY_STOPPING", ""))
    isempty(raw) && return fill(1_000_000, distance + 1)
    vals = parse.(Int, split(raw, ','))
    length(vals) == distance + 1 ||
        error("CHOPOFF_HUMAN_EARLY_STOPPING must have distance+1 comma-separated values")
    return vals
end

function require_file(path::String)
    isfile(path) || error("Required file not found: $path")
    return path
end

function validate_genome(path::String)
    require_file(path)
    if !endswith(lowercase(path), ".2bit")
        require_file(path * ".fai")
    end
    return path
end

function load_guides(path::String, motif::Motif)
    require_file(path)
    raw = String[]
    for line in eachline(path)
        s = uppercase(strip(line))
        isempty(s) && continue
        startswith(s, "#") && continue
        push!(raw, s)
    end
    isempty(raw) && error("No guides found in $path")

    expected_len = length_noPAM(motif)
    bad = findall(!=(expected_len), length.(raw))
    if !isempty(bad)
        first_bad = bad[1]
        error(
            "Guide length mismatch for motif $(motif.alias): expected $expected_len, " *
            "got $(length(raw[first_bad])) at guide $first_bad ($(raw[first_bad]))",
        )
    end
    return LongDNA{4}.(raw)
end

function count_result_rows(path::String)
    isfile(path) || return 0
    rows = 0
    for _ in eachline(path)
        rows += 1
    end
    return max(rows - 1, 0)
end

filesize_bytes(path::String) = isfile(path) ? stat(path).size : 0

function count_ambiguous_reference_rows(path::String)
    isfile(path) || return 0
    rows = 0
    for row in CSV.Rows(path)
        ref = getproperty(row, :alignment_reference)
        if any(c -> !(c in ('A', 'C', 'G', 'T', '-', 'a', 'c', 'g', 't')), ref)
            rows += 1
        end
    end
    return rows
end

function normalize_core_parity(path::String)
    df = DataFrame(CSV.File(path))
    cols = Set(Symbol.(names(df)))
    for c in CORE_PARITY_COLS
        c in cols || error("Missing required parity column '$c' in $path")
    end
    core = select(df, CORE_PARITY_COLS)
    core.guide = String.(core.guide)
    core.distance = Int.(core.distance)
    core.chromosome = String.(core.chromosome)
    core.start = Int.(core.start)
    core.strand = String.(core.strand)
    sort!(core, CORE_PARITY_COLS)
    return core
end

function write_parity_diffs(sassy_path::String, prefix_path::String, output_dir::String)
    return write_named_parity_diffs(sassy_path, prefix_path, output_dir, "sassy", "prefix")
end

function write_named_parity_diffs(lhs_path::String, rhs_path::String, output_dir::String, lhs_name::String, rhs_name::String)
    lhs_core = normalize_core_parity(lhs_path)
    rhs_core = normalize_core_parity(rhs_path)
    lhs_only = antijoin(lhs_core, rhs_core, on = CORE_PARITY_COLS)
    rhs_only = antijoin(rhs_core, lhs_core, on = CORE_PARITY_COLS)
    if nrow(lhs_only) == 0
        lhs_only = DataFrame(
            guide = String[],
            distance = Int[],
            chromosome = String[],
            start = Int[],
            strand = String[],
        )
    end
    if nrow(rhs_only) == 0
        rhs_only = DataFrame(
            guide = String[],
            distance = Int[],
            chromosome = String[],
            start = Int[],
            strand = String[],
        )
    end
    CSV.write(joinpath(output_dir, "parity_$(lhs_name)_only.csv"), lhs_only)
    CSV.write(joinpath(output_dir, "parity_$(rhs_name)_only.csv"), rhs_only)
    return nrow(lhs_only), nrow(rhs_only)
end


function valid_dna_or_gap(s::AbstractString)
    return all(c -> c in ('A', 'C', 'G', 'T', '-', 'a', 'c', 'g', 't'), s)
end

function pam_matches_string(ref::AbstractString, pam::AbstractString)
    length(ref) == length(pam) || return false
    for (r, p) in zip(ref, pam)
        iscompatible(DNA(r), DNA(p)) || return false
    end
    return true
end

function ref_ambig_within_limit(seq::LongDNA{4}, start_pos::Int, end_pos::Int, limit::Int)
    limit >= end_pos - start_pos + 1 && return true
    n_ambig = 0
    for i in start_pos:end_pos
        if isambiguous(seq[i])
            n_ambig += 1
            n_ambig > limit && return false
        end
    end
    return true
end

function revalidate_prefix_row(row, chrom_seq::LongDNA{4}, motif::Motif, distance::Int)
    guide = LongDNA{4}(String(row.guide))
    guide_len = length(guide)
    pam_len = length(motif.pam_loci_fwd)
    pos = Int(row.start)
    is_plus = String(row.strand) == "+"
    chrom_len = length(chrom_seq)

    if motif.extends5
        if is_plus
            motif_end = pos
            motif_start = motif_end - guide_len - pam_len + 1
            guide_start = motif_start
            guide_end = motif_end - pam_len
            pam_seq = String(motif.fwd[motif.pam_loci_fwd])
            pam_start = guide_end + 1
            pam_end = motif_end
        else
            motif_start = pos
            motif_end = motif_start + guide_len + pam_len - 1
            guide_start = motif_start + pam_len
            guide_end = motif_end
            pam_seq = String(motif.rve[motif.pam_loci_rve])
            pam_start = motif_start
            pam_end = motif_start + pam_len - 1
        end
    else
        return false, "unsupported_non_extends5", missing
    end

    if motif_start < 1 || motif_end > chrom_len || guide_start < 1 || guide_end > chrom_len
        return false, "out_of_bounds", missing
    end
    if !ref_ambig_within_limit(chrom_seq, motif_start, motif_end, motif.ambig_max)
        return false, "ambiguous_motif_window", missing
    end
    if !pam_matches_string(String(chrom_seq[pam_start:pam_end]), pam_seq)
        return false, "pam_mismatch", missing
    end

    guide_slice = chrom_seq[guide_start:guide_end]
    query_for_aln = reverse(guide)
    if is_plus
        ext = CHOPOFF.getExt5(chrom_seq, guide_start - 1, distance)
        ref_for_aln = reverse(ext * guide_slice)
    else
        ext = CHOPOFF.getExt3(chrom_seq, chrom_len, guide_end + 1, distance)
        ref_for_aln = complement(guide_slice * ext)
    end
    aln = align(query_for_aln, ref_for_aln, distance)
    if aln.dist > distance
        return false, "distance_recomputed_gt_k", missing
    end
    if !valid_dna_or_gap(aln.ref)
        return false, "ambiguous_consumed_reference", missing
    end
    return true, "", aln.dist
end

function empty_rejected_df()
    return DataFrame(
        guide = String[],
        distance = Int[],
        chromosome = String[],
        start = Int[],
        strand = String[],
        reason = String[],
    )
end

function write_validated_prefix_oracle(prefix_path::String, genome_path::String, motif::Motif, distance::Int, output_dir::String)
    prefix_df = DataFrame(CSV.File(prefix_path))
    if nrow(prefix_df) == 0
        validated_path = joinpath(output_dir, "prefixhash_validated.csv")
        CSV.write(validated_path, prefix_df)
        CSV.write(joinpath(output_dir, "prefixhash_rejected.csv"), empty_rejected_df())
        return validated_path, 0, 0
    end

    grouped_rows = Dict{String, Vector{Int}}()
    for (i, chrom) in enumerate(String.(prefix_df.chromosome))
        push!(get!(grouped_rows, chrom, Int[]), i)
    end

    accepted = falses(nrow(prefix_df))
    recomputed_dist = Vector{Union{Missing, Int}}(missing, nrow(prefix_df))
    reject_reason = fill("", nrow(prefix_df))

    open(genome_path, "r") do io
        reader = FASTA.Reader(io, index = genome_path * ".fai")
        for chrom in keys(grouped_rows)
            record = reader[chrom]
            chrom_seq = FASTA.sequence(LongDNA{4}, record)
            for row_idx in grouped_rows[chrom]
                ok, reason, dist = revalidate_prefix_row(prefix_df[row_idx, :], chrom_seq, motif, distance)
                accepted[row_idx] = ok
                reject_reason[row_idx] = reason
                recomputed_dist[row_idx] = dist
            end
        end
    end

    validated = prefix_df[accepted, :]
    validated.distance = Int.(recomputed_dist[accepted])
    rejected = prefix_df[.!accepted, CORE_PARITY_COLS]
    rejected.reason = reject_reason[.!accepted]

    validated_path = joinpath(output_dir, "prefixhash_validated.csv")
    rejected_path = joinpath(output_dir, "prefixhash_rejected.csv")
    CSV.write(validated_path, validated)
    CSV.write(rejected_path, rejected)
    return validated_path, nrow(validated), nrow(rejected)
end

function warmup_compile(distance::Int)
    isfile(SAMPLE_GENOME) || return
    isfile(SAMPLE_GENOME * ".fai") || return
    isfile(SAMPLE_GUIDES) || return

    motif = Motif("Cas9"; distance = min(distance, 3))
    guides = LongDNA{4}.(readlines(SAMPLE_GUIDES)[1:1])
    tmp = mktempdir(prefix = "chopoff_human_warmup_")
    try
        search_sassy(
            guides,
            SAMPLE_GENOME,
            motif,
            joinpath(tmp, "warmup.csv");
            distance = min(distance, 3),
            early_stopping = fill(10, min(distance, 3) + 1),
        )
    finally
        rm(tmp; recursive = true, force = true)
    end
end

function main()
    genome = validate_genome(String(strip(get(ENV, "CHOPOFF_HUMAN_GENOME", DEFAULT_GENOME))))
    guides_path = String(strip(get(ENV, "CHOPOFF_HUMAN_GUIDES", DEFAULT_GUIDES)))
    motif_name = String(strip(get(ENV, "CHOPOFF_HUMAN_MOTIF", "Cas9")))
    distance = parse_int_env("CHOPOFF_HUMAN_DISTANCE", 3)
    compare_prefix = parse_bool_env("CHOPOFF_HUMAN_COMPARE_PREFIX", true)
    compare_scan = parse_bool_env("CHOPOFF_HUMAN_COMPARE_SCAN", true)
    scan_query_variant = Symbol(strip(get(ENV, "CHOPOFF_HUMAN_SCAN_QUERY_VARIANT", "auto")))
    scan_backend = Symbol(strip(get(ENV, "CHOPOFF_HUMAN_SCAN_BACKEND", "auto")))
    scan_bucket_bases = parse_int_env("CHOPOFF_HUMAN_SCAN_BUCKET_BASES", 11)
    scan_threads = parse_int_env("CHOPOFF_HUMAN_SCAN_THREADS", Threads.nthreads())
    scan_prefilter_bits = parse_int_env("CHOPOFF_HUMAN_SCAN_PREFILTER_BITS", 26)
    scan_verify_variant = Symbol(strip(get(ENV, "CHOPOFF_HUMAN_SCAN_VERIFY_VARIANT", "auto")))
    rebuild_prefix = parse_bool_env("CHOPOFF_HUMAN_REBUILD_PREFIX", false)
    keep_outputs = parse_bool_env("CHOPOFF_HUMAN_KEEP_OUTPUTS", true)
    use_avx512 = parse_bool_env("CHOPOFF_HUMAN_USE_AVX512", false)
    force_safe_minima = parse_bool_env("CHOPOFF_HUMAN_FORCE_SAFE_MINIMA", false)

    motif = Motif(motif_name; distance = distance)
    guides = load_guides(guides_path, motif)
    early_stopping = parse_early_stopping(distance)

    stamp = Dates.format(now(Dates.UTC), dateformat"yyyymmdd_HHMMSS")
    output_dir = joinpath(@__DIR__, "outputs", stamp)
    genome_key = replace(basename(genome), r"[^A-Za-z0-9_.-]" => "_")
    index_dir = joinpath(@__DIR__, "indexes", "$(genome_key)_$(motif_name)_d$(distance)")
    mkpath(output_dir)

    println("Human SASSY benchmark")
    println("time_utc: ", Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println("genome: ", genome)
    println("guides: ", guides_path)
    println("guide_count: ", length(guides))
    println("motif: ", motif_name)
    println("distance: ", distance)
    println("threads: ", Threads.nthreads())
    println("compare_prefix: ", compare_prefix)
    println("compare_scan: ", compare_scan)
    println("scan_query_variant: ", scan_query_variant)
    println("scan_backend: ", scan_backend)
    println("scan_bucket_bases: ", scan_bucket_bases)
    println("scan_threads: ", scan_threads)
    println("scan_prefilter_bits: ", scan_prefilter_bits)
    println("scan_verify_variant: ", scan_verify_variant)
    println("output_dir: ", output_dir)

    warmup_compile(distance)

    sassy_path = joinpath(output_dir, "sassy.csv")
    sassy_elapsed = @elapsed search_sassy(
        guides,
        genome,
        motif,
        sassy_path;
        distance = distance,
        early_stopping = early_stopping,
        use_avx512 = use_avx512,
        force_safe_minima = force_safe_minima,
    )
    sassy_rows = count_result_rows(sassy_path)
    sassy_bytes = filesize_bytes(sassy_path)
    sassy_ambig_ref_rows = count_ambiguous_reference_rows(sassy_path)

    prefix_elapsed = missing
    prefix_rows = missing
    prefix_bytes = missing
    prefix_ambig_ref_rows = missing
    parity_sassy_only = missing
    parity_prefix_only = missing
    validated_prefix_path = ""
    validated_prefix_rows = missing
    rejected_prefix_rows = missing
    scan_elapsed = missing
    scan_rows = missing
    scan_bytes = missing
    scan_ambig_ref_rows = missing
    scan_path_source = missing
    scan_query_variant_used = missing
    scan_backend_used = missing
    scan_alignment_calls = missing
    scan_distance_calls = missing
    scan_traceback_calls = missing
    scan_metadata_s = missing
    scan_query_build_s = missing
    scan_path_load_s = missing
    scan_record_io_s = missing
    scan_sequence_convert_s = missing
    scan_scan_s = missing
    scan_findguides_s = missing
    scan_candidate_hash_s = missing
    scan_align_s = missing
    parity_scan_only = missing
    parity_prefix_only_vs_scan = missing

    if compare_prefix
        mkpath(index_dir)
        db_file = joinpath(index_dir, "prefixHashDB.bin")
        if rebuild_prefix || !isfile(db_file)
            build_prefixHashDB("human_$(motif_name)_d$(distance)", genome, motif, index_dir)
        else
            println("Reusing prefixHashDB index: ", index_dir)
        end

        prefix_path = joinpath(output_dir, "prefixhash.csv")
        prefix_elapsed = @elapsed search_prefixHashDB(
            index_dir,
            guides,
            prefix_path;
            distance = distance,
            early_stopping = early_stopping,
        )
        prefix_rows = count_result_rows(prefix_path)
        prefix_bytes = filesize_bytes(prefix_path)
        prefix_ambig_ref_rows = count_ambiguous_reference_rows(prefix_path)
        validated_prefix_path, validated_prefix_rows, rejected_prefix_rows =
            write_validated_prefix_oracle(prefix_path, genome, motif, distance, output_dir)
        parity_sassy_only, parity_prefix_only =
            write_parity_diffs(sassy_path, validated_prefix_path, output_dir)
    end

    if compare_scan
        scan_path = joinpath(output_dir, "prefixhashscan.csv")
        scan_stats = CHOPOFF.PrefixHashScanStats()
        scan_elapsed = @elapsed CHOPOFF.search_prefixHashScan(
            guides,
            genome,
            motif,
            scan_path;
            distance = distance,
            early_stopping = early_stopping,
            query_variant = scan_query_variant,
            scan_backend = scan_backend,
            bucket_bases = scan_bucket_bases,
            scan_threads = scan_threads,
            prefilter_bits = scan_prefilter_bits,
            verify_variant = scan_verify_variant,
            stats = scan_stats,
        )
        scan_rows = count_result_rows(scan_path)
        scan_bytes = filesize_bytes(scan_path)
        scan_ambig_ref_rows = count_ambiguous_reference_rows(scan_path)
        scan_path_source = scan_stats.path_source
        scan_query_variant_used = scan_stats.query_variant
        scan_backend_used = scan_stats.scan_backend
        scan_alignment_calls = scan_stats.alignment_calls
        scan_distance_calls = scan_stats.distance_calls
        scan_traceback_calls = scan_stats.traceback_calls
        scan_metadata_s = scan_stats.metadata_ns / 1e9
        scan_query_build_s = scan_stats.query_build_ns / 1e9
        scan_path_load_s = scan_stats.path_load_ns / 1e9
        scan_record_io_s = scan_stats.record_io_ns / 1e9
        scan_sequence_convert_s = scan_stats.sequence_convert_ns / 1e9
        scan_scan_s = scan_stats.scan_ns / 1e9
        scan_findguides_s = scan_stats.findguides_ns / 1e9
        scan_candidate_hash_s = scan_stats.candidate_hash_ns / 1e9
        scan_align_s = scan_stats.align_ns / 1e9
        scan_summary = DataFrame([scan_stats])
        CSV.write(joinpath(output_dir, "prefixhashscan_stats.csv"), scan_summary)
        if !isempty(validated_prefix_path)
            parity_scan_only, parity_prefix_only_vs_scan =
                write_named_parity_diffs(scan_path, validated_prefix_path, output_dir, "scan", "prefix_vs_scan")
        end
    end

    summary = DataFrame([(
        timestamp = Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
        genome = genome,
        guides = guides_path,
        guide_count = length(guides),
        motif = motif_name,
        distance = distance,
        threads = Threads.nthreads(),
        use_avx512 = use_avx512,
        force_safe_minima = force_safe_minima,
        early_stopping = join(early_stopping, ","),
        sassy_elapsed_s = sassy_elapsed,
        sassy_rows = sassy_rows,
        sassy_bytes = sassy_bytes,
        sassy_ambig_ref_rows = sassy_ambig_ref_rows,
        compare_prefix = compare_prefix,
        compare_scan = compare_scan,
        scan_query_variant = scan_query_variant,
        scan_backend = scan_backend,
        scan_bucket_bases = scan_bucket_bases,
        scan_threads = scan_threads,
        scan_prefilter_bits = scan_prefilter_bits,
        scan_verify_variant = scan_verify_variant,
        prefix_elapsed_s = prefix_elapsed,
        prefix_rows = prefix_rows,
        prefix_bytes = prefix_bytes,
        prefix_ambig_ref_rows = prefix_ambig_ref_rows,
        validated_prefix_rows = validated_prefix_rows,
        rejected_prefix_rows = rejected_prefix_rows,
        parity_sassy_only = parity_sassy_only,
        parity_prefix_only = parity_prefix_only,
        scan_elapsed_s = scan_elapsed,
        scan_rows = scan_rows,
        scan_bytes = scan_bytes,
        scan_ambig_ref_rows = scan_ambig_ref_rows,
        scan_path_source = scan_path_source,
        scan_query_variant_used = scan_query_variant_used,
        scan_backend_used = scan_backend_used,
        scan_alignment_calls = scan_alignment_calls,
        scan_distance_calls = scan_distance_calls,
        scan_traceback_calls = scan_traceback_calls,
        scan_metadata_s = scan_metadata_s,
        scan_query_build_s = scan_query_build_s,
        scan_path_load_s = scan_path_load_s,
        scan_record_io_s = scan_record_io_s,
        scan_sequence_convert_s = scan_sequence_convert_s,
        scan_scan_s = scan_scan_s,
        scan_findguides_s = scan_findguides_s,
        scan_candidate_hash_s = scan_candidate_hash_s,
        scan_align_s = scan_align_s,
        parity_scan_only = parity_scan_only,
        parity_prefix_only_vs_scan = parity_prefix_only_vs_scan,
        output_dir = output_dir,
        index_dir = compare_prefix ? index_dir : "",
    )])

    summary_path = joinpath(output_dir, "summary.csv")
    CSV.write(summary_path, summary)

    @printf("SASSY elapsed: %.3fs | rows=%d | bytes=%d
", sassy_elapsed, sassy_rows, sassy_bytes)
    @printf("SASSY ambiguous-reference rows: %d
", sassy_ambig_ref_rows)
    if compare_scan
        @printf(
            "prefixHashScan elapsed: %.3fs | rows=%d | bytes=%d | ambiguous_reference_rows=%d | path_source=%s | query_variant=%s | backend=%s | verify=%s | verification_calls=%d | distance_calls=%d | tracebacks=%d
",
            scan_elapsed,
            scan_rows,
            scan_bytes,
            scan_ambig_ref_rows,
            string(scan_path_source),
            string(scan_query_variant_used),
            string(scan_backend_used),
            string(scan_verify_variant),
            scan_alignment_calls,
            scan_distance_calls,
            scan_traceback_calls,
        )
        @printf(
            "prefixHashScan timers: metadata=%.3fs | query_build=%.3fs | path_load=%.3fs | record_io=%.3fs | sequence_convert=%.3fs | scan=%.3fs | findguides=%.3fs | candidate_hash=%.3fs | align=%.3fs
",
            scan_metadata_s,
            scan_query_build_s,
            scan_path_load_s,
            scan_record_io_s,
            scan_sequence_convert_s,
            scan_scan_s,
            scan_findguides_s,
            scan_candidate_hash_s,
            scan_align_s,
        )
        if !ismissing(parity_scan_only)
            @printf(
                "Parity prefixHashScan vs validated prefixHashDB: scan_only=%s | prefix_only=%s
",
                string(parity_scan_only),
                string(parity_prefix_only_vs_scan),
            )
        end
    end
    if compare_prefix
        @printf(
            "prefixHashDB elapsed: %.3fs | raw_rows=%d | validated_rows=%s | rejected_rows=%s | bytes=%d | ambiguous_reference_rows=%d
",
            prefix_elapsed,
            prefix_rows,
            string(validated_prefix_rows),
            string(rejected_prefix_rows),
            prefix_bytes,
            prefix_ambig_ref_rows,
        )
        @printf(
            "Parity vs validated prefixHashDB: sassy_only=%s | prefix_only=%s
",
            string(parity_sassy_only),
            string(parity_prefix_only),
        )
    end
    println("summary: ", summary_path)

    if !keep_outputs
        rm(output_dir; recursive = true, force = true)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
