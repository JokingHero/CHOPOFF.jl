#!/usr/bin/env julia

using CHOPOFF
using BioSequences
using CSV
using DataFrames
using Dates
using FASTX
using SHA
using TwoBit

include(joinpath(@__DIR__, "..", "helpers", "prefix_parity.jl"))
using .PrefixParity

const ROOT_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_GENOME = "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
const STANDARD_CHROMS = vcat(string.(1:22), ["X", "Y"])
const DETAIL_HEADER =
    "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n"

struct GenericParityCase
    id::String
    motif::Motif
    reference_paths::Vector{String}
    guide_count::Int
    distances::Vector{Int}
    ambiguities::Vector{Int}
    build_index::Bool
    baseline_id::String
    baseline_strand::Union{Nothing, String}
    chunk_bases::Int
end

function parse_bool_env(name::String, default::Bool)
    value = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    value in ("1", "true", "yes") && return true
    value in ("0", "false", "no") && return false
    error("Invalid boolean $name=$value")
end

function output_dir()
    configured = strip(get(ENV, "CHOPOFF_GENERIC_PARITY_OUT", ""))
    if !isempty(configured)
        return abspath(configured)
    end
    stamp = Dates.format(now(Dates.UTC), dateformat"yyyymmdd_HHMMSS")
    return joinpath(@__DIR__, "outputs", "generic_parity_" * stamp)
end

mode() = Symbol(lowercase(strip(get(ENV, "CHOPOFF_GENERIC_PARITY_MODE", "smoke"))))
genome_path() = abspath(strip(get(ENV, "CHOPOFF_HUMAN_GENOME", DEFAULT_GENOME)))
case_dir(out::String, case::GenericParityCase) = joinpath(out, "cases", case.id)
guide_path(out::String, case::GenericParityCase) = joinpath(case_dir(out, case), "guides.txt")
source_path(out::String, case::GenericParityCase) = joinpath(case_dir(out, case), "guide_sources.csv")
index_dir(out::String, case::GenericParityCase) = joinpath(out, "indexes", case.id)

function motif_25n_ngg()
    return Motif(
        "25N_NGG", repeat("N", 25) * "XXX",
        repeat("X", 25) * "NGG", true, true, 4, true, 3)
end

function motif_internal()
    return Motif(
        "10N_AGG_10N", repeat("N", 10) * "XXX" * repeat("N", 10),
        repeat("X", 10) * "AGG" * repeat("X", 10),
        true, true, 4, true, 3)
end

function motif_pamless()
    return Motif(
        "20N_pamless", repeat("N", 20), repeat("X", 20),
        true, true, 4, true, 3)
end

function motif_right_ngg(length_::Int)
    return Motif(
        "$(length_)N_NGG", repeat("N", length_) * "XXX",
        repeat("X", length_) * "NGG", true, true, 4, true, 3)
end

function strand_motif(forward::Bool)
    return Motif(
        forward ? "Cas9_NGA_forward" : "Cas9_NGA_reverse",
        repeat("N", 20) * "XXX", repeat("X", 20) * "NGA",
        forward, !forward, 4, true, 3)
end

function write_fasta(path::String, name::String, seq::LongDNA{4})
    mkpath(dirname(path))
    width = 60
    open(path, "w") do io
        write(io, ">", name, "\n")
        text = string(seq)
        for first_idx in 1:width:length(text)
            write(io, text[first_idx:min(first_idx + width - 1, end)], "\n")
        end
    end
    open(path * ".fai", "w") do io
        write(io, name, "\t", string(length(seq)), "\t",
            string(ncodeunits(name) + 2), "\t", string(width), "\t",
            string(width + 1), "\n")
    end
    return path
end

function write_twobit(path::String, name::String, seq::LongDNA{4})
    writer = TwoBit.Writer(open(path, "w"), [name])
    write(writer, TwoBit.Record(name, seq))
    close(writer)
    return path
end

function read_chrom(path::String, chrom::String)
    ref = open(path, "r")
    try
        if endswith(lowercase(path), ".2bit")
            return TwoBit.sequence(LongDNA{4}, TwoBit.Reader(ref)[chrom])
        end
        reader = FASTA.Reader(ref, index = path * ".fai")
        return FASTA.sequence(LongDNA{4}, reader[chrom])
    finally
        close(ref)
    end
end

function write_bounded_reference(out::String, genome::String)
    ref_dir = joinpath(out, "references")
    fasta = joinpath(ref_dir, "grch38_chr21_4m.fa")
    twobit = joinpath(ref_dir, "grch38_chr21_4m.2bit")
    if !isfile(fasta) || !isfile(twobit)
        chrom = read_chrom(genome, "21")
        first_pos = 10_000_001
        seq = chrom[first_pos:(first_pos + 4 * 1024 * 1024 - 1)]
        write_fasta(fasta, "21_slice", seq)
        write_twobit(twobit, "21_slice", seq)
    end
    return fasta, twobit
end

function case_fingerprint(case::GenericParityCase, reference_crc32c::String)
    payload = join((
        case.id,
        string(case.motif.fwd),
        string(case.motif.rve),
        string(case.motif.pam_loci_fwd),
        string(case.motif.pam_loci_rve),
        string(case.motif.extends5),
        string(maximum(case.distances)),
        "3",
        string(min(length_noPAM(case.motif) - maximum(case.distances), 16)),
        reference_crc32c,
    ), "|")
    return bytes2hex(sha256(payload))
end

function write_iupac_reference(out::String)
    fasta = joinpath(out, "references", "iupac_generic.fa")
    isfile(fasta) && return fasta
    guide = collect("ACGTACGTACGTACGTACGT")
    codes = [('R', 1), ('Y', 2), ('S', 3), ('W', 4), ('K', 3), ('M', 1),
        ('B', 2), ('D', 1), ('H', 1), ('V', 1), ('N', 1)]
    sites = String[]
    cursor = 1
    for count in (0, 1, 2, 3, 4, 1)
        candidate = copy(guide)
        for offset in 1:count
            code, position = codes[cursor]
            cursor += 1
            candidate[(offset - 1) * 4 + position] = code
        end
        pam = count == 3 ? "RGA" : "AGA"
        push!(sites, String(candidate) * pam)
    end
    # Same-length rotation has edit distance two through one deletion/insertion.
    push!(sites, String(guide[2:end]) * string(first(guide)) * "AGA")
    # Force one candidate to cross a 64-base chunk boundary.
    sequence = sites[1] * repeat("C", 59 - length(sites[1])) *
        join(sites[2:end], repeat("C", 41))
    write_fasta(fasta, "iupac", LongDNA{4}(sequence))
    return fasta
end

function write_smoke_reference(out::String)
    fasta = joinpath(out, "references", "smoke_generic.fa")
    isfile(fasta) && return fasta
    guide = "ACGTACGTACGTACGTACGT"
    ambiguous = guide[1:8] * "N" * guide[10:end]
    seq = LongDNA{4}(repeat("C", 40) * guide * "AGA" * repeat("C", 40) *
        ambiguous * "AGA" * repeat("C", 40))
    write_fasta(fasta, "smoke", seq)
    return fasta
end

function qualification_cases(out::String, genome::String)
    bounded_fasta, bounded_twobit = write_bounded_reference(out, genome)
    iupac = write_iupac_reference(out)
    full = [
        GenericParityCase("Cas9_NGA", Motif("Cas9_NGA"; distance = 4, ambig_max = 3),
            [genome], 65, collect(0:4), collect(0:3), true, "Cas9_NGA", nothing, 2 * 1024 * 1024),
        GenericParityCase("CasX", Motif("CasX"; distance = 4, ambig_max = 3),
            [genome], 65, collect(0:4), collect(0:3), true, "CasX", nothing, 2 * 1024 * 1024),
        GenericParityCase("25N_NGG", motif_25n_ngg(),
            [genome], 65, collect(0:4), collect(0:3), true, "25N_NGG", nothing, 2 * 1024 * 1024),
    ]
    subsets = [
        GenericParityCase("Cas9_NGA_forward", strand_motif(true), [genome], 65,
            collect(0:4), [0, 3], false, "Cas9_NGA", "+", 2 * 1024 * 1024),
        GenericParityCase("Cas9_NGA_reverse", strand_motif(false), [genome], 65,
            collect(0:4), [0, 3], false, "Cas9_NGA", "-", 2 * 1024 * 1024),
    ]
    bounded = [
        GenericParityCase("internal_pam", motif_internal(), [bounded_fasta, bounded_twobit],
            8, collect(0:4), collect(0:3), true, "internal_pam", nothing, 2 * 1024 * 1024),
        GenericParityCase("pamless", motif_pamless(), [bounded_fasta, bounded_twobit],
            8, collect(0:4), collect(0:3), true, "pamless", nothing, 2 * 1024 * 1024),
        GenericParityCase("guide16", motif_right_ngg(16), [bounded_fasta, bounded_twobit],
            8, collect(0:4), collect(0:3), true, "guide16", nothing, 2 * 1024 * 1024),
        GenericParityCase("iupac", Motif("Cas9_NGA"; distance = 4, ambig_max = 3),
            [iupac], 1, collect(0:4), collect(0:3), true, "iupac", nothing, 64),
    ]
    return vcat(full, subsets, bounded)
end

function smoke_cases(out::String)
    fasta = write_smoke_reference(out)
    return [GenericParityCase(
        "smoke", Motif("Cas9_NGA"; distance = 1, ambig_max = 1),
        [fasta], 1, [0, 1], [0, 1], true, "smoke", nothing, 64)]
end

function cases(out::String)
    selected = mode()
    selected == :smoke && return smoke_cases(out)
    selected == :qualification || error(
        "CHOPOFF_GENERIC_PARITY_MODE must be smoke or qualification")
    genome = genome_path()
    isfile(genome) || error("Genome not found: $genome")
    isfile(genome * ".fai") || error("FASTA index not found: $(genome).fai")
    return qualification_cases(out, genome)
end

function exact_pattern_match(pattern::LongDNA{4}, seq::LongDNA{4}, first_pos::Int)
    last_pos = first_pos + length(pattern) - 1
    1 <= first_pos && last_pos <= length(seq) || return false
    @inbounds for offset in 0:(length(pattern) - 1)
        base = seq[first_pos + offset]
        isambiguous(base) && return false
        iscompatible(pattern[offset + 1], base) || return false
    end
    return true
end

function guide_at(seq::LongDNA{4}, motif::Motif, first_pos::Int, isplus::Bool)
    pattern = isplus ? motif.fwd : motif.rve
    pam = isplus ? motif.pam_loci_fwd : motif.pam_loci_rve
    candidate = seq[first_pos:(first_pos + length(pattern) - 1)]
    guide = CHOPOFF.removepam(candidate, pam)
    return isplus ? guide : reverse_complement(guide)
end

function nearest_guide(
    seq::LongDNA{4}, motif::Motif, target::Int, isplus::Bool,
    used::Set{String})

    pattern = isplus ? motif.fwd : motif.rve
    isempty(pattern) && return nothing
    max_start = length(seq) - length(pattern) + 1
    for delta in 0:max_start
        positions = delta == 0 ? (target,) : (target + delta, target - delta)
        for first_pos in positions
            1 <= first_pos <= max_start || continue
            exact_pattern_match(pattern, seq, first_pos) || continue
            guide = guide_at(seq, motif, first_pos, isplus)
            guide_string = string(guide)
            guide_string in used && continue
            return (guide = guide_string, first_pos = first_pos,
                strand = isplus ? "+" : "-")
        end
    end
    return nothing
end

function sample_guides(
    reference::String, motif::Motif, count::Int, chromosomes::Vector{String})

    quotas = fill(div(count, length(chromosomes)), length(chromosomes))
    for idx in 1:rem(count, length(chromosomes))
        quotas[idx] += 1
    end
    used = Set{String}()
    rows = NamedTuple[]
    ordinal = 0
    for (chrom_idx, chrom) in enumerate(chromosomes)
        quota = quotas[chrom_idx]
        quota == 0 && continue
        seq = read_chrom(reference, chrom)
        for slot in 1:quota
            ordinal += 1
            target = clamp(round(Int, slot * length(seq) / (quota + 1)), 1, length(seq))
            desired_plus = isodd(ordinal)
            found = nearest_guide(seq, motif, target, desired_plus, used)
            found === nothing &&
                (found = nearest_guide(seq, motif, target, !desired_plus, used))
            found === nothing && error(
                "Could not sample guide for $(motif.alias) on $chrom")
            push!(used, found.guide)
            push!(rows, (
                guide = found.guide, chromosome = chrom,
                candidate_start = found.first_pos, strand = found.strand))
        end
    end
    length(rows) == count || error("Sampled $(length(rows))/$count guides")
    return DataFrame(rows)
end

function write_guides(out::String, case::GenericParityCase, rows::DataFrame)
    dir = case_dir(out, case)
    mkpath(dir)
    open(guide_path(out, case), "w") do io
        for guide in rows.guide
            println(io, guide)
        end
    end
    CSV.write(source_path(out, case), rows)
    return nothing
end

function prepare_guides!(out::String, all_cases::Vector{GenericParityCase})
    by_id = Dict(case.id => case for case in all_cases)
    for case in all_cases
        isfile(guide_path(out, case)) && continue
        if !case.build_index
            base = by_id[case.baseline_id]
            write_guides(out, case, DataFrame(CSV.File(source_path(out, base))))
        elseif case.id == "smoke" || case.id == "iupac"
            guide = "ACGTACGTACGTACGTACGT"
            write_guides(out, case, DataFrame(
                guide = [guide], chromosome = [case.id],
                candidate_start = [case.id == "smoke" ? 41 : 1], strand = ["+"]))
        else
            reference = first(case.reference_paths)
            chroms = occursin("chr21_4m", reference) ? ["21_slice"] : STANDARD_CHROMS
            write_guides(out, case,
                sample_guides(reference, case.motif, case.guide_count, chroms))
        end
    end
end

function run_prepare(out::String)
    mkpath(out)
    all_cases = cases(out)
    prepare_guides!(out, all_cases)
    checksums = Dict{String, String}()
    for reference in unique(first(case.reference_paths) for case in all_cases)
        checksums[reference] = string(open(CHOPOFF.crc32c, reference))
    end
    rows = [(
        case_id = case.id,
        motif = case.motif.alias,
        fwd = string(case.motif.fwd),
        rve = string(case.motif.rve),
        extends5 = case.motif.extends5,
        guide_bases = length_noPAM(case.motif),
        guide_count = case.guide_count,
        references = join(case.reference_paths, ";"),
        distances = join(case.distances, ","),
        ambiguities = join(case.ambiguities, ","),
        build_index = case.build_index,
        baseline_id = case.baseline_id,
        baseline_strand = something(case.baseline_strand, "both"),
        reference_crc32c = checksums[first(case.reference_paths)],
        fingerprint = case_fingerprint(
            case, checksums[first(case.reference_paths)]),
    ) for case in all_cases]
    CSV.write(joinpath(out, "manifest.csv"), DataFrame(rows))
    println("Prepared generic parity matrix: $out")
end

function load_guides(out::String, case::GenericParityCase)
    return LongDNA{4}.([strip(line) for line in eachline(guide_path(out, case))
        if !isempty(strip(line))])
end

function directory_size(path::String)
    return sum(stat(joinpath(root, file)).size
        for (root, _, files) in walkdir(path) for file in files; init = 0)
end

function run_build(out::String)
    mode() == :qualification && Threads.nthreads() != 1 &&
        error("Qualification prefixHashDB builds require one Julia thread")
    rebuild = parse_bool_env("CHOPOFF_GENERIC_PARITY_REBUILD", false)
    manifest = DataFrame(CSV.File(joinpath(out, "manifest.csv")))
    rows = NamedTuple[]
    for case in cases(out)
        case.build_index || continue
        dbdir = index_dir(out, case)
        dbfile = joinpath(dbdir, "prefixHashDB.bin")
        manifest_idx = findfirst(==(case.id), String.(manifest.case_id))
        manifest_idx === nothing && error("Missing manifest row for $(case.id)")
        fingerprint = String(manifest[manifest_idx, :fingerprint])
        db_manifest = joinpath(dbdir, "qualification_manifest.csv")
        reusable = false
        if isfile(dbfile) && isfile(db_manifest)
            stored = DataFrame(CSV.File(db_manifest))
            reusable = nrow(stored) == 1 &&
                String(stored[1, :fingerprint]) == fingerprint
        end
        if isfile(dbfile) && !reusable && !rebuild
            error("Index metadata mismatch for $(case.id); set CHOPOFF_GENERIC_PARITY_REBUILD=1")
        end
        reused = reusable && !rebuild
        elapsed = 0.0
        if !reused
            if isdir(dbdir)
                stale = dbdir * ".stale_" *
                    Dates.format(now(Dates.UTC), dateformat"yyyymmdd_HHMMSS")
                mv(dbdir, stale)
            end
            mkpath(dbdir)
            motif = setambig(setdist(case.motif, maximum(case.distances)), 3)
            elapsed = @elapsed build_prefixHashDB(
                "generic_parity_$(case.id)", first(case.reference_paths),
                motif, dbdir)
            CSV.write(db_manifest, DataFrame(
                case_id = [case.id], fingerprint = [fingerprint]))
        end
        push!(rows, (case_id = case.id, reused = reused,
            elapsed_s = elapsed, index_bytes = directory_size(dbdir), index_dir = dbdir))
        CSV.write(joinpath(out, "builds.csv"), DataFrame(rows))
        println("build case=$(case.id) reused=$reused elapsed=$(round(elapsed; digits=3))s")
        flush(stdout)
    end
end

function atomic_search(search::Function, path::String)
    mkpath(dirname(path))
    temporary = path * ".tmp.$(getpid())"
    search(temporary)
    mv(temporary, path; force = true)
    return nothing
end

reference_label(path::String) = replace(basename(path), r"[^A-Za-z0-9]+" => "_")

function run_scan(
    guides, reference, motif, output_path, distance, ambiguity, output_mode,
    chunk_bases)

    configured = setambig(setdist(motif, distance), ambiguity)
    limits = fill(typemax(Int), distance + 1)
    if chunk_bases == 2 * 1024 * 1024 || length(guides) > 64
        return search_prefixHashScan(
            guides, reference, output_path; motif = configured,
            distance = distance, early_stopping = limits,
            output = output_mode, verbose = output_mode == :detail)
    end
    hash_len = min(length_noPAM(configured) - distance, 16)
    return CHOPOFF.search_prefixHashScan(
        guides, reference, configured, output_path;
        distance = distance, hash_len = hash_len,
        early_stopping = limits, output = output_mode,
        stream_chunk_bases = chunk_bases, verbose = output_mode == :detail)
end

function expected_scan_backend(motif::Motif, reference::String, distance::Int)
    hash_len = min(length_noPAM(motif) - distance, 16)
    geometry = CHOPOFF.resolve_prefix_scan_geometry(motif, distance, hash_len)
    geometry === nothing && return :legacy
    CHOPOFF.prefix_scan_kind(geometry) == :generic ||
        error("Expected generic geometry for $(motif.alias) at distance $distance")
    CHOPOFF.can_use_prefix_hash_scan_simd() || return :fused_directory
    return endswith(lowercase(reference), ".2bit") ?
        :streaming_2bit_simd : :streaming_fasta_simd
end

function run_search(out::String)
    isfile(joinpath(out, "manifest.csv")) || error("Run prepare stage first")
    timing_rows = NamedTuple[]
    all_cases = cases(out)
    by_id = Dict(case.id => case for case in all_cases)
    for case in all_cases
        guides = load_guides(out, case)
        if case.build_index
            for distance in case.distances
                path = joinpath(case_dir(out, case), "prefix_d$(distance).csv")
                elapsed = @elapsed atomic_search(path) do temporary
                    search_prefixHashDB(
                        index_dir(out, case), guides, temporary;
                        distance = distance,
                        early_stopping = fill(typemax(Int), distance + 1))
                end
                push!(timing_rows, (case_id = case.id, reference = first(case.reference_paths),
                    distance = distance, ambiguity = 3, output = "detail",
                    algorithm = "prefixHashDB", expected_backend = "n/a",
                    elapsed_s = elapsed, path = path))
            end
        else
            base = by_id[case.baseline_id]
            isfile(joinpath(case_dir(out, base), "prefix_d$(first(case.distances)).csv")) ||
                error("Baseline search must precede $(case.id)")
        end
        for reference in case.reference_paths, distance in case.distances,
                ambiguity in case.ambiguities, output_mode in (:detail, :counts)
            label = reference_label(reference)
            path = joinpath(case_dir(out, case),
                "scan_$(label)_d$(distance)_a$(ambiguity)_$(output_mode).csv")
            backend = expected_scan_backend(
                setambig(setdist(case.motif, distance), ambiguity),
                reference, distance)
            elapsed = @elapsed atomic_search(path) do temporary
                run_scan(guides, reference, case.motif, temporary,
                    distance, ambiguity, output_mode, case.chunk_bases)
            end
            push!(timing_rows, (case_id = case.id, reference = reference,
                distance = distance, ambiguity = ambiguity,
                output = string(output_mode), algorithm = "prefixHashScan",
                expected_backend = string(backend),
                elapsed_s = elapsed, path = path))
            CSV.write(joinpath(out, "timings.csv"), DataFrame(timing_rows))
            println("search case=$(case.id) ref=$label d=$distance a=$ambiguity mode=$output_mode elapsed=$(round(elapsed; digits=3))s")
            flush(stdout)
        end
    end
end

mutable struct ReferenceCache
    path::String
    sequences::Dict{String, LongDNA{4}}
end

ReferenceCache(path::String) = ReferenceCache(path, Dict{String, LongDNA{4}}())
cached_chrom(cache::ReferenceCache, chrom::String) =
    get!(cache.sequences, chrom) do
        read_chrom(cache.path, chrom)
    end

function candidate_range(row, motif::Motif)
    span = length(motif)
    isplus = String(row.strand) == "+"
    anchor_end = isplus == motif.extends5
    anchor = Int(row.start)
    return anchor_end ? ((anchor - span + 1):anchor) :
        (anchor:(anchor + span - 1))
end

function normalized_candidate(
    chrom::LongDNA{4}, range::UnitRange{Int}, motif::Motif,
    is_antisense::Bool, distance::Int)

    pam = is_antisense ? motif.pam_loci_rve : motif.pam_loci_fwd
    guide = CHOPOFF.removepam(chrom[range], pam)
    if distance > 0
        if motif.extends5 && is_antisense
            guide *= CHOPOFF.getExt3(chrom, length(chrom), last(range) + 1, distance)
        elseif motif.extends5
            guide = CHOPOFF.getExt5(chrom, first(range) - 1, distance) * guide
        elseif is_antisense
            guide = CHOPOFF.getExt5(chrom, first(range) - 1, distance) * guide
        else
            guide *= CHOPOFF.getExt3(chrom, length(chrom), last(range) + 1, distance)
        end
    end
    if motif.extends5 && is_antisense
        return complement(guide)
    elseif motif.extends5
        return reverse(guide)
    elseif is_antisense
        return reverse_complement(guide)
    end
    return guide
end

strip_gaps(value) = replace(String(value), "-" => "")

function reported_alignment_cost(guide::String, reference::String)
    length(guide) == length(reference) || return typemax(Int)
    cost = 0
    for (guide_base, ref_base) in zip(guide, reference)
        if guide_base == '-' || ref_base == '-'
            guide_base == ref_base && return typemax(Int)
            cost += 1
        else
            g = LongDNA{4}(string(guide_base))[1]
            r = LongDNA{4}(string(ref_base))[1]
            iscompatible(g, r) || (cost += 1)
        end
    end
    return cost
end

function row_inspector(cache::ReferenceCache, motif::Motif, distance::Int)
    return function(row)
        try
            chrom = cached_chrom(cache, String(row.chromosome))
            range = candidate_range(row, motif)
            first(range) >= 1 && last(range) <= length(chrom) ||
                return (ambiguity_count = 0, valid = false, optimal = false)
            is_antisense = String(row.strand) == "-"
            pattern = is_antisense ? motif.rve : motif.fwd
            isempty(pattern) &&
                return (ambiguity_count = 0, valid = false, optimal = false)
            window = chrom[range]
            ambiguity = count(isambiguous, window)
            motif_valid = all(iscompatible(pattern[idx], window[idx])
                for idx in eachindex(pattern))
            ot = normalized_candidate(chrom, range, motif, is_antisense, distance)
            guide = LongDNA{4}(String(row.guide))
            query = motif.extends5 ? reverse(guide) : guide
            oracle = CHOPOFF.align(query, ot, distance, iscompatible)
            output_ot = motif.extends5 ? reverse(ot) : ot
            aligned_guide = String(row.alignment_guide)
            aligned_ref = String(row.alignment_reference)
            source_valid = strip_gaps(aligned_guide) == String(row.guide) &&
                occursin(strip_gaps(aligned_ref), string(output_ot))
            cost_valid = reported_alignment_cost(aligned_guide, aligned_ref) ==
                Int(row.distance)
            optimal = Int(row.distance) == oracle.dist
            return (ambiguity_count = ambiguity,
                valid = motif_valid && ambiguity <= 3 && source_valid && cost_valid,
                optimal = optimal)
        catch
            return (ambiguity_count = 0, valid = false, optimal = false)
        end
    end
end

function filtered_prefix_path(
    out::String, case::GenericParityCase, distance::Int,
    by_id::Dict{String, GenericParityCase})

    base = by_id[case.baseline_id]
    source = joinpath(case_dir(out, base), "prefix_d$(distance).csv")
    case.baseline_strand === nothing && return source
    filtered = joinpath(case_dir(out, case), "prefix_d$(distance)_filtered.csv")
    table = DataFrame(CSV.File(source))
    CSV.write(filtered, table[String.(table.strand) .== case.baseline_strand, :])
    return filtered
end

function normalize_counts(path::String)
    table = DataFrame(CSV.File(path))
    table.guide = String.(table.guide)
    for column in names(table)
        startswith(column, "D") && (table[!, column] = Int.(table[!, column]))
    end
    table.complete = Bool.(table.complete)
    return table
end

function reason_count(result, reason::String)
    return count(==(reason), result.differences.reason)
end

function run_compare(out::String)
    all_cases = cases(out)
    by_id = Dict(case.id => case for case in all_cases)
    summary_rows = NamedTuple[]
    difference_tables = DataFrame[]
    count_rows = NamedTuple[]
    passed = true
    for case in all_cases
        guides = string.(load_guides(out, case))
        for reference in case.reference_paths, distance in case.distances,
                ambiguity in case.ambiguities
            label = reference_label(reference)
            scan_detail = joinpath(case_dir(out, case),
                "scan_$(label)_d$(distance)_a$(ambiguity)_detail.csv")
            prefix = filtered_prefix_path(out, case, distance, by_id)
            difference_path = joinpath(case_dir(out, case),
                "differences_$(label)_d$(distance)_a$(ambiguity).csv")
            inspector = row_inspector(ReferenceCache(reference),
                setdist(case.motif, distance), distance)
            result = classify_semantic_parity(
                scan_detail, prefix, ambiguity, inspector;
                output_path = difference_path)
            if nrow(result.differences) > 0
                table = copy(result.differences)
                insertcols!(table, 1,
                    :case_id => fill(case.id, nrow(table)),
                    :reference => fill(reference, nrow(table)),
                    :requested_distance => fill(distance, nrow(table)),
                    :requested_ambiguity => fill(ambiguity, nrow(table)))
                push!(difference_tables, table)
            end
            push!(summary_rows, (
                case_id = case.id, reference = reference,
                distance = distance, ambiguity = ambiguity,
                exact = result.exact, scan_only = result.scan_only,
                prefix_only = result.prefix_only, accepted = result.accepted,
                rejected = result.rejected,
                above_limit = reason_count(result, "baseline_above_ambiguity_limit"),
                legacy_duplicate = reason_count(result, "legacy_duplicate"),
                legacy_false_negative = reason_count(result, "legacy_false_negative"),
                legacy_false_positive = reason_count(result, "legacy_false_positive"),
                traceback_tie = reason_count(result, "optimal_traceback_tie"),
                passed = result.passed,
            ))
            passed &= result.passed

            expected = expected_counts_from_detail(scan_detail, guides, distance)
            count_path = joinpath(case_dir(out, case),
                "scan_$(label)_d$(distance)_a$(ambiguity)_counts.csv")
            actual = normalize_counts(count_path)
            count_ok = actual == expected
            push!(count_rows, (case_id = case.id, reference = reference,
                distance = distance, ambiguity = ambiguity, passed = count_ok,
                expected_rows = nrow(expected), actual_rows = nrow(actual)))
            passed &= count_ok
            println("compare case=$(case.id) ref=$label d=$distance a=$ambiguity parity=$(result.passed ? "PASS" : "FAIL") counts=$(count_ok ? "PASS" : "FAIL")")
            flush(stdout)
        end
    end
    CSV.write(joinpath(out, "parity_summary.csv"), DataFrame(summary_rows))
    CSV.write(joinpath(out, "count_parity.csv"), DataFrame(count_rows))
    differences = isempty(difference_tables) ? DataFrame() :
        reduce((left, right) -> vcat(left, right; cols = :union), difference_tables)
    CSV.write(joinpath(out, "parity_differences.csv"), differences)
    passed || error("Generic human parity matrix failed; see $out")
    println("Generic human parity matrix PASS: $out")
end

function run_child(stage::String, threads::Int, out::String)
    script = abspath(@__FILE__)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    env = copy(ENV)
    env["CHOPOFF_GENERIC_PARITY_STAGE"] = stage
    env["CHOPOFF_GENERIC_PARITY_OUT"] = out
    run(setenv(`$julia --project=$ROOT_DIR --threads=$threads $script`, env))
end

function main()
    stage = lowercase(strip(get(ENV, "CHOPOFF_GENERIC_PARITY_STAGE", "all")))
    out = output_dir()
    if stage == "all"
        run_prepare(out)
        if mode() == :qualification
            search_threads = parse(Int, get(ENV, "CHOPOFF_GENERIC_PARITY_THREADS", "24"))
            run_child("build", 1, out)
            run_child("search", search_threads, out)
            run_child("compare", 1, out)
        else
            run_build(out)
            run_search(out)
            run_compare(out)
        end
    elseif stage == "prepare"
        run_prepare(out)
    elseif stage == "build"
        run_build(out)
    elseif stage == "search"
        run_search(out)
    elseif stage == "compare"
        run_compare(out)
    else
        error("CHOPOFF_GENERIC_PARITY_STAGE must be all, prepare, build, search, or compare")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
