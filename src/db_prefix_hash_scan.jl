struct PrefixScanGeometry{Kind, M}
    guide_bases::Int
    pam_bases::Int
    prefix_bases::Int
    distance::Int
    matcher::M
end

PrefixScanGeometry{Kind}(
    guide_bases::Int, pam_bases::Int, prefix_bases::Int, distance::Int) where Kind =
    PrefixScanGeometry{Kind, Nothing}(
        guide_bases, pam_bases, prefix_bases, distance, nothing)

const CAS9_D3_PREFIX_SCAN_GEOMETRY =
    PrefixScanGeometry{:cas9}(20, 3, 16, 3)
const CAS12A_D3_PREFIX_SCAN_GEOMETRY =
    PrefixScanGeometry{:cas12a}(21, 4, 16, 3)

@inline prefix_scan_kind(::PrefixScanGeometry{Kind}) where Kind = Kind

function matches_prefix_scan_motif(motif::Motif, name::String)
    template = Motif(name)
    return motif.fwd == template.fwd &&
        motif.rve == template.rve &&
        motif.pam_loci_fwd == template.pam_loci_fwd &&
        motif.pam_loci_rve == template.pam_loci_rve &&
        motif.extends5 == template.extends5
end

function resolve_prefix_scan_geometry(
    motif::Motif, distance::Int, hash_len::Int)

    distance in 0:4 || return nothing
    1 <= hash_len <= 16 || return nothing
    hash_len == 16 && matches_prefix_scan_motif(motif, "Cas9") &&
        return PrefixScanGeometry{:cas9}(20, 3, 16, distance)
    hash_len == 16 && matches_prefix_scan_motif(motif, "Cas12a") &&
        return PrefixScanGeometry{:cas12a}(21, 4, 16, distance)
    return resolve_generic_prefix_scan_geometry(motif, distance, hash_len)
end

@inline prefix_scan_candidate_bases(geometry::PrefixScanGeometry) =
    geometry.guide_bases + geometry.pam_bases

@inline prefix_scan_candidate_last_offset(geometry::PrefixScanGeometry) =
    prefix_scan_candidate_bases(geometry) - 1

mutable struct PrefixHashScanStats
    motif_candidates::Int
    ambiguous_prefixes::Int
    prefix_hits::Int
    guide_pairs::Int
    alignment_calls::Int
    distance_calls::Int
    traceback_calls::Int
    emitted_rows::Int
    path_rows::Int
    query_hashes::Int
    bruteforce_guide_pairs::Int
    metadata_ns::UInt64
    query_build_ns::UInt64
    path_load_ns::UInt64
    query_hash_ns::UInt64
    query_format_ns::UInt64
    query_fold_ns::UInt64
    query_dedup_ns::UInt64
    query_insert_ns::UInt64
    query_lookup_ns::UInt64
    chrom_load_ns::UInt64
    record_io_ns::UInt64
    sequence_convert_ns::UInt64
    findguides_ns::UInt64
    candidate_prefix_ns::UInt64
    candidate_hash_ns::UInt64
    candidate_materialize_ns::UInt64
    align_ns::UInt64
    emit_ns::UInt64
    scan_ns::UInt64
    verify_ns::UInt64
    work_items_total::Int
    work_items_claimed::Int
    retired_guides::Int
    path_source::Symbol
    query_variant::Symbol
    scan_backend::Symbol
    simd_backend::Symbol
end

# Built from the field types rather than a positional list of 34 zeros and four
# symbols, so adding a counter to the struct needs no edit here or in reset!.
PrefixHashScanStats() = PrefixHashScanStats(
    (T === Symbol ? :none : zero(T)
     for T in fieldtypes(PrefixHashScanStats))...)

@inline prefix_hash_scan_timer(::Nothing) = UInt64(0)
@inline prefix_hash_scan_timer(::PrefixHashScanStats) = time_ns()
@inline prefix_hash_scan_worker_stats(::Nothing) = nothing
@inline prefix_hash_scan_worker_stats(::PrefixHashScanStats) = PrefixHashScanStats()

function reset!(stats::PrefixHashScanStats)
    for (name, T) in zip(
            fieldnames(PrefixHashScanStats),
            fieldtypes(PrefixHashScanStats))
        setfield!(stats, name, T === Symbol ? :none : zero(T))
    end
    return stats
end

struct PrefixHashScanBitmaskQuery{H <: Unsigned}
    masks::Dict{H, UInt64}
end

struct PrefixHashScanDirectory
    offsets::Vector{UInt32}
    suffixes::Vector{UInt16}
    masks::Vector{UInt64}
    hash_len::UInt8
    bucket_bases::UInt8
end

struct PrefixHashScanPrefilteredDirectory{D}
    directory::D
    presence::Vector{UInt64}
    prefix_bits::UInt8
end

struct PrefixHashScanHit
    start::Int
    mask::UInt64
end

struct PrefixHashScanLookupScratch
    plus_candidates::Vector{UInt64}
    minus_candidates::Vector{UInt64}
    plus_radix::Vector{UInt64}
    minus_radix::Vector{UInt64}
    radix_counts::Vector{Int}
end

PrefixHashScanLookupScratch() = PrefixHashScanLookupScratch(
    UInt64[], UInt64[], UInt64[], UInt64[], zeros(Int, 2048))

struct PrefixHashScanVerifiedHit
    guide_idx::Int
    pos::Int
    dist::Int
    is_antisense::Bool
    aln_guide::String
    aln_ref::String
end

struct PrefixHashScanChromResult{S}
    plus::Vector{PrefixHashScanVerifiedHit}
    minus::Vector{PrefixHashScanVerifiedHit}
    stats::S
end

struct PrefixHashScanCountResult{S}
    counts::Matrix{Int}
    stats::S
end

struct PrefixHashScanChunkWork
    chrom_idx::Int
    core_first::Int
    core_last::Int
end

struct PrefixHashScanMyersProfile
    eq_by_iupac::NTuple{16, UInt64}
    length::UInt8
    final_bit::UInt64
end

mutable struct PrefixHashScanEarlyStopState
    limits::Vector{Int}
    counts::Matrix{Int}
    active_mask::Threads.Atomic{UInt64}
    retired_guides::Threads.Atomic{Int}
    work_items_total::Int
    work_items_claimed::Threads.Atomic{Int}
    merge_lock::ReentrantLock
end

function PrefixHashScanEarlyStopState(
    guide_count::Int, limits::Vector{Int})

    active_mask = guide_count == 64 ? typemax(UInt64) :
        (UInt64(1) << guide_count) - UInt64(1)
    return PrefixHashScanEarlyStopState(
        copy(limits),
        zeros(Int, guide_count, length(limits)),
        Threads.Atomic{UInt64}(active_mask),
        Threads.Atomic{Int}(0),
        0,
        Threads.Atomic{Int}(0),
        ReentrantLock(),
    )
end

@inline prefix_hash_scan_active_mask(::Nothing) = typemax(UInt64)
@inline prefix_hash_scan_active_mask(state::PrefixHashScanEarlyStopState) =
    state.active_mask[]

@inline function retire_prefix_hash_scan_guide!(
    state::PrefixHashScanEarlyStopState, guide_idx::Int)

    bit = UInt64(1) << (guide_idx - 1)
    old = Threads.atomic_and!(state.active_mask, ~bit)
    old & bit != 0 && Threads.atomic_add!(state.retired_guides, 1)
    return nothing
end

@inline reserve_prefix_hash_scan_hit!(
    ::Nothing, ::Union{Nothing, Matrix{Int}}, ::Int, ::Int) = true

@inline function reserve_prefix_hash_scan_hit!(
    state::PrefixHashScanEarlyStopState,
    chunk_counts::Union{Nothing, Matrix{Int}},
    guide_idx::Int,
    dist::Int)

    bit = UInt64(1) << (guide_idx - 1)
    state.active_mask[] & bit == 0 && return false
    counts = chunk_counts::Matrix{Int}
    counts[guide_idx, dist + 1] += 1
    return counts[guide_idx, dist + 1] <= state.limits[dist + 1]
end

function merge_prefix_hash_scan_chunk_counts!(
    state::PrefixHashScanEarlyStopState, counts::Matrix{Int})

    lock(state.merge_lock) do
        for guide_idx in axes(counts, 1), dist_idx in axes(counts, 2)
            count = counts[guide_idx, dist_idx]
            count == 0 && continue
            old = state.counts[guide_idx, dist_idx]
            state.counts[guide_idx, dist_idx] += count
            old <= state.limits[dist_idx] < old + count &&
                retire_prefix_hash_scan_guide!(state, guide_idx)
        end
    end
    return nothing
end

function prefix_hash_scan_early_stop_counts(state::PrefixHashScanEarlyStopState)
    return copy(state.counts)
end

function merge_prefix_hash_scan_early_stop_stats!(
    stats::Union{Nothing, PrefixHashScanStats},
    state::Union{Nothing, PrefixHashScanEarlyStopState})

    stats === nothing && return nothing
    state === nothing && return nothing
    stats.work_items_total += state.work_items_total
    stats.work_items_claimed += state.work_items_claimed[]
    stats.retired_guides += state.retired_guides[]
    return nothing
end


include("prefix_hash_scan/query.jl")
include("prefix_hash_scan/verification.jl")
include("prefix_hash_scan/kernel_common.jl")
include("prefix_hash_scan/cas9.jl")
include("prefix_hash_scan/cas12a.jl")
include("prefix_hash_scan/generic.jl")
include("prefix_hash_scan/twobit.jl")
include("prefix_hash_scan/streaming.jl")


# One verified candidate, aligned and written out. The legacy engine reaches
# this from three different guide iterations (bitmask, bruteforce, candidate
# list) that used to carry three byte-identical copies of this body.
function emit_prefix_hash_scan_legacy_hit!(
    out, guide_idx::Int, ot, pos, chrom_name, strand,
    guides, guides_, distance::Int, motif::Motif,
    seen, es_acc, is_es, early_stopping, stats)

    is_es[guide_idx] && return nothing
    if stats !== nothing
        stats.alignment_calls += 1
        stats.traceback_calls += 1
    end
    align_start = prefix_hash_scan_timer(stats)
    aln = align(guides_[guide_idx], ot, distance, iscompatible)
    stats === nothing || (stats.align_ns += time_ns() - align_start)
    aln.dist > distance && return nothing

    if motif.extends5
        aln_guide = reverse(aln.guide)
        aln_ref = reverse(aln.ref)
    else
        aln_guide = aln.guide
        aln_ref = aln.ref
    end

    key = (
        string(guides[guide_idx]),
        aln.dist,
        chrom_name,
        Int(pos),
        strand,
        aln_guide,
        aln_ref,
    )
    key in seen[guide_idx] && return nothing
    push!(seen[guide_idx], key)
    reserve_prefix_hash_scan_detail_hit!(
        es_acc, is_es, guide_idx, aln.dist, early_stopping) || return nothing

    emit_start = prefix_hash_scan_timer(stats)
    print(out, guides[guide_idx], ",", aln_guide, ",", aln_ref, ",",
        aln.dist, ",", chrom_name, ",", pos, ",", strand, "\n")
    if stats !== nothing
        stats.emit_ns += time_ns() - emit_start
        stats.emitted_rows += 1
    end
    return nothing
end


# Shared by both `search_prefixHashScan` methods so the two entry points cannot
# drift apart on what they accept. The motif/distance check is load bearing:
# when no precomputed path asset matches, `load_prefix_hash_scan_paths` falls
# back to `build_PathTemplates`, which enumerates paths for `motif.distance`
# and ignores the requested `distance`. A motif carrying a smaller distance
# would silently yield too few paths and under-report off-targets.
function validate_prefix_hash_scan_query(
    guides::Vector{LongDNA{4}},
    motif::Motif,
    distance::Int,
    early_stopping::Vector{Int},
    output::Symbol)

    distance in 0:4 ||
        error("search_prefixHashScan supports distances 0 through 4.")
    isempty(guides) &&
        error("search_prefixHashScan requires at least one guide.")
    distance <= motif.distance ||
        error("For this motif maximum distance is " * string(motif.distance))
    motif.ambig_max in 0:3 ||
        error("search_prefixHashScan supports ambig_max values 0 through 3.")
    guide_bases = length_noPAM(motif)
    all(==(guide_bases), length.(guides)) ||
        error("Guide queries are not of the correct length to use with this Motif: " * string(motif))
    any(isambig.(guides)) &&
        error("search_prefixHashScan does not support ambiguous query guides.")
    length(early_stopping) == (distance + 1) ||
        error("Specify one early stopping condition for each distance from 0 to $distance.")
    all(>=(0), early_stopping) ||
        error("early_stopping limits must be nonnegative.")
    output in (:detail, :counts) ||
        error("output must be :detail or :counts.")
    return nothing
end


"""
```
search_prefixHashScan(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif,
    output_file::String;
    distance::Int = 3,
    hash_len::Int = min(length_noPAM(motif) - distance, 16),
    early_stopping::Vector{Int} = fill(1_000_000, distance + 1),
    scan_threads::Int = Threads.nthreads(),
    simd_backend::Symbol = :auto,
    output::Symbol = :detail,
    verbose::Bool = false,
    kwargs...)
```

Low-level, fully configurable form of the indexless prefixHashScan search. It
takes an explicit `motif` and exposes the scan tuning knobs. Prefer the
three-argument method for ordinary use; this one is for benchmarking and for
callers that need to pin a specific engine.

Reuses prefixHashDB's symbolic edit-distance prefix paths. Standard Cas9 and
Cas12a distance-0-through-4/prefix-16 motifs use separate specialized scan
kernels. Other compatible 16-base-prefix motifs use a typed generic scan
kernel; configurations outside its size envelope use the exact legacy engine.
Reference windows may contain zero through three IUPAC ambiguity symbols, as
specified by `motif.ambig_max`. Query guides must be unambiguous.

`distance` must not exceed `motif.distance`: when no precomputed path asset
matches, the symbolic paths are generated from the motif, so a motif carrying a
smaller distance would silently produce an incomplete result. This is checked.

# Arguments

`motif` - Motif to search for. `distance` must be at most `motif.distance`, and
every guide must be `length_noPAM(motif)` bases long.

`output_file` - Path for the comma separated output table, so a `.csv`
extension is preferred.

`distance` - Maximum levenshtein distance (insertions, deletions, mismatches),
0 through 4.

`hash_len` - Symbolic prefix length in `1:16`. The default keeps
`guide length - distance` bases, capped at 16. Only 16 reaches the specialized
and typed generic SIMD kernels; shorter prefixes use the exact legacy engine.

`early_stopping` - One limit per distance, from 0 through `distance`. A guide
stops accumulating hits at a given distance once its limit is reached.

`scan_threads` - Threads used for the genome scan.

`simd_backend` - `:auto`, `:avx2`, or `:avx512` for raw SIMD scans. `:auto`
resolves from the CPU and the scan geometry.

`output` - `:detail` writes aligned loci with the usual CHOPOFF columns;
`:counts` writes one row per unique guide as `guide,D0,...,Dk,complete`,
without traceback. Note this is not the same schema as `summarize_offtargets`,
which has no `complete` column.

`verbose` - Print a summary of the resolved engine and timings.

## Tuning arguments

These select between equivalent implementations. They exist for benchmarking
and are not needed for correct results; every combination is expected to
produce identical output. They are not covered by semantic versioning.

`query_variant` - `:auto`, `:baseline`, `:columnwise`, `:bitmask64`, or
`:bruteforce`. `:bruteforce` skips the prefix-hash filter and verifies every
motif candidate, so it is the reference for checking the filter has no false
negatives. `:bitmask64` supports at most 64 guides.

`scan_backend` - Pins the scan engine, e.g. `:legacy`, `:fused_directory`,
`:streaming_fasta_simd`, `:streaming_2bit_simd`. `:auto` picks the fastest
applicable one.

`query_build_backend` - `:auto`, `:serial`, or `:parallel` query construction.

`lookup_variant` - `:auto`, `:inline`, or `:bucketed` directory lookup.

`verify_variant` - `:auto`, `:align`, `:distance_first`, or `:myers_raw`.

`bucket_bases` - Directory bucket width in bases.

`prefilter_bits` - Presence-prefilter width: 0 (disabled), 22, 24, or 26.

`stream_chunk_bases` - Reference chunk size for the streaming scheduler, at
least 64.

`stats` - Optional `CHOPOFF.PrefixHashScanStats` to fill with counters and
per-stage timings.

# Examples
```julia
$(make_prefix_hash_scan_example_doc())
```
"""
function search_prefixHashScan(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    motif::Motif,
    output_file::String;
    distance::Int = 3,
    hash_len::Int = min(length_noPAM(motif) - distance, 16),
    early_stopping::Vector{Int} = fill(1_000_000, distance + 1),
    query_variant::Symbol = prefix_hash_scan_query_variant(),
    scan_backend::Symbol = :auto,
    simd_backend::Symbol = :auto,
    bucket_bases::Int = 11,
    scan_threads::Int = Threads.nthreads(),
    stream_chunk_bases::Int = 2 * 1024 * 1024,
    prefilter_bits::Int = 26,
    query_build_backend::Symbol = :auto,
    lookup_variant::Symbol = :auto,
    verify_variant::Symbol = :auto,
    output::Symbol = :detail,
    verbose::Bool = false,
    stats::Union{Nothing, PrefixHashScanStats} = nothing,
    _append_output::Bool = false,
    _collect_counts::Bool = false,
    _paths = nothing)

    validate_prefix_hash_scan_query(
        guides, motif, distance, early_stopping, output)
    if hash_len < 1 || hash_len > 16
        error("hash_len must be in 1:16.")
    end
    scan_threads >= 1 || error("scan_threads must be positive.")
    stream_chunk_bases >= 64 ||
        error("stream_chunk_bases must be at least 64.")
    prefilter_bits in (0, 22, 24, 26) ||
        error("prefilter_bits must be 0, 22, 24, or 26.")
    query_build_backend in (:auto, :serial, :parallel) ||
        error("query_build_backend must be :auto, :serial, or :parallel.")
    lookup_variant in (:auto, :inline, :bucketed) ||
        error("lookup_variant must be :auto, :inline, or :bucketed.")
    verify_variant in (:auto, :align, :distance_first, :myers_raw) ||
        error("verify_variant must be :auto, :align, :distance_first, or :myers_raw.")

    if stats !== nothing
        reset!(stats)
    end

    metadata_start = prefix_hash_scan_timer(stats)
    dbi, reference_lengths = prefix_hash_scan_dbinfo(genome_path, motif)
    if stats !== nothing
        stats.metadata_ns += time_ns() - metadata_start
    end
    hash_type = smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
    resolved_query_variant =
        resolve_prefix_hash_scan_query_variant(query_variant, length(guides))
    scan_backend in (
        :auto, :legacy, :fused_dict, :fused_directory, :fused_fasta_simd,
        :streaming_fasta_simd, :streaming_fasta_simd_fused,
        :streaming_2bit_simd) ||
        error("scan_backend must be :auto, :legacy, :fused_dict, :fused_directory, :fused_fasta_simd, :streaming_fasta_simd, :streaming_fasta_simd_fused, or :streaming_2bit_simd.")
    geometry = resolve_prefix_scan_geometry(motif, distance, hash_len)
    scan_kind = geometry === nothing ? :none : prefix_scan_kind(geometry)
    resolved_simd_backend = resolve_prefix_hash_scan_simd_backend(
        simd_backend; scan_kind)
    supports_fused = geometry !== nothing &&
        resolved_query_variant == :bitmask64
    supports_raw_simd = supports_fused &&
        hash_len == (geometry::PrefixScanGeometry).prefix_bases &&
        resolved_simd_backend != :none
    supports_fasta_simd = supports_raw_simd && dbi.gi.is_fa
    supports_twobit_simd = supports_raw_simd && !dbi.gi.is_fa
    resolved_scan_backend = if scan_backend != :auto
        scan_backend
    elseif supports_fasta_simd
        :streaming_fasta_simd
    elseif supports_twobit_simd
        :streaming_2bit_simd
    elseif supports_fused
        :fused_directory
    else
        :legacy
    end
    if output == :counts
        resolved_scan_backend == :streaming_fasta_simd_fused &&
            (resolved_scan_backend = :streaming_fasta_simd)
        resolved_scan_backend in (
            :streaming_fasta_simd, :streaming_2bit_simd) ||
            (resolved_scan_backend = :legacy)
    end
    if resolved_scan_backend != :legacy && !supports_fused
        error("Fused scan backends require an optimized distance-0-through-4, 16-base-prefix geometry and at most 64 guides.")
    end
    if resolved_scan_backend in (
            :fused_fasta_simd, :streaming_fasta_simd,
            :streaming_fasta_simd_fused) &&
            !supports_fasta_simd
        resolved_scan_backend = :fused_directory
    end
    if resolved_scan_backend == :streaming_2bit_simd && !supports_twobit_simd
        resolved_scan_backend = :fused_directory
    end
    # Derive the early-stopping state only once the backend is final. The
    # downgrades above used to run afterwards, so a requested streaming backend
    # that fell back to :fused_directory left a live state attached to an
    # engine that never consumes it. Only the two streaming backends thread
    # early stopping through their chunk workers; :streaming_fasta_simd_fused
    # does not (stream_prefix_hash_scan_chunk's :fused branch takes no
    # early_stop_state), so it is deliberately excluded here.
    early_stop_state = all(==(typemax(Int)), early_stopping) ||
        !(resolved_scan_backend in (
            :streaming_fasta_simd, :streaming_2bit_simd)) ? nothing :
        PrefixHashScanEarlyStopState(length(guides), early_stopping)
    uses_raw_simd = resolved_scan_backend in (
        :fused_fasta_simd, :streaming_fasta_simd,
        :streaming_fasta_simd_fused, :streaming_2bit_simd)
    if simd_backend != :auto && !uses_raw_simd
        error("simd_backend=$(repr(simd_backend)) requires an applicable raw SIMD scan backend.")
    end
    effective_simd_backend = uses_raw_simd ? resolved_simd_backend : :none
    resolved_lookup_variant = if lookup_variant == :auto
        resolved_scan_backend in (:streaming_fasta_simd, :streaming_2bit_simd) &&
            prefilter_bits != 0 && bucket_bases == 11 &&
            stream_chunk_bases + prefix_scan_candidate_last_offset(
                geometry::PrefixScanGeometry) + distance <= typemax(UInt32) ?
            :bucketed : :inline
    else
        lookup_variant
    end
    if resolved_lookup_variant == :bucketed
        resolved_scan_backend in (:streaming_fasta_simd, :streaming_2bit_simd) ||
            error("lookup_variant=:bucketed requires a buffered streaming SIMD backend.")
        prefilter_bits != 0 ||
            error("lookup_variant=:bucketed requires a nonzero prefilter.")
        bucket_bases == 11 ||
            error("lookup_variant=:bucketed currently requires bucket_bases=11.")
        stream_chunk_bases + prefix_scan_candidate_last_offset(
            geometry::PrefixScanGeometry) + distance <= typemax(UInt32) ||
            error("lookup_variant=:bucketed requires chunks smaller than 4 GiB.")
    end
    if verify_variant == :myers_raw &&
            !(resolved_scan_backend in (
                :fused_fasta_simd, :streaming_fasta_simd,
                :streaming_fasta_simd_fused, :streaming_2bit_simd))
        error("verify_variant=:myers_raw requires a raw SIMD backend.")
    end
    if resolved_scan_backend in (
            :streaming_fasta_simd, :streaming_fasta_simd_fused,
            :streaming_2bit_simd) &&
            !(verify_variant in (:auto, :myers_raw))
        error("Streaming SIMD backends require verify_variant=:auto or :myers_raw.")
    end
    if verbose
        worker_count = min(scan_threads, Threads.nthreads(), length(guides))
        resolved_query_build_backend = resolve_prefix_hash_scan_query_build_backend(
            query_build_backend, worker_count)
        scheduler = resolved_scan_backend in (
            :streaming_fasta_simd, :streaming_fasta_simd_fused,
            :streaming_2bit_simd) ?
            :chunk : :record
        @info(
            "prefixHashScan execution",
            scan_geometry = geometry === nothing ? :generic : prefix_scan_kind(geometry),
            reference_format = dbi.gi.is_fa ? :fasta : :twobit,
            scan_backend = resolved_scan_backend,
            simd_backend = effective_simd_backend,
            lookup_variant = resolved_lookup_variant,
            query_build_backend = resolved_query_build_backend,
            scheduler = scheduler,
            scan_threads = scan_threads,
            chunk_bases = stream_chunk_bases,
        )
    end

    query_start = prefix_hash_scan_timer(stats)
    if resolved_scan_backend in (
            :fused_directory, :fused_fasta_simd, :streaming_fasta_simd,
            :streaming_fasta_simd_fused, :streaming_2bit_simd)
        query, guides_ = build_prefix_hash_scan_compact_query(
            guides,
            motif,
            distance,
            hash_len,
            stats;
            bucket_bases = bucket_bases,
            prefilter_bits = resolved_scan_backend in (
                :fused_fasta_simd, :streaming_fasta_simd,
                :streaming_fasta_simd_fused, :streaming_2bit_simd) ?
                prefilter_bits : 0,
            query_build_backend = query_build_backend,
            query_threads = scan_threads,
            paths = _paths,
        )
    elseif resolved_query_variant == :bruteforce
        query = nothing
        guides_ = oriented_prefix_hash_scan_guides(guides, motif)
        if stats !== nothing
            stats.query_variant = :bruteforce
        end
    else
        query, guides_ = build_prefix_hash_scan_map(
            guides,
            motif,
            distance,
            hash_len,
            hash_type,
            stats;
            query_variant = resolved_query_variant,
            paths = _paths,
        )
    end
    if stats !== nothing
        stats.scan_backend = resolved_scan_backend
        stats.simd_backend = effective_simd_backend
        stats.query_build_ns += time_ns() - query_start
    end

    es_acc = zeros(Int, length(guides), length(early_stopping))
    is_es = falses(length(guides))
    seen = [Set{Tuple{String, Int, String, Int, String, String, String}}() for _ in guides]
    candidate_guides = Int[]
    use_bruteforce_query = resolved_query_variant == :bruteforce
    use_bitmask_query = query isa PrefixHashScanBitmaskQuery
    use_direct_specialized_hash = geometry !== nothing
    use_fused_scan = resolved_scan_backend in (
        :fused_dict, :fused_directory, :fused_fasta_simd,
        :streaming_fasta_simd, :streaming_fasta_simd_fused,
        :streaming_2bit_simd)
    use_myers_raw = verify_variant == :myers_raw ||
        (verify_variant == :auto && resolved_scan_backend in (
            :fused_fasta_simd, :streaming_fasta_simd,
            :streaming_fasta_simd_fused, :streaming_2bit_simd))
    use_distance_first = verify_variant == :distance_first ||
        (verify_variant == :auto && use_fused_scan && !use_myers_raw)
    myers_profiles = use_myers_raw || output == :counts &&
            resolved_scan_backend in (
                :streaming_fasta_simd, :streaming_2bit_simd) ?
        build_prefix_hash_scan_myers_profiles(guides_) : nothing

    if output == :counts
        scan_start = prefix_hash_scan_timer(stats)
        counts = if resolved_scan_backend in (
                :streaming_fasta_simd, :streaming_2bit_simd)
            chunk_results = stream_prefix_hash_scan_counts(
                geometry::PrefixScanGeometry,
                dbi.gi.filepath,
                reference_lengths,
                query,
                dbi,
                guides_,
                myers_profiles::Vector{PrefixHashScanMyersProfile},
                distance,
                stream_chunk_bases,
                scan_threads,
                Val(resolved_lookup_variant == :bucketed ?
                    :bucketed_reuse : :buffered_reuse),
                stats,
                early_stop_state,
                simd_backend = Val(effective_simd_backend),
            )
            merged = early_stop_state === nothing ?
                zeros(Int, length(guides), distance + 1) :
                prefix_hash_scan_early_stop_counts(early_stop_state)
            for result_ in chunk_results
                result_ === nothing && continue
                early_stop_state === nothing && (merged .+= result_.counts)
                stats === nothing || merge_prefix_hash_scan_worker_stats!(
                    stats, result_.stats)
            end
            merged
        else
            search_prefix_hash_scan_legacy_counts(
                dbi, reference_lengths, query, guides_, geometry, distance,
                hash_len, hash_type,
                use_bruteforce_query, use_bitmask_query,
                use_direct_specialized_hash, stats)
        end
        merge_prefix_hash_scan_early_stop_stats!(stats, early_stop_state)
        stats === nothing || (stats.scan_ns += time_ns() - scan_start)
        _collect_counts && return counts
        write_prefix_hash_scan_counts(
            output_file, guides, counts, early_stopping, stats)
        return nothing
    end

    mkpath(dirname(output_file))
    open(output_file, _append_output ? "a" : "w") do out
        _append_output ||
            write(out, "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand\n")

        if resolved_scan_backend in (
                :streaming_fasta_simd, :streaming_fasta_simd_fused,
                :streaming_2bit_simd)
            scan_start = prefix_hash_scan_timer(stats)
            chunk_results, chrom_chunk_ranges = stream_prefix_hash_scan(
                geometry::PrefixScanGeometry,
                dbi.gi.filepath,
                reference_lengths,
                query,
                dbi,
                guides_,
                myers_profiles,
                distance,
                stream_chunk_bases,
                scan_threads,
                Val(resolved_scan_backend == :streaming_fasta_simd_fused ?
                    :fused : (resolved_lookup_variant == :bucketed ?
                        :bucketed_reuse : :buffered_reuse)),
                stats,
                Val(:chunk),
                early_stop_state,
                simd_backend = Val(effective_simd_backend),
            )
            for chrom_idx in eachindex(chrom_chunk_ranges)
                chrom_name = dbi.gi.chrom[chrom_idx]
                chunk_range = chrom_chunk_ranges[chrom_idx]
                if stats !== nothing
                    for chunk_idx in chunk_range
                        result_ = chunk_results[chunk_idx]
                        result_ === nothing && continue
                        merge_prefix_hash_scan_worker_stats!(stats, result_.stats)
                    end
                end
                for strand in (:plus, :minus)
                    for chunk_idx in chunk_range
                        result_ = chunk_results[chunk_idx]
                        result_ === nothing && continue
                        hits = getfield(result_, strand)
                        for hit in hits
                            commit_prefix_hash_scan_verified!(
                                out,
                                hit,
                                guides[hit.guide_idx],
                                chrom_name,
                                early_stopping,
                                es_acc,
                                is_es,
                                seen,
                                stats,
                                prelimited = early_stop_state !== nothing,
                            )
                        end
                    end
                end
            end
            if stats !== nothing
                stats.scan_ns += time_ns() - scan_start
            end
            return nothing
        end

        ref = open(dbi.gi.filepath, "r")
        try
            reader = dbi.gi.is_fa ? FASTA.Reader(ref; copy = false) : TwoBit.Reader(ref)
            records = dbi.gi.is_fa ?
                PrefixHashScanFASTARecords(reader, dbi.gi.chrom, reference_lengths) :
                PrefixHashScanIndexedRecords(reader, dbi.gi.chrom)
            scan_start = time_ns()
            for (chrom_name, record, record_io_ns) in records
                if resolved_scan_backend == :fused_fasta_simd
                    raw_part = FASTX.seq_data_part(record, 1:FASTX.seqsize(record))
                    raw = @view record.data[raw_part]
                    if stats !== nothing
                        stats.record_io_ns += record_io_ns
                        stats.chrom_load_ns += record_io_ns
                    end
                    plus_hits, minus_hits = scan_prefix_hits_raw(
                        geometry::PrefixScanGeometry, raw, dbi, query, stats;
                        scan_threads = scan_threads,
                        simd_backend = Val(effective_simd_backend))
                    for (is_antisense, hits) in ((false, plus_hits), (true, minus_hits))
                        if stats !== nothing
                            stats.prefix_hits += length(hits)
                            stats.guide_pairs +=
                                sum(hit -> count_ones(hit.mask), hits; init = 0)
                        end
                        for hit in hits
                            verify_prefix_hash_scan_bitmask_candidate!(
                                out,
                                raw,
                                hit.start:(hit.start + prefix_scan_candidate_last_offset(
                                    geometry::PrefixScanGeometry)),
                                geometry::PrefixScanGeometry,
                                dbi,
                                is_antisense,
                                hit.mask,
                                guides,
                                guides_,
                                chrom_name,
                                distance,
                                early_stopping,
                                es_acc,
                                is_es,
                                seen,
                                stats;
                                distance_first = use_distance_first,
                                myers_profiles = myers_profiles,
                            )
                        end
                    end
                    continue
                end

                convert_start = time_ns()
                chrom_seq = dbi.gi.is_fa ?
                    FASTA.sequence(LongDNA{4}, record) :
                    TwoBit.sequence(LongDNA{4}, record)
                sequence_convert_ns = time_ns() - convert_start
                if stats !== nothing
                    stats.record_io_ns += record_io_ns
                    stats.sequence_convert_ns += sequence_convert_ns
                    stats.chrom_load_ns += record_io_ns + sequence_convert_ns
                end

                if use_fused_scan
                    plus_hits, minus_hits = scan_prefix_hits(
                        geometry::PrefixScanGeometry, chrom_seq, dbi, query,
                        hash_len, stats; scan_threads = scan_threads)
                    for (is_antisense, hits) in ((false, plus_hits), (true, minus_hits))
                        if stats !== nothing
                            stats.prefix_hits += length(hits)
                            stats.guide_pairs += sum(hit -> count_ones(hit.mask), hits; init = 0)
                        end
                        for hit in hits
                            candidate_range = hit.start:(hit.start +
                                prefix_scan_candidate_last_offset(
                                    geometry::PrefixScanGeometry))
                            verify_prefix_hash_scan_bitmask_candidate!(
                                out,
                                chrom_seq,
                                candidate_range,
                                geometry::PrefixScanGeometry,
                                dbi,
                                is_antisense,
                                hit.mask,
                                guides,
                                guides_,
                                chrom_name,
                                distance,
                                early_stopping,
                                es_acc,
                                is_es,
                                seen,
                                stats;
                                distance_first = use_distance_first,
                            )
                        end
                    end
                    continue
                end

                for is_antisense in (false, true)
                    findguides_start = time_ns()
                    positions = findguides(dbi, chrom_seq, is_antisense)
                    if stats !== nothing
                        stats.findguides_ns += time_ns() - findguides_start
                    end
                    isempty(positions) && continue
                    strand = is_antisense ? "-" : "+"

                    for candidate_range in positions
                        if stats !== nothing
                            stats.motif_candidates += 1
                        end

                        candidate_mask = zero(UInt64)
                        if use_bruteforce_query
                            guide_count = count(!, is_es)
                            has_candidate_guides = guide_count != 0
                        else
                            hash_start = time_ns()
                            if use_direct_specialized_hash
                                hashes = candidate_prefix_hashes_direct(
                                    geometry::PrefixScanGeometry, chrom_seq,
                                    candidate_range, is_antisense, hash_len, hash_type)
                            else
                                hashes = nothing
                            end
                            if hashes === nothing
                                if stats !== nothing
                                    stats.candidate_hash_ns += time_ns() - hash_start
                                end
                                prefix_start = time_ns()
                                prefix = normalized_candidate_prefix(chrom_seq, candidate_range, dbi, is_antisense, hash_len)
                                if stats !== nothing
                                    stats.candidate_prefix_ns += time_ns() - prefix_start
                                end
                                hash_start = time_ns()
                                hashes = candidate_prefix_hashes(prefix, hash_type, stats)
                            end
                            if stats !== nothing
                                stats.candidate_hash_ns += time_ns() - hash_start
                            end
                            lookup_start = time_ns()
                            if use_bitmask_query
                                candidate_mask = prefix_hash_scan_candidate_mask(query, hashes)
                                has_candidate_guides = candidate_mask != 0
                                guide_count = count_ones(candidate_mask)
                            else
                                has_candidate_guides = append_prefix_hash_scan_guides!(candidate_guides, query, hashes)
                                guide_count = length(candidate_guides)
                            end
                            if stats !== nothing
                                stats.query_lookup_ns += time_ns() - lookup_start
                            end
                        end
                        has_candidate_guides || continue

                        if stats !== nothing
                            stats.prefix_hits += 1
                            stats.guide_pairs += guide_count
                            if use_bruteforce_query
                                stats.bruteforce_guide_pairs += guide_count
                            end
                        end

                        materialize_start = time_ns()
                        ot, pos = materialize_normalized_candidate(chrom_seq, candidate_range, dbi, is_antisense)
                        if stats !== nothing
                            stats.candidate_materialize_ns += time_ns() - materialize_start
                        end
                        verify_start = time_ns()
                        if use_bitmask_query
                            mask = candidate_mask
                            while mask != 0
                                guide_idx = trailing_zeros(mask) + 1
                                mask &= mask - 1
                                emit_prefix_hash_scan_legacy_hit!(
                                    out, guide_idx, ot, pos, chrom_name,
                                    strand, guides, guides_, distance, motif,
                                    seen, es_acc, is_es, early_stopping, stats)
                            end
                        elseif use_bruteforce_query
                            for guide_idx in eachindex(guides_)
                                emit_prefix_hash_scan_legacy_hit!(
                                    out, guide_idx, ot, pos, chrom_name,
                                    strand, guides, guides_, distance, motif,
                                    seen, es_acc, is_es, early_stopping, stats)
                            end
                        else
                            for guide_idx in candidate_guides
                                emit_prefix_hash_scan_legacy_hit!(
                                    out, guide_idx, ot, pos, chrom_name,
                                    strand, guides, guides_, distance, motif,
                                    seen, es_acc, is_es, early_stopping, stats)
                            end
                        end
                        if stats !== nothing
                            stats.verify_ns += time_ns() - verify_start
                        end
                    end
                end
            end
            if stats !== nothing
                stats.scan_ns += time_ns() - scan_start
            end
        finally
            close(ref)
        end
    end
    merge_prefix_hash_scan_early_stop_stats!(stats, early_stop_state)
    return
end

function search_prefix_hash_scan_legacy_counts(
    dbi::DBInfo,
    reference_lengths,
    query,
    guides_::Vector{LongDNA{4}},
    geometry,
    distance::Int,
    hash_len::Int,
    hash_type::Type{<:Unsigned},
    use_bruteforce_query::Bool,
    use_bitmask_query::Bool,
    use_direct_specialized_hash::Bool,
    stats::Union{Nothing, PrefixHashScanStats})

    counts = zeros(Int, length(guides_), distance + 1)
    candidate_guides = Int[]
    ref = open(dbi.gi.filepath, "r")
    try
        reader = dbi.gi.is_fa ? FASTA.Reader(ref; copy = false) : TwoBit.Reader(ref)
        records = dbi.gi.is_fa ?
            PrefixHashScanFASTARecords(
                reader, dbi.gi.chrom, reference_lengths) :
            PrefixHashScanIndexedRecords(reader, dbi.gi.chrom)
        for (_, record, record_io_ns) in records
            convert_start = prefix_hash_scan_timer(stats)
            chrom_seq = dbi.gi.is_fa ?
                FASTA.sequence(LongDNA{4}, record) :
                TwoBit.sequence(LongDNA{4}, record)
            if stats !== nothing
                sequence_convert_ns = time_ns() - convert_start
                stats.record_io_ns += record_io_ns
                stats.sequence_convert_ns += sequence_convert_ns
                stats.chrom_load_ns += record_io_ns + sequence_convert_ns
            end
            for is_antisense in (false, true)
                findguides_start = prefix_hash_scan_timer(stats)
                positions = findguides(dbi, chrom_seq, is_antisense)
                stats === nothing ||
                    (stats.findguides_ns += time_ns() - findguides_start)
                for candidate_range in positions
                    stats === nothing || (stats.motif_candidates += 1)
                    candidate_mask = zero(UInt64)
                    if use_bruteforce_query
                        candidate_indices = eachindex(guides_)
                        guide_count = count_ones(
                            (length(guides_) == 64 ? typemax(UInt64) :
                                (UInt64(1) << length(guides_)) - UInt64(1)))
                    else
                        hash_start = prefix_hash_scan_timer(stats)
                        hashes = use_direct_specialized_hash ?
                            candidate_prefix_hashes_direct(
                                geometry::PrefixScanGeometry, chrom_seq,
                                candidate_range, is_antisense, hash_len,
                                hash_type) : nothing
                        if hashes === nothing
                            stats === nothing ||
                                (stats.candidate_hash_ns += time_ns() - hash_start)
                            prefix_start = prefix_hash_scan_timer(stats)
                            prefix = normalized_candidate_prefix(
                                chrom_seq, candidate_range, dbi,
                                is_antisense, hash_len)
                            stats === nothing ||
                                (stats.candidate_prefix_ns +=
                                    time_ns() - prefix_start)
                            hash_start = prefix_hash_scan_timer(stats)
                            hashes = candidate_prefix_hashes(
                                prefix, hash_type, stats)
                        end
                        stats === nothing ||
                            (stats.candidate_hash_ns += time_ns() - hash_start)
                        lookup_start = prefix_hash_scan_timer(stats)
                        if use_bitmask_query
                            candidate_mask = prefix_hash_scan_candidate_mask(
                                query, hashes)
                            has_candidate_guides = candidate_mask != 0
                            candidate_indices = nothing
                            guide_count = count_ones(candidate_mask)
                        else
                            has_candidate_guides =
                                append_prefix_hash_scan_guides!(
                                    candidate_guides, query, hashes)
                            candidate_indices = candidate_guides
                            guide_count = length(candidate_guides)
                        end
                        stats === nothing ||
                            (stats.query_lookup_ns += time_ns() - lookup_start)
                        has_candidate_guides || continue
                    end
                    if stats !== nothing
                        stats.prefix_hits += 1
                        stats.guide_pairs += guide_count
                        use_bruteforce_query &&
                            (stats.bruteforce_guide_pairs += guide_count)
                    end

                    materialize_start = prefix_hash_scan_timer(stats)
                    ot, _ = materialize_normalized_candidate(
                        chrom_seq, candidate_range, dbi, is_antisense)
                    stats === nothing ||
                        (stats.candidate_materialize_ns +=
                            time_ns() - materialize_start)
                    verify_start = prefix_hash_scan_timer(stats)
                    if use_bitmask_query
                        while candidate_mask != 0
                            guide_idx = trailing_zeros(candidate_mask) + 1
                            candidate_mask &= candidate_mask - 1
                            align_start = prefix_hash_scan_timer(stats)
                            if stats !== nothing
                                stats.alignment_calls += 1
                                stats.distance_calls += 1
                            end
                            dist = levenshtein(
                                guides_[guide_idx], ot, distance, iscompatible)
                            stats === nothing ||
                                (stats.align_ns += time_ns() - align_start)
                            if dist <= distance
                                counts[guide_idx, dist + 1] += 1
                            end
                        end
                    else
                        for guide_idx in candidate_indices
                            align_start = prefix_hash_scan_timer(stats)
                            if stats !== nothing
                                stats.alignment_calls += 1
                                stats.distance_calls += 1
                            end
                            dist = levenshtein(
                                guides_[guide_idx], ot, distance, iscompatible)
                            stats === nothing ||
                                (stats.align_ns += time_ns() - align_start)
                            if dist <= distance
                                counts[guide_idx, dist + 1] += 1
                            end
                        end
                    end
                    stats === nothing ||
                        (stats.verify_ns += time_ns() - verify_start)
                end
            end
        end
    finally
        close(ref)
    end
    return counts
end

function write_prefix_hash_scan_counts(
    output_file::String,
    guides::Vector{LongDNA{4}},
    counts::Matrix{Int},
    early_stopping::Vector{Int},
    stats::Union{Nothing, PrefixHashScanStats} = nothing)

    guide_rows = String[]
    row_by_guide = Dict{String, Int}()
    grouped_counts = Vector{Vector{Int}}()
    complete = Bool[]
    for guide_idx in eachindex(guides)
        guide = string(guides[guide_idx])
        row_idx = get(row_by_guide, guide, 0)
        if row_idx == 0
            push!(guide_rows, guide)
            push!(grouped_counts, zeros(Int, size(counts, 2)))
            push!(complete, true)
            row_idx = length(guide_rows)
            row_by_guide[guide] = row_idx
        end
        for dist_idx in axes(counts, 2)
            raw_count = counts[guide_idx, dist_idx]
            grouped_counts[row_idx][dist_idx] +=
                min(raw_count, early_stopping[dist_idx])
            raw_count > early_stopping[dist_idx] &&
                (complete[row_idx] = false)
        end
    end

    mkpath(dirname(output_file))
    emit_start = prefix_hash_scan_timer(stats)
    open(output_file, "w") do out
        print(out, "guide")
        for distance in 0:(size(counts, 2) - 1)
            print(out, ",D", distance)
        end
        print(out, ",complete\n")
        for row_idx in eachindex(guide_rows)
            print(out, guide_rows[row_idx])
            for count in grouped_counts[row_idx]
                print(out, ",", count)
            end
            print(out, ",", complete[row_idx], "\n")
        end
    end
    if stats !== nothing
        stats.emit_ns += time_ns() - emit_start
        stats.emitted_rows += length(guide_rows)
    end
    return nothing
end

function with_atomic_prefix_hash_scan_output(
    f::Function, output_file::String)

    output_dir = dirname(output_file)
    mkpath(output_dir)
    temporary, temporary_io = mktemp(output_dir; cleanup = false)
    close(temporary_io)
    try
        f(temporary)
        Base.Filesystem.rename(temporary, output_file; force = true)
    finally
        ispath(temporary) && rm(temporary; force = true)
    end
    return nothing
end

function prefix_hash_scan_pamless_motif(
    guide_bases::Int, distance::Int, ambig_max::Int = 0)

    return Motif(
        "PAMless_$(guide_bases)",
        repeat("N", guide_bases),
        repeat("X", guide_bases),
        true, true, distance, true, ambig_max)
end

"""
```
search_prefixHashScan(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    output_file::String;
    motif::Union{Nothing, String, Motif} = nothing,
    pamless::Bool = false,
    ambig_max::Union{Nothing, Int} = nothing,
    distance::Int = 3,
    early_stopping::Vector{Int} = fill(1_000_000, distance + 1),
    scan_threads::Int = Threads.nthreads(),
    simd_backend::Symbol = :auto,
    output::Symbol = :detail,
    verbose::Bool = false)
```

Find all off-targets for `guides` within `distance` by scanning `genome_path`
directly. Unlike prefixHashDB this builds no CHOPOFF database: it needs only an
indexed FASTA (with the standard `.fai`) or a `.2bit` reference. Internally it
reuses prefixHashDB's symbolic prefix paths and scans the reference with SIMD
kernels.

Assumes your guides do not contain PAM, and are all in the same direction as
you would order for the lab e.g.:

```
5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM

3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide
```

# Arguments

`genome_path` - Path to an indexed FASTA or a `.2bit` reference.

`output_file` - Path and name for the output file, this will be comma
separated table, therefore `.csv` extension is preferred. Guide lists larger
than 64 are searched as sequential 64-guide batches; the requested path is
replaced only after every batch completes successfully.

`motif` - Registered motif name, a custom `Motif`, or `nothing` for Cas9. The
motif's distance is set from `distance`, so it need not be preset.

`pamless` - Set `true` to infer a homogeneous, both-strand PAMless motif from
the guide length. Every guide in one call must then be the same length. Cannot
be combined with `motif`.

`ambig_max` - Overrides the inferred or supplied motif's `ambig_max`, that is,
how many ambiguous IUPAC reference bases a candidate guide/PAM window may
contain. Supported range is 0 through 3. Query guides must be unambiguous.

`distance` - Defines maximum levenshtein distance (insertions, deletions,
mismatches) for which off-targets are considered. Supported range is 0
through 4.

`early_stopping` - Integer vector. Early stopping condition. For example for
distance 2, we need vector with 3 values e.g. [1, 1, 5]. Which means we will
search with "up to 1 offtargets within distance 0", "up to 1 offtargets within
distance 1"...

`scan_threads` - Threads used for the genome scan.

`simd_backend` - `:auto`, `:avx2`, or `:avx512` for raw SIMD scans. `:auto`
resolves from the CPU and the scan geometry.

`output` - `:detail` writes aligned loci with the usual CHOPOFF columns.
`:counts` writes one row per unique guide as `guide,D0,...,Dk,complete`,
without traceback; counts above an early-stopping limit are capped and the row
is marked incomplete. This is not the same schema as `summarize_offtargets`,
which has no `complete` column.

See the four-argument method for the scan tuning options.

# Examples
```julia
$(make_prefix_hash_scan_example_doc())
```
"""
function search_prefixHashScan(
    guides::Vector{LongDNA{4}},
    genome_path::String,
    output_file::String;
    motif::Union{Nothing, String, Motif} = nothing,
    pamless::Bool = false,
    ambig_max::Union{Nothing, Int} = nothing,
    distance::Int = 3,
    early_stopping::Vector{Int} = fill(1_000_000, distance + 1),
    scan_threads::Int = Threads.nthreads(),
    simd_backend::Symbol = :auto,
    output::Symbol = :detail,
    verbose::Bool = false)

    distance in 0:4 ||
        error("search_prefixHashScan supports distances 0 through 4.")
    isempty(guides) &&
        error("search_prefixHashScan requires at least one guide.")
    guide_bases = length(first(guides))
    if pamless
        motif === nothing ||
            error("pamless=true cannot be combined with motif.")
        all(==(guide_bases), length.(guides)) ||
            error("search_prefixHashScan requires homogeneous PAMless guide lengths.")
        motif_ = prefix_hash_scan_pamless_motif(
            guide_bases, distance, something(ambig_max, 0))
    else
        motif_ = motif === nothing ? Motif("Cas9"; distance = distance) :
            (motif isa String ? Motif(motif; distance = distance) :
                setdist(motif, distance))
        ambig_max === nothing || (motif_ = setambig(motif_, ambig_max))
    end
    guide_bases = length_noPAM(motif_)
    all(==(guide_bases), length.(guides)) ||
        error("search_prefixHashScan requires $guide_bases-base guides for $(motif_.alias).")
    validate_prefix_hash_scan_query(
        guides, motif_, distance, early_stopping, output)
    hash_len = min(guide_bases - distance, 16)
    hash_len >= 1 ||
        error("search_prefixHashScan requires guide length greater than distance.")
    prepared_paths, _ = load_prefix_hash_scan_paths(
        motif_, distance, hash_len)
    batch_size = 64
    batch_count = cld(length(guides), batch_size)
    if verbose && batch_count > 1
        @info(
            "prefixHashScan guide batching",
            guides = length(guides),
            batch_size = batch_size,
            batches = batch_count,
        )
    end
    counts = output == :counts ? zeros(Int, length(guides), distance + 1) : nothing
    function run_batches(batch_output_file::String)
        for (batch_idx, first_idx) in enumerate(1:batch_size:length(guides))
            last_idx = min(first_idx + batch_size - 1, length(guides))
            result = search_prefixHashScan(
                guides[first_idx:last_idx],
                genome_path,
                motif_,
                batch_output_file;
                distance = distance,
                hash_len = hash_len,
                early_stopping = early_stopping,
                query_variant = :bitmask64,
                scan_backend = :auto,
                bucket_bases = 11,
                scan_threads = scan_threads,
                simd_backend = simd_backend,
                stream_chunk_bases = 2 * 1024 * 1024,
                prefilter_bits = 26,
                query_build_backend = :auto,
                lookup_variant = :auto,
                verify_variant = :auto,
                output = output,
                verbose = verbose && batch_idx == 1,
                _append_output = batch_idx > 1,
                _collect_counts = output == :counts,
                _paths = prepared_paths,
            )
            output == :counts && (counts[first_idx:last_idx, :] .= result)
        end
        return nothing
    end
    if output == :detail && batch_count > 1
        with_atomic_prefix_hash_scan_output(output_file) do temporary
            run_batches(temporary)
        end
    else
        run_batches(output_file)
    end
    output == :counts && write_prefix_hash_scan_counts(
        output_file, guides, counts, early_stopping)
    return nothing
end
