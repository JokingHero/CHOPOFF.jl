#!/usr/bin/env julia

using CHOPOFF
using DataFrames
using Dates
using Printf
using Statistics

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = joinpath(ROOT_DIR, "data")

function parse_int_env(name::String, default::Int)
    raw = strip(get(ENV, name, string(default)))
    try
        return parse(Int, raw)
    catch
        error("Invalid integer in environment variable $name: '$raw'")
    end
end

function summarize(values::Vector{Float64})
    return (median = median(values), mean = mean(values), min = minimum(values), max = maximum(values))
end

function time_once(f)
    t = time_ns()
    value = f()
    return value, Float64(time_ns() - t) / 1e9
end

function profile_full_asset(asset::String, distance::Int, hash_len::Int)
    pfile1 = joinpath(DATA_DIR, "$(asset)_d4_p16_paths_part1.bin")
    pfile2 = joinpath(DATA_DIR, "$(asset)_d4_p16_paths_part2.bin")
    dfile = joinpath(DATA_DIR, "$(asset)_d4_p16_distances.bin")

    paths1, load_p1 = time_once(() -> CHOPOFF.load(pfile1))
    paths2, load_p2 = time_once(() -> CHOPOFF.load(pfile2))
    distances, load_dist = time_once(() -> CHOPOFF.load(dfile))
    paths, vcat_s = time_once(() -> vcat(paths1, paths2))
    paths, slice_s = time_once(() -> paths[:, 1:hash_len])
    not_dups, dedup_s = time_once(() -> map(!, BitVector(nonunique(DataFrame(paths, :auto)))))
    keep, mask_s = time_once(() -> not_dups .& BitVector(distances .<= distance))
    paths, filter_s = time_once(() -> paths[keep, :])
    paths, convert_s = time_once(() -> convert.(CHOPOFF.smallestutype(maximum(paths)), paths))
    return (
        rows = size(paths, 1),
        load_p1_s = load_p1,
        load_p2_s = load_p2,
        load_dist_s = load_dist,
        vcat_s = vcat_s,
        slice_s = slice_s,
        dedup_s = dedup_s,
        mask_s = mask_s,
        filter_s = filter_s,
        convert_s = convert_s,
    )
end

function profile_exact_asset(asset::String, distance::Int, hash_len::Int; need_distances::Bool)
    motif = Motif(asset; distance = distance)
    stats = Float64[]
    rows = 0
    for _ in 1:parse_int_env("CHOPOFF_PATH_PROFILE_RUNS", 5)
        GC.gc()
        t = @elapsed begin
            loaded = CHOPOFF.load_precomputed_prefix_paths(motif, distance, hash_len; need_distances = need_distances)
            loaded === nothing && error("Missing exact/fallback asset for $asset d$distance p$hash_len")
            paths, _, _ = loaded
            rows = size(paths, 1)
        end
        push!(stats, t)
    end
    return merge((rows = rows, need_distances = need_distances), summarize(stats))
end

function main()
    asset = get(ENV, "CHOPOFF_PATH_PROFILE_ASSET", "Cas9")
    distance = parse_int_env("CHOPOFF_PATH_PROFILE_DISTANCE", 3)
    hash_len = parse_int_env("CHOPOFF_PATH_PROFILE_HASH_LEN", 16)

    println("prefix path loading profile")
    println("time_utc: ", Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println("asset: $asset distance: $distance hash_len: $hash_len")

    full = profile_full_asset(asset, distance, hash_len)
    @printf("full d4/p16 fallback | rows=%d load_p1=%.3fs load_p2=%.3fs load_dist=%.3fs vcat=%.3fs slice=%.3fs dedup=%.3fs mask=%.3fs filter=%.3fs convert=%.3fs\n",
        full.rows, full.load_p1_s, full.load_p2_s, full.load_dist_s, full.vcat_s, full.slice_s, full.dedup_s, full.mask_s, full.filter_s, full.convert_s)

    exact_paths = joinpath(DATA_DIR, "$(asset)_d$(distance)_p$(hash_len)_paths.bin")
    if isfile(exact_paths)
        scan = profile_exact_asset(asset, distance, hash_len; need_distances = false)
        db = profile_exact_asset(asset, distance, hash_len; need_distances = true)
        @printf("exact paths-only | rows=%d median=%.3fs mean=%.3fs min=%.3fs max=%.3fs\n",
            scan.rows, scan.median, scan.mean, scan.min, scan.max)
        @printf("exact with distances | rows=%d median=%.3fs mean=%.3fs min=%.3fs max=%.3fs\n",
            db.rows, db.median, db.mean, db.min, db.max)
    else
        println("exact asset missing: $exact_paths")
    end
end

main()
