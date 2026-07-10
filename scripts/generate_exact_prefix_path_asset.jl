#!/usr/bin/env julia

using CHOPOFF
using DataFrames
using Printf

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

function generate_exact_asset(asset::String, distance::Int, hash_len::Int)
    pfile1 = joinpath(DATA_DIR, "$(asset)_d4_p16_paths_part1.bin")
    pfile2 = joinpath(DATA_DIR, "$(asset)_d4_p16_paths_part2.bin")
    dfile = joinpath(DATA_DIR, "$(asset)_d4_p16_distances.bin")
    isfile(pfile1) || error("Missing $pfile1")
    isfile(pfile2) || error("Missing $pfile2")
    isfile(dfile) || error("Missing $dfile")

    @info "Loading full asset" asset distance hash_len
    paths = vcat(CHOPOFF.load(pfile1), CHOPOFF.load(pfile2))
    distances = CHOPOFF.load(dfile)
    paths = paths[:, 1:hash_len]

    @info "Filtering full asset"
    keep = BitVector(distances .<= distance)
    paths = paths[keep, :]
    distances = distances[keep]
    paths = convert.(CHOPOFF.smallestutype(maximum(paths)), paths)

    stem = joinpath(DATA_DIR, "$(asset)_d$(distance)_p$(hash_len)")
    CHOPOFF.save(paths, stem * "_paths.bin")
    CHOPOFF.save(distances, stem * "_distances.bin")

    paths2 = CHOPOFF.load(stem * "_paths.bin")
    distances2 = CHOPOFF.load(stem * "_distances.bin")
    paths == paths2 || error("Saved paths differ after reload")
    distances == distances2 || error("Saved distances differ after reload")
    @printf("saved %s rows=%d cols=%d\n", stem, size(paths, 1), size(paths, 2))
end

asset = get(ENV, "CHOPOFF_EXACT_PATH_ASSET", "Cas9")
distance = parse_int_env("CHOPOFF_EXACT_PATH_DISTANCE", 3)
hash_len = parse_int_env("CHOPOFF_EXACT_PATH_HASH_LEN", 16)
generate_exact_asset(asset, distance, hash_len)
