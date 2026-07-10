using CHOPOFF
using DataFrames

function compute_default_path_template(name::String; distance::Int = 4, hash_len::Int = 16)
    motif = Motif(name; distance = distance)
    mpt = build_PathTemplates(motif; restrict_to_len = hash_len, withPAM = false)

    paths = mpt.paths[:, 1:hash_len]
    not_dups = map(!, BitVector(nonunique(DataFrame(paths, :auto))))
    paths = paths[not_dups, :]
    distances = mpt.distances[not_dups]

    split = Int(floor(length(distances) / 2))
    paths1 = paths[1:split, :]
    paths2 = paths[(split + 1):end, :]
    stem = "./data/$(name)_d$(distance)_p$(hash_len)"

    CHOPOFF.save(distances, stem * "_distances.bin")
    CHOPOFF.save(paths1, stem * "_paths_part1.bin")
    CHOPOFF.save(paths2, stem * "_paths_part2.bin")

    d2 = CHOPOFF.load(stem * "_distances.bin")
    p1 = CHOPOFF.load(stem * "_paths_part1.bin")
    p2 = CHOPOFF.load(stem * "_paths_part2.bin")
    if (paths != vcat(p1, p2)) | (distances != d2)
        @warn "Failed to successfully save $name path templates. CHOPOFF will still work, but will be slower in some cases."
    end
    return size(paths, 1)
end

for name in ("Cas9", "Cas12a")
    n = compute_default_path_template(name)
    @info "Saved precomputed path templates" name rows=n
end
