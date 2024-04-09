using CHOPOFF
using DataFrames

motif = Motif("Cas9"; distance = 4)
mpt = build_PathTemplates(motif; restrict_to_len = 16)

paths = mpt.paths[:, 1:16]
not_dups = map(!, BitVector(nonunique(DataFrame(paths, :auto))))
paths = paths[not_dups, :]
distances = mpt.distances[not_dups]

split = Int(floor(length(distances) / 2))
paths1 = paths[1:split, :]
paths2 = paths[(split + 1):end, :]

CHOPOFF.save(distances, "./data/Cas9_d4_p16_distances.bin")
CHOPOFF.save(paths1, "./data/Cas9_d4_p16_paths_part1.bin")
CHOPOFF.save(paths2, "./data/Cas9_d4_p16_paths_part2.bin")

d2 = CHOPOFF.load("./data/Cas9_d4_p16_distances.bin")
p1 = CHOPOFF.load("./data/Cas9_d4_p16_paths_part1.bin")
p2 = CHOPOFF.load("./data/Cas9_d4_p16_paths_part2.bin")
p2 = vcat(p1, p2)
if (paths != p2) | (distances != d2)
    @warn "Failed to sucessfully save the path templates. CHOPOFF will still work, but will be slower in some cases."
    rm("./data/Cas9_d4_p16_paths.bin")
    rm("./data/Cas9_d4_p16_distances.bin")
end