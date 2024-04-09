using CHOPOFF

motif = Motif("Cas9"; distance = 4)
mpt = build_PathTemplates(motif; restrict_to_len = 16)
paths = UInt8.(mpt.paths)
distances = UInt8.(mpt.distances)

CHOPOFF.save(distances, "./data/Cas9_d4_p16_distances.bin")
CHOPOFF.save(paths, "./data/Cas9_d4_p16_paths.bin")

d2 = CHOPOFF.load("./data/Cas9_d4_p16_distances.bin")
p2 = CHOPOFF.load("./data/Cas9_d4_p16_paths.bin")
if (paths != p2) | (distances != d2)
    @warn "Failed to sucessfully save the path templates. CHOPOFF will still work, but will be slower in some cases."
    rm("./data/Cas9_d4_p16_paths.bin")
    rm("./data/Cas9_d4_p16_distances.bin")
end