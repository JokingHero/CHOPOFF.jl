function buildvpDB(name::String,
    genomepath::String,
    motif::Motif,
    prefix::Int)

    gi = GenomeInfo(genomepath, name)
    # "kmer" => "guide" => [Loc, Loc, Loc]
    #db = IdDict{String, IdDict{LongDNASeq, Vector{Loc}}}()
    #gatherofftargets!(gi, motif, db)
    return #db
end
