struct PAMinFMI
    motif::Motif
    genome::String
    pam_loc_fwd::IdDict{UInt8, Vector{Int}} # sorted
    pam_loc_rve::IdDict{UInt8, Vector{Int}}
end

