struct PAMinFMI
    motif::Motif
    gi::GenomeInfo
    pam_loc_fwd::IdDict{Int, Vector{Int}} # sorted
    pam_loc_rve::IdDict{Int, Vector{Int}}
end 


function build_pamDB(fmidbdir::String, motif::Motif; storagedir::String = "")
    
    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    pam = motif.fwd[motif.pam_loci_fwd]
    if motif.extends5
        pam = reverse(pam)
    end

    pam = expand_ambiguous(pam)
    pam_loc_fwd = IdDict{Int, Vector{Int}}()
    pam_loc_rve = IdDict{Int, Vector{Int}}()

    for (ic, chrom) in enumerate(gi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        pam_loc_fwd_chrom = Vector{gi.pos_type}()
        pam_loc_rve_chrom = Vector{gi.pos_type}()
        for pam_i in pam
            # fwd
            append!(pam_loc_fwd_chrom, ARTEMIS.locateall(reverse(pam_i), fmi))
            # rve
            append!(pam_loc_rve_chrom, ARTEMIS.locateall(reverse_complement(pam_i), fmi))
        end

        sort!(pam_loc_fwd_chrom)
        sort!(pam_loc_rve_chrom)

        pam_loc_rve_chrom = pam_loc_rve_chrom .+ length(pam) .- 1
        pam_loc_fwd[ic] = pam_loc_fwd_chrom
        pam_loc_rve[ic] = pam_loc_rve_chrom
    end

    pamDB = PAMinFMI(motif, gi, pam_loc_fwd, pam_loc_rve)

    if storagedir != ""
        save(pamDB, storagedir)
    end
    return pamDB
end


function build_fmiDB(
    genomepath::String,
    storagedir::String)

    gi = GenomeInfo(genomepath)
    ref = open(gi.filepath, "r")
    reader = gi.is_fa ? FASTA.Reader(ref, index = gi.filepath * ".fai") : TwoBit.Reader(ref)
    @showprogress 60 for chrom in gi.chrom
        fmi = FMIndex(getchromseq(gi.is_fa, reader[chrom]), 16, r = 32)
        p = joinpath(storagedir, chrom * ".bin")
        save(fmi, p)
    end
    close(ref)
    save(gi, joinpath(storagedir, "genomeInfo.bin"))
    @info "FMindex is build."
    return storagedir
end


"""
Filters out super close and neighbouring off-targets, 
but only within small distance.
"""
function filter_overlaps(pos)
    pos = Set.(pos)
    pos_filtered = copy(pos)
    for i in 2:length(pos)
        pos[i - 1] = union(pos[i - 1], pos[i - 1] .+ 1, pos[i - 1] .- 1)
        pos_filtered[i] = setdiff(pos[i], pos[i - 1])
    end
    return pos_filtered
end