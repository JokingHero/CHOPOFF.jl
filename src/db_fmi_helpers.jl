struct PAMinFMI
    motif::Motif
    gi::GenomeInfo
    pam_loc_fwd::IdDict{Int, Vector{Int}} # sorted
    pam_loc_rve::IdDict{Int, Vector{Int}}
end 

# this stores location of the PAMs
# its important that PAMs positions here
# are NGG  CCN TTTN  NAAA
#     ^      ^    ^  ^
function build_pamDB(fmidbdir::String, motif::Motif; storage_dir::String = "")
    
    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    pam = motif.fwd[motif.pam_loci_fwd]

    pam = expand_ambiguous(pam)
    pam_loc_fwd = IdDict{Int, Vector{Int}}()
    pam_loc_rve = IdDict{Int, Vector{Int}}()

    for (ic, chrom) in enumerate(gi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        pam_loc_fwd_chrom = Vector{gi.pos_type}()
        pam_loc_rve_chrom = Vector{gi.pos_type}()
        for pam_i in pam
            pam_pos = ARTEMIS.locateall(pam_i, fmi) # NGG on N
            if !motif.extends5 # TTTN on N
                pam_pos = pam_pos .+ (length(pam_i) - 1)
            end
            append!(pam_loc_fwd_chrom, pam_pos)
            
            pam_pos = ARTEMIS.locateall(reverse_complement(pam_i), fmi) # NAAA on N
            if motif.extends5 # CCN on N
                pam_pos = pam_pos .+ (length(pam_i) - 1)
            end
            append!(pam_loc_rve_chrom, pam_pos)
        end

        sort!(pam_loc_fwd_chrom)
        sort!(pam_loc_rve_chrom)

        pam_loc_fwd[ic] = pam_loc_fwd_chrom
        pam_loc_rve[ic] = pam_loc_rve_chrom
    end

    pamDB = PAMinFMI(motif, gi, pam_loc_fwd, pam_loc_rve)

    if storage_dir != ""
        save(pamDB, storage_dir)
    end
    return pamDB
end


function build_fmiDB(
    genomepath::String,
    storage_dir::String)

    gi = GenomeInfo(genomepath)
    ref = open(gi.filepath, "r")
    reader = gi.is_fa ? FASTA.Reader(ref, index = gi.filepath * ".fai") : TwoBit.Reader(ref)
    @showprogress 60 for chrom in gi.chrom # no need for paralelization as this is super fast
        fmi = FMIndex(getchromseq(gi.is_fa, reader[chrom]), 16, r = 32)
        p = joinpath(storage_dir, chrom * ".bin")
        save(fmi, p)
    end
    close(ref)
    save(gi, joinpath(storage_dir, "genomeInfo.bin"))
    @info "FMindex is build."
    return storage_dir
end