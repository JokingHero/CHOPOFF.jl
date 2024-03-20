struct PAMinFMI
    motif::Motif
    gi::GenomeInfo
    pam_loc_fwd::IdDict{Int, Vector{Int}} # sorted
    pam_loc_rve::IdDict{Int, Vector{Int}}
end


"""
```
build_pamDB(
    fmidbdir::String, 
    motif::Motif; 
    storage_path::String = "")
```

Find locations of all the PAM for a given motif in the genome.

Example of what position we store for traditional Cas9 (NGG, CCN) and Cas12a (TTTV, BAAA):
```
NGG  CCN TTTN  NAAA
^      ^    ^  ^
```

# Arguments

`fmidbdir` - Path to directory with the FM-index build with `build_fmiDB`.

`motif` - `Motif` defines which PAM we will locate in the genome.

`storage_path`  - Path to the DIRECTORY where index will be saved.


# Examples
```julia-repl
# prepare libs
using CHOPOFF, BioSequences

# make a temporary directory
tdir = tempname()
fmi_dir = joinpath(tdir, "fmi")
mkpath(fmi_dir)

# use CHOPOFF example genome
genome = joinpath(
    vcat(
        splitpath(dirname(pathof(CHOPOFF)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))
# build FM-index
build_fmiDB(genome, fmi_dir)

# build a pamDB
pamDB = build_pamDB(fmi_dir, Motif("Cas9"))
```
"""
function build_pamDB(fmidbdir::String, motif::Motif; storage_path::String = "")
    
    gi = load(joinpath(fmidbdir, "genomeInfo.bin"))
    pam = motif.fwd[motif.pam_loci_fwd]

    pam, idx = expand_ambiguous(pam)
    pam_loc_fwd = IdDict{Int, Vector{Int}}()
    pam_loc_rve = IdDict{Int, Vector{Int}}()

    for (ic, chrom) in enumerate(gi.chrom)
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        pam_loc_fwd_chrom = Vector{gi.pos_type}()
        pam_loc_rve_chrom = Vector{gi.pos_type}()
        for pam_i in pam
            pam_pos = CHOPOFF.locateall(pam_i, fmi) # NGG on N
            if !motif.extends5 # TTTN on N
                pam_pos = pam_pos .+ (length(pam_i) - 1)
            end
            append!(pam_loc_fwd_chrom, pam_pos)
            
            pam_pos = CHOPOFF.locateall(reverse_complement(pam_i), fmi) # NAAA on N
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

    if storage_path != ""
        save(pamDB, storage_path)
    end
    return pamDB
end


"""
```
build_fmiDB(
    genomepath::String,
    storage_dir::String)
```

Prepare FM-index for future searches.


# Arguments

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`storage_dir`  - Path to the DIRECTORY where index with many files will be saved.


# Examples
```julia-repl
# prepare libs
using CHOPOFF, BioSequences

# make a temporary directory
tdir = tempname()
fmi_dir = joinpath(tdir, "fmi")
mkpath(fmi_dir)

# use CHOPOFF example genome
genome = joinpath(vcat(splitpath(dirname(pathof(CHOPOFF)))[1:end-1], 
    "test", "sample_data", "genome", "semirandom.fa"))

# build a fmiDB!
build_fmiDB(genome, fmi_dir)
```
"""
function build_fmiDB(
    genomepath::String,
    storage_dir::String)

    gi = GenomeInfo(genomepath)
    ref = open(gi.filepath, "r")
    reader = gi.is_fa ? FASTA.Reader(ref, index = gi.filepath * ".fai") : TwoBit.Reader(ref)
    @showprogress dt=60 for chrom in gi.chrom # no need for parallelization as this is super fast
        fmi = FMIndex(getchromseq(gi.is_fa, reader[chrom]), 16, r = 32)
        p = joinpath(storage_dir, chrom * ".bin")
        save(fmi, p)
    end
    close(ref)
    save(gi, joinpath(storage_dir, "genomeInfo.bin"))
    @info "FMindex is build."
    return storage_dir
end