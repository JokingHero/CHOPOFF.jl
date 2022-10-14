"
Contians information about which genome was used.

date - date and time of genome access
filepath - path used to access genome
genomechecksum - checksum of that file
chrom - vector of chromosome names
pos_type - Smallest possible DataType that can be used to mark positions on the genome.
"
struct GenomeInfo
    date::DateTime
    filepath::String
    genomechecksum::UInt32
    chrom::Vector{String}
    chrom_type::DataType
    pos_type::DataType
    is_fa::Bool
end


function is_fasta(filepath::String)
    ext = extension(filepath)
    is_fa = false
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        error("Wrong extension of the genome.",
              "We can parse only .fa/.fna/.fasta or .2bit references.")
    end
    return is_fa
end


"Instantiate GenomeInfo object from genom file path."
function GenomeInfo(filepath::String)
    checksum = open(crc32c, filepath)

    is_fa = is_fasta(filepath)
    ref = open(filepath, "r")
    maxchromlen = 0

    if !is_fa
        reader = TwoBit.Reader(ref)
        chrom = TwoBit.seqnames(reader)
        for record in reader
            if TwoBit.hassequence(record)
                tl = record.dnasize
                maxchromlen = maxchromlen >  tl ? maxchromlen : tl
            end
        end
    else
        chrom = Vector{String}()
        for record in FASTA.Reader(ref)
            if FASTA.hasidentifier(record) & FASTA.hassequence(record)
                push!(chrom, FASTA.identifier(record))
                tl = FASTA.seqlen(record)
                maxchromlen = maxchromlen > tl ? maxchromlen : tl
            end
        end
    end
    close(ref)
    pos_type = smallestutype(unsigned(maxchromlen))
    chrom_type = smallestutype(unsigned(length(chrom)))
    return GenomeInfo(now(Dates.UTC), filepath, checksum, chrom, chrom_type, pos_type, is_fa)
end

function Base.isequal(gi::GenomeInfo, gi2::GenomeInfo)
    return isequal(gi.genomechecksum, gi2.genomechecksum) & isequal(Set(gi.chrom), Set(gi2.chrom))
end


struct DBInfo
    name::String
    date::DateTime
    gi::GenomeInfo
    vcf_filepath::String
    motif::Motif
end


"""
`DBInfo(filepath::String, name::String, motif::Motif; vcf_filepath::String = "")`


Motif defines what genome file is being used for searches.


# Arguments
`filepath` - Path to the genome file, if file is fasta (ends with .fa or .fasta or .fna) 
make sure you also have fasta index file with extension .fai. Alternatively, you can use .2bit 
genome file.  

`name` - Your name for this instance of DBInfo: the genome with connection to the motif and vcf file.

`motif`   - `Motif` object defining search parameters

`vcf_filepath`  - Optional. Path to the VCF file to include in the searches.


Alignments will be performed from opposite to the extension direction (which is deifned by extend5).

# Examples
```julia-repl
# use ARTEMIS example genome
genome = joinpath(
    vcat(
        splitpath(dirname(pathof(ARTEMIS)))[1:end-1], 
        "test", "sample_data", "genome", "semirandom.fa"))
# construct example DBInfo
DBInfo(genome, "Cas9_semirandom_noVCF", Motif("Cas9"))
```
"""
function DBInfo(filepath::String, name::String, motif::Motif; vcf_filepath::String = "")
    gi = GenomeInfo(filepath::String)
    return DBInfo(name, now(Dates.UTC), gi, vcf_filepath, motif)
end


struct Loc{T<:Unsigned,K<:Unsigned}
    chrom::T
    pos::K
    isplus::Bool
end


function decode(loc::Loc, dbi::DBInfo)
    strand = loc.isplus ? "+" : "-"
    return dbi.gi.chrom[loc.chrom] * "," * string(loc.pos) * "," * strand
end