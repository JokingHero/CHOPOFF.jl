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

"
Contains basic information about the
genome for which databases were indexed.

vcf_filepath - path to the file with snp annotations, or ''
"
struct DBInfo
    name::String
    date::DateTime
    gi::GenomeInfo
    vcf_filepath::String
    motif::Motif
end

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
