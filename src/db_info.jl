"
Contains basic information about the
genome for which databases were indexed.

date - date and time of genome access
filepath - path used to access genome
genomechecksum - checksum of that file
genome - user name of the genome
chrom - vector of chromosome names
vcf_filepath - path to the file with snp annotations, or ''
"
struct DBInfo
    name::String
    date::DateTime
    filepath::String
    genomechecksum::UInt32
    vcf_filepath::String
    chrom::Vector{String}
    chrom_type::DataType
    pos_type::DataType
    is_fa::Bool
    motif::Motif
end

function DBInfo(filepath::String, name::String, motif::Motif; vcf_filepath::String = "")
    checksum = open(crc32c, filepath)

    ext = extension(filepath)
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        error("Wrong extension of the genome.",
              "We can parse only .fa/.fna/.fasta or .2bit references.")
    end

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
    return DBInfo(name, now(Dates.UTC), filepath, checksum,
        vcf_filepath, chrom, chrom_type, pos_type, is_fa, motif)
end

struct Loc{T<:Unsigned,K<:Unsigned}
    chrom::T
    pos::K
    isplus::Bool
end

function decode(loc::Loc, dbi::DBInfo)
    strand = loc.isplus ? "+" : "-"
    return dbi.chrom[loc.chrom] * "," * string(loc.pos) * "," * strand
end
