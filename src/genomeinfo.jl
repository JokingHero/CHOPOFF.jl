"
Contains basic information about the
genome for which databases were indexed.

date - date and time of genome access
filepath - path used to access genome
genomechecksum - checksum of that file
genome - user name of the genome
chrom - vector of chromosome names
snpfilepath - path to the file with snp annotations, or ''
"
struct GenomeInfo
    date::DateTime
    filepath::String
    genomechecksum::UInt32
    genome::String
    snpfilepath::String
    chrom::Vector{String}
    maxchromlen::BigInt
    is_fa::Bool
end

function GenomeInfo(filepath::String, name::String)
    checksum = open(crc32c, filepath)

    ext = extension(filepath)
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        throw(
        "Wrong extension of the genome.",
        "We can parse only .fa/.fna/.fasta or .2bit references.")
    end

    ref = open(filepath, "r")
    maxchromlen = 0

    if !is_fa
        reader = TwoBit.Reader(ref)
        chrom = copy(reader.names)
        for record in reader
            maxchromlen = maxchromlen > record.dnasize ? maxchromlen : record.dnasize
        end
    else
        chrom = Vector{String}()
        for record in FASTA.Reader(ref)
            if FASTA.hasidentifier(record)
                push!(chrom, FASTA.identifier(record))
                this_len = last(record.sequence) - first(record.sequence) + 1
                maxchromlen = maxchromlen > this_len ? maxchromlen : this_len
            end
        end
    end
    close(ref)
    return GenomeInfo(now(Dates.UTC), filepath, checksum, name, "",
                      chrom, maxchromlen, is_fa)
end

struct Locus{T<:Unsigned,K<:Unsigned}
    chrom::T
    pos::K
    isplus::Bool
end

function decode(locus::Locus, gI::GenomeInfo)
    strand = locus.isplus ? "+" : "-"
    return gI[locus.chrom] * ":" * string(locus.pos) * ":" * strand
end
