"
Contains information about which genome was used.

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


"Instantiate GenomeInfo object from genome file path."
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
                tl = length(TwoBit.sequence(record))
                maxchromlen = maxchromlen >  tl ? maxchromlen : tl
            end
        end
    else
        chrom = Vector{String}()
        for record in FASTA.Reader(ref)
            seqlen = length(FASTA.sequence(LongDNA{4}, record))
            if (FASTA.identifier(record) != "") & (seqlen > 0)
                push!(chrom, FASTA.identifier(record))
                maxchromlen = maxchromlen > seqlen ? maxchromlen : seqlen
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


Motif defines what genome file is being used for the searches.


# Arguments
`filepath` - Path to the genome file, if file is fasta (ends with .fa or .fasta or .fna) 
make sure you also have fasta index file with extension .fai. Alternatively, you can use .2bit 
genome file.  

`name` - Your name for this instance of DBInfo: the genome with connection to the motif and vcf file.

`motif`   - `Motif` object defining search parameters

`vcf_filepath`  - Optional. Path to the VCF file to include in the searches.


Alignments will be performed from opposite to the extension direction (which is defined by extend5).

# Examples
```julia
# use CHOPOFF example genome
genome = joinpath(vcat(splitpath(dirname(pathof(CHOPOFF)))[1:end-1], 
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


struct Offtarget{T<:Unsigned,K<:Unsigned}
    loc::Loc{T,K}
    dist::Int
    aln_guide::String
    aln_ref::String
end

function Base.show(io::IO, m::Offtarget)
    show(io, "$(m.loc.chrom):$(m.loc.isplus ? '+' : '-'):$(m.loc.pos); $(m.dist)")
end

import Base: isless
isless(x::Offtarget, y::Offtarget) = isless(
    (x.loc.chrom, x.loc.isplus, x.loc.pos, x.dist), 
    (y.loc.chrom, y.loc.isplus, y.loc.pos, y.dist))

"""x has to be >= to y"""
function is_in_range(x::Offtarget, y::Offtarget, range::Int)
    return x.loc.chrom == y.loc.chrom && x.loc.isplus == y.loc.isplus && ((x.loc.pos - y.loc.pos) <= range)
end

"""
Insert to vector checking which offtarget is overlapping which, so that 
we can filter out those almost same site off-targets, we take the one with smallest distance.
This returns two values, first indicates distance of the offtarget that it replaced (or nothing)
second returns distance of the offtraget it added (or nothing).
"""
function insert_offtarget!(x::Vector{Offtarget}, offt::Offtarget, range::Int)
    if isempty(x)
        push!(x, offt)
        return nothing, offt.dist
    end
    # idx of first value in a >= x, if x >= all in a lastindex(a) + 1
    idx = searchsortedfirst(x, offt)

    if idx == 1
        le = x[1]
        if is_in_range(le, offt, range)
            if le.dist > offt.dist
                x[1] = offt
                return le.dist, offt.dist
            else
                return nothing, nothing
            end
        end
    elseif idx == (lastindex(x) + 1)
        le = x[end]
        if is_in_range(offt, le, range)
            if le.dist > offt.dist
                x[end] = offt
                return le.dist, offt.dist
            else
                return nothing, nothing
            end
        end
    else
        le = x[idx - 1] # smaller than offt
        le2 = x[idx] # larger or equal
        in_range_le = is_in_range(offt, le, range)
        in_range_le2 = is_in_range(le2, offt, range)
        if in_range_le && le.dist > offt.dist
            x[idx - 1] = offt
            return le.dist, offt.dist
        elseif in_range_le2 && le2.dist > offt.dist
            x[idx] = offt
            return le2.dist, offt.dist
        elseif in_range_le || in_range_le2 # in range but distance is no good
            return nothing, nothing
        end
    end

    # when not in range of anything we just insert
    insert!(x, idx, offt)
    return nothing, offt.dist
end