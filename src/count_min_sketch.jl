using Probably
using BioSequences
using FASTX
using Serialization

using BenchmarkTools
# what is the time cost of hashing Uint128 vs dna
# not that important?!
@benchmark hash(UInt128(1231383193))


# Step 1 - Can we use HyperLogLog for estimating
# unique count of guides?
# hll = HyperLogLog{18}()

# SPEEDUP strategy
# Paralelize HLL on each chromosome + separate HLL and merge with union
#query = ExactSearchQuery(dna"NNGG")
#findfirst(query, dna"ATTGGTTACTGCATCGACTAGG")
## can be run with the positions consecutively

# using CodecZlib.jl
# FASTA.Reader(GzipDecompressorStream(open("my-reads.fasta.gz")))

#reader = FASTA.Reader(open("/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa", "r"))
#record = FASTA.Record()
#read!(reader, record)
#close(reader)

#for record in reader
#    @show FASTA.identifier(record)

#    for x in eachmatch(biore"NNNNNNNNNNNNNNNNNNNNNGG"d, FASTA.sequence(record))
#        seq = matched(x)
#        if !hasambiguity(seq)
#            push!(hll, seq[1:end-3])
#        end
#    end

    # reverse
#    for x in eachmatch(biore"CCNNNNNNNNNNNNNNNNNNNNN"d, FASTA.sequence(record))
#        seq = matched(x)
#        if !hasambiguity(seq)
#            push!(hll, seq[4:end])
#        end
#    end
#end
#close(reader)

max_len = 229145118 # length(hll)
# very close to the real value of 230132199 # not bad!


# Step 2
# build sketch for human genome
# test it works on 0 mm
# and 1 MM

# 0MM
probability_of_error = 0.01 # lets try this

# we compute how many hashing functions we need
depth = ceil(log(1/probability_of_error))
# estimate our error E based on the max_len
# we add 0.5% error of estimation for HLL
E = 1/(max_len + 0.005 * max_len)

# we need to simulate errors
width = ceil(Base.ℯ/E)

# UInt16 -> 65535
sketch0 = CountMinSketch{UInt8}(width, depth)
Base.summarysize(sketch0)/1e+6 # in mb

# make size (in mb) vs error (maximum error boxplots)
# code faster FASTA scan
# copy & adjust CMS code
# figure out transformations for deletes
# figure out hashes - are we using correct hashing method?!

# TwoBit.Reader
# seqnames(reader)
# TwoBit.Record()

reader = FASTA.Reader(open("/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa", "r"))
for record in reader
    @show FASTA.identifier(record)
    chrom = FASTA.sequence(record)

    for x in eachmatch(biore"NNNNNNNNNNNNNNNNNNNNNGG"d, chrom)
        seq = matched(x)
        if !hasambiguity(seq)
            push!(sketch0, seq[1:end-3])
        end
    end

    # reverse
    for x in eachmatch(biore"CCNNNNNNNNNNNNNNNNNNNNN"d, chrom)
        seq = matched(x)
        if !hasambiguity(seq)
            push!(sketch0, seq[4:end])
        end
    end
end
close(reader)

# test
function fprof_fixed(sketch::CountMinSketch)
    rate = 1
    for col in 1:sketch.width
        full_in_row = 0
        for row in 1:sketch.len
            full_in_row += sketch.matrix[row, col] > zero(eltype(sketch))
        end
        rate *= full_in_row / sketch.len
    end
    return rate
end

fprof_fixed(sketch0)
# how filled up is our sketch
# FIll ScaleE ProbErr MB
# 0.88 10     0.01    312
#      10


# 0 MM guides
# CGGTGTTCCGAGGCAGGACGGGG
# GTAGAAGTCCGCCGCCACGCTGG
# GGCGGACTTCTACCACCGGCTGG
# GCTCTACGCCATGCACCCGTGGG
# GCAAGGGCGCTATCATGGTATGG

# CCTGAGGAAGTTGATCTCGTCGG chr12:52900603  + 0       1       1       4
# GCGGAATGAATGGGGTGAGCTGG chr12:52926432  - 0       0       0       11
# GGGAGGCATCACCGCAGTTACGG chr12:52904780  - 0       2       4       8
# AGAAGATCGAGACACGTGATGGG chr12:52897467  - 1       3       6       16

# Other
# TTTTGCAGGTGAGAACACCAAGG 0	3 47 >=259
# TTTTTTTTTTTGTCTTTTTTAGG 1	142	674	>=0
# ACCTGTAATCCCAGCTACTCGGG 1200 0 0 >=0
Int(sketch0[dna"CGGTGTTCCGAGGCAGGACG"])
Int(sketch0[dna"GTAGAAGTCCGCCGCCACGC"])
Int(sketch0[dna"GGCGGACTTCTACCACCGGC"])
Int(sketch0[dna"GCTCTACGCCATGCACCCGT"])
Int(sketch0[dna"GCAAGGGCGCTATCATGGTA"])
Int(sketch0[dna"GGGAGGCATCACCGCAGTTA"])
Int(sketch0[dna"ACCTGTAATCCCAGCTACTC"])
Int(sketch0[dna"TTTTTTTTTTTGTCTTTTTT"])
