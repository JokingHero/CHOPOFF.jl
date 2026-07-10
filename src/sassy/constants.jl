# --- Constants ---
const LANES = 4
const BLOCK_SIZE = 64

# IUPAC encoding 
const IUPAC_SASSY = let
    enc = zeros(UInt8, 32)
    A = UInt8(1 << 0)
    C = UInt8(1 << 1)
    T = UInt8(1 << 2)
    G = UInt8(1 << 3)
    enc[Int('A') & 0x1F + 1] = A
    enc[Int('C') & 0x1F + 1] = C
    enc[Int('T') & 0x1F + 1] = T
    enc[Int('U') & 0x1F + 1] = T
    enc[Int('G') & 0x1F + 1] = G
    enc[Int('N') & 0x1F + 1] = A | C | T | G
    enc[Int('R') & 0x1F + 1] = A | G
    enc[Int('Y') & 0x1F + 1] = C | T
    enc[Int('S') & 0x1F + 1] = G | C
    enc[Int('W') & 0x1F + 1] = A | T
    enc[Int('K') & 0x1F + 1] = G | T
    enc[Int('M') & 0x1F + 1] = A | C
    enc[Int('B') & 0x1F + 1] = C | G | T
    enc[Int('D') & 0x1F + 1] = A | G | T
    enc[Int('H') & 0x1F + 1] = A | C | T
    enc[Int('V') & 0x1F + 1] = A | C | G
    enc[Int('X') & 0x1F + 1] = 0
    enc
end

# 256-element Branchless LUT
const IUPAC_TABLE = let
    t = zeros(UInt8, 256)
    for i in 0:255
        c = UInt8(i)
        if (UInt8('A') <= c <= UInt8('Z')) || (UInt8('a') <= c <= UInt8('z'))
            t[i+1] = IUPAC_SASSY[(c & 0x1F) + 1]
        end
    end
    t
end

@inline function get_iupac_mask(c::UInt8)
    return @inbounds IUPAC_TABLE[c + 1]
end

const REF_AMBIG_TABLE = let
    t = ones(UInt8, 256)
    for c in (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T'),
              UInt8('a'), UInt8('c'), UInt8('g'), UInt8('t'))
        t[c + 1] = UInt8(0)
    end
    t
end

@inline function is_ref_ambig(c::UInt8)
    return @inbounds REF_AMBIG_TABLE[c + 1] != 0
end

# Pre-expanded per-byte → per-base match LUT (256 × 4 UInt8)
# Avoids repeated get_iupac_mask + bitwise AND in the 4-base fast path.
const BASE_MATCH = let
    t = zeros(UInt8, 256, 4)
    masks = (UInt8(0x01), UInt8(0x02), UInt8(0x04), UInt8(0x08))  # A, C, T, G
    for i in 0:255
        m = IUPAC_TABLE[i + 1]
        for j in 1:4
            t[i + 1, j] = UInt8((m & masks[j]) != 0)
        end
    end
    t
end

# BioSequences DNAAlphabet{4} nibble → ASCII lookup
# BioSymbols encoding: A=0x01, C=0x02, G=0x04, T=0x08
const NIBBLE_TO_ASCII = UInt8[
    UInt8('-'),  # 0x00 gap
    UInt8('A'),  # 0x01
    UInt8('C'),  # 0x02
    UInt8('M'),  # 0x03 A|C
    UInt8('G'),  # 0x04
    UInt8('R'),  # 0x05 A|G
    UInt8('S'),  # 0x06 C|G
    UInt8('V'),  # 0x07 A|C|G
    UInt8('T'),  # 0x08
    UInt8('W'),  # 0x09 A|T
    UInt8('Y'),  # 0x0A C|T
    UInt8('H'),  # 0x0B A|C|T
    UInt8('K'),  # 0x0C G|T
    UInt8('D'),  # 0x0D A|G|T
    UInt8('B'),  # 0x0E C|G|T
    UInt8('N'),  # 0x0F A|C|G|T
]