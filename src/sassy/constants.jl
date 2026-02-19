# --- Constants ---
const LANES = 4
const BLOCK_SIZE = 64

# IUPAC encoding (Matches Rust sassy/profiles/iupac.rs)
const IUPAC_SASSY = let
    enc = zeros(UInt8, 32)
    A = UInt8(1 << 0)
    C = UInt8(1 << 1)
    T = UInt8(1 << 2)
    G = UInt8(1 << 3)
    # ASCII & 0x1F gives 1-based index (A=1, C=3, etc.)
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

@inline function get_iupac_mask(c::UInt8)
    # Ensure only letters are mapped to IUPAC (ASCII 65-90, 97-122)
    # This prevents non-letters like '-' from colliding with IUPAC codes.
    if (UInt8('A') <= c <= UInt8('Z')) || (UInt8('a') <= c <= UInt8('z'))
        return @inbounds IUPAC_SASSY[(c & 0x1F) + 1]
    end
    return UInt8(0)
end