using Test
using CHOPOFF

# Access sassy internal functions via qualified name
const TextBlockProfile = CHOPOFF.TextBlockProfile
const encode_text_block = CHOPOFF.encode_text_block
const get_iupac_mask = CHOPOFF.get_iupac_mask

"""
Helper function to get match positions from a u64 mask.
Returns vector of bit positions (0-indexed) that are set to 1.
"""
function get_match_positions_u64(mask::UInt64)::Vector{Int}
    positions = Int[]
    for pos in 0:63
        if (mask & (UInt64(1) << pos)) != 0
            push!(positions, pos)
        end
    end
    return positions
end

"""
Helper function to extract match positions from all masks in a TextBlockProfile.
"""
function get_match_positions(profile::TextBlockProfile)::Vector{Vector{Int}}
    [get_match_positions_u64(mask) for mask in profile.masks]
end

@testset "IUPAC encoding - Basic ATGC" begin
    # Test from Rust test_just_atgc
    # seq[1] = 'a', seq[2] = 'y' (C or T), rest = 'g'
    seq = UInt8[UInt8('g'), UInt8('a'), UInt8('y'), UInt8('G'), UInt8('G'), UInt8('G'), UInt8('G'), UInt8('G')]
    # Pad to 64 with 'g'
    while length(seq) < 64
        push!(seq, UInt8('g'))
    end

    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])
    positions = get_match_positions(profile)

    # A matches at position 1 (0-indexed: 1)
    @test positions[1] == [1]
    # T matches at position 2 (Y is C|T, so T matches Y) (0-indexed: 2)
    @test positions[3] == [2]
    # G matches at positions 0, 3-63 (64 - 2 = 62 positions)
    @test length(positions[4]) == 62
    @test 0 in positions[4]  # First position (0-indexed)
    # C matches at position 2 (Y is C|T) (0-indexed: 2)
    @test positions[2] == [2]
end

@testset "IUPAC encoding - Extra bases NY" begin
    # Test from Rust test_extra_bases_ny
    # seq[1] = 'a', seq[2] = 'y', seq[3] = 'C', rest = 'g'
    seq = UInt8[UInt8('g'), UInt8('a'), UInt8('y'), UInt8('C'), UInt8('G'), UInt8('G')]
    while length(seq) < 64
        push!(seq, UInt8('g'))
    end

    # Add N and Y as extra bases
    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G'), UInt8('N'), UInt8('Y')])
    positions = get_match_positions(profile)

    # N matches all positions
    @test length(positions[5]) == 64
    # Y matches positions 2 and 3 (0-indexed): y at pos 2, C at pos 3
    @test positions[6] == [2, 3]
end

@testset "IUPAC encoding - Full 64-byte block" begin
    # Test from Rust test_just_atgc_64
    seq = fill(UInt8('g'), 64)
    seq[2] = UInt8('a')   # 1-indexed: position 2, 0-indexed bit: 1
    seq[3] = UInt8('y')   # 1-indexed: position 3, 0-indexed bit: 2
    seq[35] = UInt8('y')  # 1-indexed: position 35, 0-indexed bit: 34

    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])
    positions = get_match_positions(profile)

    # A at seq[2] = bit 1 (0-indexed)
    @test positions[1] == [1]
    # T at positions 3 and 35 (0-indexed bits: 2, 34)
    @test Set(positions[3]) == Set([2, 34])
    # G at all other positions
    @test 0 in positions[4]
    @test 3 in positions[4]
    @test !(2 in positions[4])  # Bit 2 is y, not G
    @test !(34 in positions[4])  # Bit 34 is y, not G
end

@testset "IUPAC encoding - Extra bases NY 64" begin
    # Test from Rust test_extra_bases_ny_64
    seq = fill(UInt8('g'), 64)
    seq[2] = UInt8('a')   # bit 1
    seq[3] = UInt8('y')   # bit 2
    seq[4] = UInt8('C')   # bit 3
    seq[51] = UInt8('y')  # bit 50
    seq[64] = UInt8('y')  # bit 63

    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G'), UInt8('N'), UInt8('Y')])
    positions = get_match_positions(profile)

    # N matches all
    @test length(positions[5]) == 64
    # Y matches at bits 2, 3, 50, 63
    @test Set(positions[6]) == Set([2, 3, 50, 63])
end

@testset "IUPAC encoding - Case insensitive" begin
    # Test from Rust test_iupac_u64_case_insensitive
    seq = fill(UInt8('G'), 64)
    seq[2] = UInt8('a')   # bit 1
    seq[3] = UInt8('A')   # bit 2
    seq[5] = UInt8('r')   # bit 4 (R = A|G)
    seq[6] = UInt8('W')   # bit 5 (W = A|T)

    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])
    positions = get_match_positions(profile)

    # A should match at bits 1, 2, 4, 5
    @test Set(positions[1]) == Set([1, 2, 4, 5])
end

@testset "IUPAC encoding - X returns 0" begin
    seq = vcat(UInt8[UInt8('X'), UInt8('X'), UInt8('X'), UInt8('X')], fill(UInt8('A'), 60))

    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])

    # First 4 bits (0, 1, 2, 3) should NOT match A (X has mask 0)
    a_positions = get_match_positions_u64(profile.masks[1])
    @test !(0 in a_positions)
    @test !(1 in a_positions)
    @test !(2 in a_positions)
    @test !(3 in a_positions)
end

@testset "IUPAC encoding - R (A|G) and S (G|C)" begin
    seq = UInt8[UInt8('A'), UInt8('G'), UInt8('R'), UInt8('S'), UInt8('C')]
    while length(seq) < 64
        push!(seq, UInt8('T'))
    end

    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])
    positions = get_match_positions(profile)

    # A matches at bits 0, 2 (R includes A)
    @test Set(positions[1]) == Set([0, 2])
    # C matches at bits 3, 4 (S includes C)
    @test Set(positions[2]) == Set([3, 4])
    # G matches at bits 1, 2, 3 (R includes G, S includes G)
    @test Set(positions[4]) == Set([1, 2, 3])
end

@testset "IUPAC encoding - Short sequence" begin
    seq = UInt8[UInt8('A'), UInt8('T'), UInt8('G'), UInt8('C')]
    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])

    @test profile.masks[1] == 0b0001  # A at pos 1 (bit 0)
    @test profile.masks[2] == 0b1000  # C at pos 4 (bit 3)
    @test profile.masks[3] == 0b0010  # T at pos 2 (bit 1)
    @test profile.masks[4] == 0b0100  # G at pos 3 (bit 2)
end

@testset "IUPAC encoding - N matches all" begin
    seq = UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G'), UInt8('N'), UInt8('R')]
    profile = encode_text_block(seq, UInt8[UInt8('N')])

    # N should match all 6 positions
    @test profile.masks[1] == 0b111111
end

@testset "IUPAC encoding - Complex ambiguity codes" begin
    # W = A|T, K = G|T, M = A|C, B = C|G|T, D = A|G|T, H = A|C|T, V = A|C|G
    seq = UInt8[UInt8('W'), UInt8('K'), UInt8('M'), UInt8('B'), UInt8('D'), UInt8('H'), UInt8('V')]
    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])

    # W (pos 1) = A|T
    @test (profile.masks[1] & 0b0000001) != 0  # A at bit 0
    @test (profile.masks[3] & 0b0000001) != 0  # T at bit 0

    # K (pos 2) = G|T
    @test (profile.masks[3] & 0b0000010) != 0  # T at bit 1
    @test (profile.masks[4] & 0b0000010) != 0  # G at bit 1

    # M (pos 3) = A|C
    @test (profile.masks[1] & 0b0000100) != 0  # A at bit 2
    @test (profile.masks[2] & 0b0000100) != 0  # C at bit 2
end

@testset "IUPAC get_iupac_mask function" begin
    # Test individual IUPAC masks
    @test get_iupac_mask(UInt8('A')) == 0x01
    @test get_iupac_mask(UInt8('C')) == 0x02
    @test get_iupac_mask(UInt8('T')) == 0x04
    @test get_iupac_mask(UInt8('G')) == 0x08
    @test get_iupac_mask(UInt8('N')) == 0x0F  # All
    @test get_iupac_mask(UInt8('Y')) == 0x06  # C|T = 0x02|0x04
    @test get_iupac_mask(UInt8('R')) == 0x09  # A|G = 0x01|0x08
    @test get_iupac_mask(UInt8('X')) == 0x00  # None

    # Case insensitivity
    @test get_iupac_mask(UInt8('a')) == 0x01
    @test get_iupac_mask(UInt8('c')) == 0x02
end

@testset "IUPAC encoding - Lowercase input" begin
    seq = UInt8[UInt8('a'), UInt8('c'), UInt8('t'), UInt8('g'), UInt8('n'), UInt8('r'), UInt8('y')]
    profile = encode_text_block(seq, UInt8[UInt8('A'), UInt8('C'), UInt8('T'), UInt8('G')])

    # All lowercase should work
    # a matches A at bit 0
    # n (A|C|T|G) matches A at bit 4
    # r (A|G) matches A at bit 5
    @test profile.masks[1] == 0b110001  # bits 0, 4, 5
end
