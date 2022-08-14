using Test

using BioSequences
using ARTEMIS: shrink_and_expand!, is_in_range
using Combinatorics

@testset "db_fmi.jl" begin

    @testset "shrink_and_expand!" begin
        long = BitVector(zeros(10))
        counts = repeat([2], 5)
        short = shrink_and_expand!(long, counts)
        @test length(short) == 5
        @test length(long) == 10
        @test all(short .== 0)
        @test all(long .== 0)

        long[3] = true
        long[7] = true
        short = shrink_and_expand!(long, counts)
        @test short == [0, 1, 0, 1, 0]
        @test long == [0, 0, 1, 1, 0, 0, 1, 1, 0, 0]

        short = shrink_and_expand!(long, counts)
        @test short == [0, 1, 0, 1, 0]
        @test long == [0, 0, 1, 1, 0, 0, 1, 1, 0, 0]

        long[2] = true
        long[6] = true
        short = shrink_and_expand!(long, counts)
        @test short == [1, 1, 1, 1, 0]
        @test long == [1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
    end

    
    @testset "is_in_range" begin
        all_pos = UInt32.([100, 200, 300, 400, 500])
        pos = BitVector(zeros(5))
        pos_chrom_fwd_in_range = is_in_range(
            all_pos, pos, [1, 2, 99, 1000], false, true, 4, 1, 3, 4)
        @test all(pos_chrom_fwd_in_range) == 0
        pos = BitVector(ones(5))
        pos_chrom_fwd_in_range = is_in_range(
            all_pos, pos, [1, 2, 95, 1000], false, true, 4, 1, 3, 4)
        @test sum(pos_chrom_fwd_in_range) == 1
    end
end