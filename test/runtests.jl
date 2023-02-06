using ModularSquareRoots
using Test

@testset "ModularSquareRoots.jl" begin

    # @test sqrtmod(13, 1234566) == [95465, 232639, 533909, 563483, 671083, 700657, 1001927, 1139101]


    @testset "sqrtmodp" begin
        @test sort(sqrtmodp(0, 2)) == [0]
        @test sort(sqrtmodp(1, 2)) == [1]

        @test sort(sqrtmodp(0, 3)) == [0]
        @test sort(sqrtmodp(1, 3)) == [1, 2]
        @test isempty(sqrtmodp(2, 3))

        @test sqrtmodp(Int8(125), Int8(127)) isa Vector{Int8}
        @test sqrtmodp(Int8(125), Int8(127)) isa Vector{Int8}
        @test sqrtmodp(Int8(125), Int8(113)) isa Vector{Int8}
        @test sqrtmodp(Int8(126), Int8(113)) isa Vector{Int8}
        @test sqrtmodp(zero(Int8), Int8(113)) isa Vector{Int8}
        @test isempty(sqrtmodp(Int8(125), Int8(127)))
        @test isempty(sqrtmodp(Int8(125), Int8(113)))
        @test sort(sqrtmodp(Int8(124), Int8(127))) == [39, 88]
        @test sort(sqrtmodp(Int8(126), Int8(113))) == [37, 76]
        @test sqrtmodp(zero(Int8), Int8(113)) == [0]

        @test sqrtmodp(UInt16(1266), UInt16(65521)) isa Vector{UInt16}
        @test sqrtmodp(UInt16(1267), UInt16(65521)) isa Vector{UInt16}
        @test isempty(sqrtmodp(UInt16(1266), UInt16(65521)))
        @test sort(sqrtmodp(UInt16(1267), UInt16(65521))) == [16185, 49336]

        @test sort(sqrtmodp(13, 10^16 + 61)) == [353670011380822, 9646329988619239]
        @test sort(sqrtmodp(483918293829192838, 5920394019203941009)) == [835763758593320363, 5084630260610620646]
        @test isempty(sqrtmodp(483918293829192839, 5920394019203941009))

        @test sort(sqrtmodp(BigInt(34789074890172839471892043422), BigInt(10000000000000000000000000120989089031284027))) == [4665676551146615114818820594052041583630135, 5334323448853384885181179526937047447653892]
        @test sort(sqrtmodp(32138901234789074890172839471892043424434, 57348293847897092349345713247809314579049)) == [25770511472149243631404498078881896545120, 31577782375747848717941215168927418033929]
    end
end
