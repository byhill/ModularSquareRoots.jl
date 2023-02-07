using ModularSquareRoots
using Test

@testset "ModularSquareRoots.jl" begin

    @testset "sqrtmodprime" begin
        @test sort(sqrtmodprime(0, 2)) == [0]
        @test sort(sqrtmodprime(1, 2)) == [1]

        @test sort(sqrtmodprime(0, 3)) == [0]
        @test sort(sqrtmodprime(1, 3)) == [1, 2]
        @test isempty(sqrtmodprime(2, 3))

        @test sqrtmodprime(Int8(125), Int8(127)) isa Vector{Int8}
        @test sqrtmodprime(Int8(125), Int8(127)) isa Vector{Int8}
        @test sqrtmodprime(Int8(125), Int8(113)) isa Vector{Int8}
        @test sqrtmodprime(Int8(126), Int8(113)) isa Vector{Int8}
        @test sqrtmodprime(zero(Int8), Int8(113)) isa Vector{Int8}
        @test isempty(sqrtmodprime(Int8(125), Int8(127)))
        @test isempty(sqrtmodprime(Int8(125), Int8(113)))
        @test sort(sqrtmodprime(Int8(124), Int8(127))) == [39, 88]
        @test sort(sqrtmodprime(Int8(126), Int8(113))) == [37, 76]
        @test sqrtmodprime(zero(Int8), Int8(113)) == [0]

        @test sqrtmodprime(UInt16(1266), UInt16(65521)) isa Vector{UInt16}
        @test sqrtmodprime(UInt16(1267), UInt16(65521)) isa Vector{UInt16}
        @test isempty(sqrtmodprime(UInt16(1266), UInt16(65521)))
        @test sort(sqrtmodprime(UInt16(1267), UInt16(65521))) == [16185, 49336]

        @test sort(sqrtmodprime(13, 10^16 + 61)) == [353670011380822, 9646329988619239]
        @test sort(sqrtmodprime(483918293829192838, 5920394019203941009)) == [835763758593320363, 5084630260610620646]
        @test isempty(sqrtmodprime(483918293829192839, 5920394019203941009))

        @test sort(sqrtmodprime(BigInt(34789074890172839471892043422), BigInt(10000000000000000000000000120989089031284027))) == [4665676551146615114818820594052041583630135, 5334323448853384885181179526937047447653892]
        @test sort(sqrtmodprime(32138901234789074890172839471892043424434, 57348293847897092349345713247809314579049)) == [25770511472149243631404498078881896545120, 31577782375747848717941215168927418033929]
    end

    @testset "sqrtmod" begin
        @test_throws DomainError sqrtmod(4, 0)
        @test_throws DomainError sqrtmod(UInt64(4), zero(UInt64))
        @test_throws DomainError sqrtmod(12, -1)

        @test sqrtmod(Int16(4), one(Int32)) isa Vector{Int32}
        @test sqrtmod(4, one(Int128)) == [0]

        @test sort(sqrtmod(0, 2)) == [0]
        @test sort(sqrtmod(1, 2)) == [1]

        @test sort(sqrtmod(0, 3)) == [0]
        @test sort(sqrtmod(1, 3)) == [1, 2]
        @test isempty(sqrtmod(2, 3))

        @test sqrtmod(UInt128(1240), Int128(289032)) isa Vector{UInt128}
        @test sort(sqrtmod(UInt128(1240), Int128(289032))) == [10712, 37460, 107056, 133804, 155228, 181976, 251572, 278320]
        @test sort(sqrtmod(13, 1234566)) == [95465, 232639, 533909, 563483, 671083, 700657, 1001927, 1139101]
    end

end
