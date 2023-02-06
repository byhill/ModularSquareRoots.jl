using ModularSquareRoots
using Test

@testset "ModularSquareRoots.jl" begin

    @test sqrtmod(13, 1234566) == [95465, 232639, 533909, 563483, 671083, 700657, 1001927, 1139101]
end
