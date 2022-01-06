using Random
using Distributions
using StaticArrays
using Dolo

@testset "Tests on outer product functions" begin

    @testset "testing the function outer of ergodic.jl" begin

        λ1, λ2 = rand(Uniform(0, 1), 2)

        @test Dolo.outer(SVector(1 - λ1, λ1), SVector(1 - λ2, λ2)) == @SMatrix [(1-λ1)*(1-λ2) (1-λ1)*λ2; λ1*(1-λ2) λ1*λ2]
        @test Dolo.outer(SVector(1 - λ1, λ1)) == @SVector(1 - λ1, λ1)

    end

    @testset "testing the function outer2 of ergodic.jl" begin
        x = @SMatrix randn(4,4)
        @test Dolo.outer2(SVector(-1.,0,1.),x) == SVector(-1*x, 0*x, 1*X)
    end
end

