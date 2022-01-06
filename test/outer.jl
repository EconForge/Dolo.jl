using Random
using Distributions
using StaticArrays
using Dolo

@testset "testing the function outer of ergodic.jl" begin

    λ1, λ2 = rand(Uniform(0,1),2)

    @test Dolo.outer(SVector(1-λ1,λ1),SVector(1-λ2, λ2)) == @SMatrix[(1-λ1)*(1-λ2) (1-λ1)*λ2; λ1*(1-λ2) λ1*λ2]
    @test Dolo.outer(tuple( (SVector(1-λ1,λ1) for i in 1:1)... )...) == @SVector(1-λ1,λ1)

    
end
