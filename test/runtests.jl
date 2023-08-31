# test/runtests.jl
using Dolo
using Test

@testset verbose=true "Elements" begin
#     include("test_spaces.jl")
    include("test_grids.jl")
end

@testset "General" begin
    include("general.jl")
end

# @testset "Special" begin
#     include("test_dev.jl")
# end