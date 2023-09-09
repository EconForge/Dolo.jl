# test/runtests.jl
using Dolo
using Test

@testset verbose=true "Elements" begin
#     include("test_spaces.jl")
    include("test_grids.jl")
    include("test_interpolation.jl")
end

@testset verbose=true "Models" begin
    
    include("time_iteration.jl")
end

# @testset "General" begin
#     include("general.jl")
# end

# @testset "Special" begin
#     include("test_dev.jl")
# end