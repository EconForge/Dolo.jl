using Dolo
using Test

@testset verbose=true "Elements" begin
    # test_spaces.jl still uses the pre-unification API (CSpace with pairs)
    include("test_grids.jl")
    include("test_interpolation.jl")
    include("test_complementarities.jl")
    include("test_linsolve.jl")
    include("test_engines.jl")
end

@testset verbose=true "Models" begin
    include("time_iteration.jl")
end

