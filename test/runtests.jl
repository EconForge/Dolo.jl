using Dolo
using Test

@testset verbose=true "Elements" begin
    # test_spaces.jl still uses the pre-unification API (CSpace with pairs)
    include("test_grids.jl")
    include("test_interpolation.jl")
    include("test_prefilter.jl")
    include("test_complementarities.jl")
    include("test_linsolve.jl")
    include("test_engines.jl")
end

@testset verbose=true "Models" begin
    include("time_iteration.jl")
end

# GPU tests, only on request. They run in their own environment (test/cuda), which loads
# CUDA and hence the CUDA extension:
#
#     Pkg.test("Dolo"; test_args=["cuda"])     or     DOLO_TEST_CUDA=true
if "cuda" in ARGS || get(ENV, "DOLO_TEST_CUDA", "false") == "true"
    @testset "CUDA" begin
        env = joinpath(@__DIR__, "cuda")
        julia = Base.julia_cmd()
        withenv("JULIA_LOAD_PATH" => nothing, "JULIA_PROJECT" => nothing) do
            @test success(run(`$julia --project=$env -e 'using Pkg; Pkg.instantiate()'`))
            @test success(run(ignorestatus(`$julia --project=$env $(joinpath(env, "runtests.jl"))`)))
        end
    end
end
