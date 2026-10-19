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
    include("test_jlarrays.jl")
end

@testset verbose=true "Models" begin
    include("time_iteration.jl")
end

# GPU tests, only on request. They run in their own environment (test/cuda or test/oneapi),
# which loads the GPU package and hence Dolo's extension for it:
#
#     Pkg.test("Dolo"; test_args=["cuda"])       or     DOLO_TEST_CUDA=true
#     Pkg.test("Dolo"; test_args=["oneapi"])     or     DOLO_TEST_ONEAPI=true
for (name, var) in (("cuda", "DOLO_TEST_CUDA"), ("oneapi", "DOLO_TEST_ONEAPI"))
    if name in ARGS || get(ENV, var, "false") == "true"
        @testset "$name" begin
            env = joinpath(@__DIR__, name)
            julia = Base.julia_cmd()
            withenv("JULIA_LOAD_PATH" => nothing, "JULIA_PROJECT" => nothing) do
                @test success(run(`$julia --project=$env -e 'using Pkg; Pkg.instantiate()'`))
                @test success(run(ignorestatus(`$julia --project=$env $(joinpath(env, "runtests.jl"))`)))
            end
        end
    end
end
