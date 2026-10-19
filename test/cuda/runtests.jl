# GPU tests on CUDA, run in their own environment (test/cuda/Project.toml), which loads CUDA
# and therefore the DoloCUDAExt extension:
#
#     julia --project=test/cuda -e 'using Pkg; Pkg.instantiate()'
#     julia --project=test/cuda test/cuda/runtests.jl
#
# or from the main test suite with `Pkg.test("Dolo"; test_args=["cuda"])`.
# The tests themselves are in test/gpu/gpu_tests.jl, shared with the other GPU backends.

using CUDA
using Dolo
using Test

include(joinpath(@__DIR__, "..", "gpu", "gpu_tests.jl"))

@testset verbose=true "CUDA" begin

    @test Base.get_extension(Dolo, :DoloCUDAExt) !== nothing

    if !CUDA.functional()
        @info "CUDA is not functional on this machine: skipping GPU tests."
    else
        CUDA.allowscalar(false)   # any scalar indexing on the GPU is a bug
        gpu_tests(CuArray; precisions=(Float64, Float32))
    end

end
