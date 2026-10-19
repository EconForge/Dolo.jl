# GPU code paths, tested on the CPU with JLArrays: a reference implementation of a GPU array
# (scalar indexing disallowed, kernels through KernelAbstractions), with no GPU needed.
# Only the operations are tested here; the complete GPU tests (with solvers) also pass with
# JLArrays (`gpu_tests(JLArray)`, about 2 minutes), and run on real GPUs in test/cuda, test/oneapi.

using JLArrays

JLArrays.allowscalar(false)

# JLArrays operations are synchronous
Dolo.KernelAbstractions.synchronize(::JLArrays.JLBackend) = nothing

include(joinpath(@__DIR__, "gpu", "gpu_tests.jl"))

@testset "GPU code paths (JLArrays)" begin
    gpu_tests(JLArray; precisions=(Float64, Float32), solvers=false)
end
