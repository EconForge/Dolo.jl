# GPU tests on Intel GPUs (oneAPI), run in their own environment (test/oneapi/Project.toml),
# which loads oneAPI and therefore the DolooneAPIExt extension:
#
#     julia --project=test/oneapi -e 'using Pkg; Pkg.instantiate()'
#     julia --project=test/oneapi test/oneapi/runtests.jl
#
# or from the main test suite with `Pkg.test("Dolo"; test_args=["oneapi"])`.
# The tests themselves are in test/gpu/gpu_tests.jl, shared with the other GPU backends.
# Many Intel GPUs (e.g. integrated Xe graphics, Arc) have no Float64 support: the tests then
# run in Float32 only.

using oneAPI
using Dolo
using Test

include(joinpath(@__DIR__, "..", "gpu", "gpu_tests.jl"))

"Does the oneAPI device support Float64?"
supports_float64() = oneAPI.oneL0.module_properties(oneAPI.device()).fp64flags &
    oneAPI.oneL0.ZE_DEVICE_MODULE_FLAG_FP64 == oneAPI.oneL0.ZE_DEVICE_MODULE_FLAG_FP64

@testset verbose=true "oneAPI" begin

    @test Base.get_extension(Dolo, :DolooneAPIExt) !== nothing

    if !oneAPI.functional()
        @info "oneAPI is not functional on this machine: skipping GPU tests."
    else
        oneAPI.allowscalar(false)   # any scalar indexing on the GPU is a bug
        precisions = supports_float64() ? (Float64, Float32) : (Float32,)
        @info "oneAPI device: $(oneAPI.oneL0.properties(oneAPI.device()).name), tests in $(join(precisions, ", "))"
        gpu_tests(oneArray; precisions)
    end

end
