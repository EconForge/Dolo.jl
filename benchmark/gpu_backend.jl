# GPU backend of the benchmark scripts, chosen with `--backend=cuda` (default) or `--backend=oneapi`.
# The corresponding package must be in the active environment: benchmark/ for CUDA,
# benchmark/oneapi/ for oneAPI.
#
# Defines `GPU`, with:
#   name, AT (array type), backend (KernelAbstractions backend), functional, float64 (support),
#   sync() (wait for the GPU), reclaim() (free cached memory), description

backend_option = let a = findfirst(startswith("--backend="), ARGS)
    a === nothing ? "cuda" : split(ARGS[a], "=")[2]
end

if backend_option == "cuda"
    using CUDA
    const GPU = (;
        name = "cuda",
        AT = CuArray,
        backend = CUDA.functional() ? CUDABackend() : nothing,
        functional = CUDA.functional(),
        float64 = true,
        sync = CUDA.synchronize,
        reclaim = CUDA.reclaim,
        description = CUDA.functional() ?
            "$(CUDA.name(CUDA.device())), $(round(CUDA.totalmem(CUDA.device())/2^30, digits=1)) GiB" : "",
    )
elseif backend_option == "oneapi"
    using oneAPI
    const GPU = let functional = oneAPI.functional()
        dev = functional ? oneAPI.device() : nothing
        (;
            name = "oneapi",
            AT = oneArray,
            backend = functional ? oneAPIBackend() : nothing,
            functional,
            # many Intel GPUs (integrated Xe, Arc) have no Float64
            float64 = functional && oneAPI.oneL0.module_properties(dev).fp64flags &
                oneAPI.oneL0.ZE_DEVICE_MODULE_FLAG_FP64 == oneAPI.oneL0.ZE_DEVICE_MODULE_FLAG_FP64,
            sync = oneAPI.synchronize,
            reclaim = () -> GC.gc(true),
            description = functional ? "$(oneAPI.oneL0.properties(dev).name)" : "",
        )
    end
else
    error("Unknown backend $backend_option: use --backend=cuda or --backend=oneapi.")
end

GPU.functional || @warn "$(GPU.name) is not functional on this machine: GPU benchmarks are skipped."
GPU.functional && !GPU.float64 && @info "$(GPU.description) has no Float64 support: GPU benchmarks in Float32 only."
