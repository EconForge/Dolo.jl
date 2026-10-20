# GPU backend of the benchmark scripts, chosen with `--backend=cuda` (default) or `--backend=oneapi`.
# The corresponding package must be in the active environment: benchmark/ for CUDA,
# benchmark/oneapi/ for oneAPI.
#
# Defines `GPU`, with:
#   name, AT (array type), backend (KernelAbstractions backend), functional, float64 (support),
#   sync() (wait for the GPU), reclaim() (free cached memory), description,
#   platform and software (dictionaries saved with the results, see record.jl)

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
        platform = CUDA.functional() ? Dict{String,Any}(
            "name" => CUDA.name(CUDA.device()),
            "memory_GiB" => round(CUDA.totalmem(CUDA.device())/2^30, digits=1),
            "capability" => string(CUDA.capability(CUDA.device())),
            "float64" => true,
        ) : Dict{String,Any}(),
        software = CUDA.functional() ? Dict{String,Any}(
            "backend" => "cuda",
            "driver" => string(CUDA.driver_version()),
            "runtime" => string(CUDA.runtime_version()),
        ) : Dict{String,Any}(),
    )
elseif backend_option == "oneapi"
    using oneAPI
    const GPU = let functional = oneAPI.functional()
        dev = functional ? oneAPI.device() : nothing
        # many Intel GPUs (integrated Xe, Arc) have no Float64
        float64 = functional && oneAPI.oneL0.module_properties(dev).fp64flags &
            oneAPI.oneL0.ZE_DEVICE_MODULE_FLAG_FP64 == oneAPI.oneL0.ZE_DEVICE_MODULE_FLAG_FP64
        (;
            name = "oneapi",
            AT = oneArray,
            backend = functional ? oneAPIBackend() : nothing,
            functional,
            float64,
            sync = oneAPI.synchronize,
            reclaim = () -> GC.gc(true),
            description = functional ? "$(oneAPI.oneL0.properties(dev).name)" : "",
            platform = functional ? Dict{String,Any}(
                "name" => oneAPI.oneL0.properties(dev).name,
                "memory_GiB" => round(sum(p.totalSize for p in oneAPI.oneL0.memory_properties(dev))/2^30, digits=1),
                "float64" => float64,
            ) : Dict{String,Any}(),
            software = functional ? Dict{String,Any}(
                "backend" => "oneapi",
                "driver" => string(oneAPI.oneL0.properties(oneAPI.driver()).driverVersion),
            ) : Dict{String,Any}(),
        )
    end
else
    error("Unknown backend $backend_option: use --backend=cuda or --backend=oneapi.")
end

GPU.functional || @warn "$(GPU.name) is not functional on this machine: GPU benchmarks are skipped."
GPU.functional && !GPU.float64 && @info "$(GPU.description) has no Float64 support: GPU benchmarks in Float32 only."
