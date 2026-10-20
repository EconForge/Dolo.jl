module DoloCUDAExt

# Loaded automatically when CUDA is loaded together with Dolo.
# Data is moved to the GPU with `adapt`, e.g. `time_iteration_workspace(dmodel; dest=CuArray)`,
# and solvers run on it with the KernelAbstractions engine (`engine=:ka`).
# The GPU code is written for any KernelAbstractions backend: this only describes the device
# for `Dolo.diagnose`.

using CUDA
import Dolo

# throughput of Float64 relative to Float32, from the compute capability (approximate)
function float64_ratio(cap)
    (cap.major, cap.minor) in ((6, 0), (7, 0), (8, 0), (9, 0), (10, 0)) && return 1//2  # data center GPUs (P100, V100, A100, H100, B200)
    cap.major >= 8 && return 1//64     # consumer Ampere, Ada, Blackwell
    return 1//32                       # consumer Pascal, Turing
end

function gpu_info()
    CUDA.functional() || return (; functional=false, reason="`CUDA.functional()` is false (no NVIDIA driver or GPU?)")
    dev = CUDA.device()
    free, total = CUDA.memory_info()
    cap = CUDA.capability(dev)
    return (;
        functional = true,
        device = CUDA.name(dev),
        memory_total = Int(total),
        memory_free = Int(free),
        integrated = false,
        float64 = true,
        float64_ratio = float64_ratio(cap),
        details = "compute capability $(cap.major).$(cap.minor), CUDA driver $(CUDA.driver_version())",
        # 4 KB of kernel parameters, 32 KB from Volta on with CUDA 12.1
        param_limit = cap >= v"7.0" && CUDA.driver_version() >= v"12.1" ? 32764 : 4096,
    )
end

__init__() = Dolo.register_gpu!("CUDA", gpu_info)

end
