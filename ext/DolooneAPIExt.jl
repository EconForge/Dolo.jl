module DolooneAPIExt

# Loaded automatically when oneAPI is loaded together with Dolo (see DoloCUDAExt).
# Data is moved to the GPU with `adapt`, e.g. `time_iteration_workspace(dmodel; dest=oneArray)`,
# and solvers run on it with the KernelAbstractions engine (`engine=:ka`).
# This only describes the device for `Dolo.diagnose`.

using oneAPI
import Dolo

const L0 = oneAPI.oneL0

function gpu_info()
    oneAPI.functional() || return (; functional=false, reason="`oneAPI.functional()` is false (no supported Intel GPU or driver)")
    dev = oneAPI.device()
    props = L0.properties(dev)
    integrated = props.flags & L0.ZE_DEVICE_PROPERTY_FLAG_INTEGRATED == L0.ZE_DEVICE_PROPERTY_FLAG_INTEGRATED
    fp64 = L0.module_properties(dev).fp64flags & L0.ZE_DEVICE_MODULE_FLAG_FP64 == L0.ZE_DEVICE_MODULE_FLAG_FP64
    return (;
        functional = true,
        device = props.name,
        memory_total = sum(m -> m.totalSize, L0.memory_properties(dev); init=0),
        # an integrated GPU uses the memory of the system
        memory_free = integrated ? Int(Sys.free_memory()) : nothing,
        integrated,
        float64 = fp64,
        float64_ratio = nothing,
        details = integrated ? "integrated" : "discrete",
        param_limit = nothing,   # not reported by Level Zero
    )
end

__init__() = Dolo.register_gpu!("oneAPI", gpu_info)

end
