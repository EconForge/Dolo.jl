module DoloCUDAExt

# Loaded automatically when CUDA is loaded together with Dolo.
# Data is moved to the GPU with `adapt`, e.g. `time_iteration_workspace(dmodel; dest=CuArray)`,
# and solvers run on it with the KernelAbstractions engine (`engine=:ka`).
# Nothing is specific to CUDA so far: the GPU code is written for any KernelAbstractions backend.

end
