module DolooneAPIExt

# Loaded automatically when oneAPI is loaded together with Dolo (see DoloCUDAExt).
# Data is moved to the GPU with `adapt`, e.g. `time_iteration_workspace(dmodel; dest=oneArray)`,
# and solvers run on it with the KernelAbstractions engine (`engine=:ka`).
# Nothing is specific to oneAPI so far.

end
