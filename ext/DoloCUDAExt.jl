module DoloCUDAExt

# Loaded automatically when CUDA is loaded together with Dolo.
# Data is moved to the GPU with `adapt`, e.g. `time_iteration_workspace(dmodel; dest=CuArray)`,
# and solvers run on it with the KernelAbstractions engine (`engine=:ka`).

using CUDA: CuArray
import Dolo: GVector, distance

# the generic method iterates with `zip`, which would require scalar indexing
distance(x::GVector{G,<:CuArray}, y::GVector{G,<:CuArray}) where G = mapreduce(distance, max, x.data, y.data)

end
