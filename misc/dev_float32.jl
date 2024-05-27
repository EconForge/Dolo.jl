using Dolo

root_dir = pkgdir(Dolo)
model_32 = include("$(root_dir)/misc/rbc_float32.jl")


model32 = Dolo.convert_precision(Float32, model_32)

dm32 = Dolo.discretize(model32, Dict(:endo=>[100000]) )


typeof(dm32)

using Adapt

import oneAPI: oneArray
import Adapt: adapt_structure
import CUDA: CuArray


wk0 = Dolo.time_iteration_workspace(dm32)

wk = Dolo.time_iteration_workspace(dm32, dest=CuArray)

# wk = adapt(oneArray, wk0)

wk = adapt(CuArray, wk0)
using KernelAbstractions: get_backend

if typeof(dm32.grid)<:Dolo.CGrid
    p = dm32.grid.ranges[1][3]
    q = dm32.grid.ranges[2][3]
else
    p = length(dm32.grid.grids[1])
    q = length(dm32.grid.grids[2])
end

t_e = get_backend(wk.x0)

using KernelAbstractions
using StaticArrays


# try

    # res = fun_(wk.r0, dm32, wk.x0, wk.φ; ndrange=(p,q))
# catch err
    # nothing
# end

wk0 = Dolo.time_iteration_workspace(dm32)
wk = adapt(CuArray, wk0)
@time Dolo.time_iteration(dm32, wk, verbose=false; engine=:gpu);

@time Dolo.time_iteration(dm32, wk0, verbose=false; engine=:gpu);


Dolo.distance(adapt(Array,wk.x0) , wk0.x0)

# 357 times faster !