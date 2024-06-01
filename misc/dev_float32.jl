using Dolo

root_dir = pkgdir(Dolo)
model_32 = include("$(root_dir)/misc/rbc_float32.jl")


model32 = Dolo.convert_precision(Float32, model_32)

dm32 = Dolo.discretize(model32, Dict(:endo=>[1000]) )

typeof(dm32)


using Adapt
import oneAPI: oneArray
import Adapt: adapt_structure
# import CUDA: CuArray

import Dolo

wk0 = Dolo.time_iteration_workspace(dm32);

# wk = Dolo.time_iteration_workspace(dm32, dest=CuArray)
@time wk = Dolo.time_iteration_workspace(dm32, dest=oneArray; interp_mode=:cubic);


time_iteration(dm32, wk0; engine=:gpu, tol_η=1e-5, improve_wait=3, improve=true)
time_iteration(dm32, wk; engine=:gpu, tol_η=1e-5, improve_wait=3, improve=true)





# list available options

# @time wk = Dolo.time_iteration_workspace(dm32, dest=Array; interp_mode=:cubic);

# @time Dolo.splines.prefilter!(c0, Val(:CPU))
# @time Dolo.splines.prefilter!(c0, Val(:Threads))
# @time Dolo.splines.prefilter!(c0, Val(:KA))

# @time Dolo.splines.prefilter!(c0)
# @time Dolo.splines.prefilter!(c1,Val(:KA))

# @time Dolo.splines.prefilter!(c,Val(:KA))

# c_gpu = adapt(Array, c)

# # Dolo.splines.prefilter!(c)



# t_e = get_backend(wk.x0)

# using KernelAbstractions: CPU
# using StaticArrays


# L0 = Dolo.dF_2(dm32, wk0.x0*0, wk0.φ)
# L0_gpu = adapt(oneArray, L0)

# L = deepcopy(L0)

# @time Dolo.dF_2!(L, dm32, wk0.x0, wk0.φ, CPU())

# # @time Dolo.dd_2!(L, dm32, wk0.x0, wk0.φ, CPU())


# @time Dolo.dF_2!(L0_gpu, dm32, wk.x0, wk.φ, t_e)

# Lg = adapt(Array, L0_gpu)


# @time L*wk0.x0;

# #It works !
# using BenchmarkTools
# @benchmark Dolo.mul!(wk0.r0, L, wk0.x0, CPU())
# @benchmark Dolo.mul!(wk.r0, L0_gpu, wk.x0, t_e)
# # try

#     # res = fun_(wk.r0, dm32, wk.x0, wk.φ; ndrange=(p,q))
# # catch err
#     # nothing
# # end

# wk0 = Dolo.time_iteration_workspace(dm32; interp_mode=:cubic)
# # wk = adapt(CuArray, wk0)
# wk = adapt(oneArray, wk0)

# @time Dolo.time_iteration(dm32, wk, verbose=true, improve=true; engine=:gpu, improve_wait=10, improve_K=500, T=20);




# @time Dolo.time_iteration(dm32, wk0, verbose=true, improve=false;  T=20);




# using ForwardDiff

# ForwardDiff.derivative(u->unsafe_trunc(u), 0.1)

# ForwardDiff.jacobian(u->trunc.(Int, u), SVector(0.3, 0.2))

# ForwardDiff.jacobian(u->unsafe_trunc.(Int, u), SVector(0.3, 0.2))
