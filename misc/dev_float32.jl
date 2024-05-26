using Dolo

root_dir = pkgdir(Dolo)
model_32 = include("$(root_dir)/misc/rbc_float32.jl")


model32 = Dolo.convert_precision(Float32, model_32)

dm32 = Dolo.discretize(model32, Dict(:endo=>[100]) )


typeof(dm32)

using Adapt


wk0 = Dolo.time_iteration_workspace(dm32)

import oneAPI: oneArray
import Adapt: adapt_structure

wk = adapt(oneArray, wk0)
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

@kernel function ggg(r, dm, x, φ)
# @kernel function fun(dm, r, x)


    n = @index(Global, Linear)
    ind = @index(Global, Cartesian)
    (i,j) = ind.I


    # k = length(dm.grid.grids[1])
    # i_ = (n-1)÷k
    # j_ = (n-1) - (i)*κ
    # i = i_+1
    # j = j_+1
    
    # TODO: special function here
    s_ = dm.grid[n]
    s = Dolo.QP((i,j), s_)

    xx = x[n]

    # (i,j) = @inline Dolo.from_linear(model.grid, n)    

    rr = Dolo.F(dm, s, xx, φ)

    r[n] = rr

end



fun_ = ggg(t_e)
@time fun_(wk.r0, dm32, wk.x0, wk.φ; ndrange=(p,q))



@time Dolo.F!(wk0.r0, dm32, wk0.x0, wk0.φ)

# try

    # res = fun_(wk.r0, dm32, wk.x0, wk.φ; ndrange=(p,q))
# catch err
    # nothing
# end



# Dolo.time_iteration(dm32, wk; engine=:gpu)


s_0 = dm32.grid[10]
s0 = Dolo.QP((2,2), s_0)
xx0 = wk0.x0[10]


fon(dm32, s0, xx0) = sum(w*S.val for (w,S) in Dolo.τ(dm32, s0, xx0))

@code_warntype fon(dm32, s0, xx0)