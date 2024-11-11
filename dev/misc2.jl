using Dolo

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")

@time res = Dolo.time_iteration(model; verbose=false, discretization=Dict(:endo=>[100000]));

@time res = Dolo.time_iteration(model; verbose=false, engine=:cpu, discretization=Dict(:endo=>[100000]));

###

dm = Dolo.discretize(model)


discr_options = Dict(:endo=>[100])

dm2 = Dolo.discretize(model, discr_options)


###

dm = Dolo.discretize(model, Dict(:endo=>[1000000]))
wk = Dolo.time_iteration_workspace(dm)

(;x0, x1, J, φ) = wk

r0 = Dolo.F(dm, x0, φ)*0
r1 = Dolo.F(dm, x0, φ)*0

Dolo.F!(r0, dm, x0, φ)

using Dolo: GArray, DFun, enum

import Base: Threads

function FF!(rv, model, controls::GArray, φ::Union{GArray, DFun})

    Threads.@threads for n = 1:length(model.grid)
        
        _s = model.grid[n]

        ind = Dolo.from_linear(model.grid, n)

        s = Dolo.QP(ind, _s)
        x = controls[n]

        r = Dolo.F(model, s, x, φ)
        
        rv[n] = r

    end
end

using Dolo: CPU


@time Dolo.F!(r0, dm, x0, φ)
@time Dolo.F!(r0, dm, x0, φ, CPU())
@time FF!(r1, dm, x0, φ)


using SIMD

v0 = GArray(
    dm.grid,
    [Vec(0.1, 0.2, 0.4, 0.2) for i=1:length(dm.grid)]
)

v1 = GArray(
    dm.grid,
    [0.1 for i=1:length(dm.grid)*4]
)

@time v0+v0;
@time v1+v1;