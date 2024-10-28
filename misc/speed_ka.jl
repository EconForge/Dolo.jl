using Dolo

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")

dm = Dolo.discretize(model, Dict(:endo=>[10000000]) )


wk = Dolo.time_iteration_workspace(dm)
(;x0, φ, r0) = wk

r1 = deepcopy(r0)

# s = Dolo.QP((1,1), dm.grid[1,1])
# x = x0[1]

using LoopVectorization

function residual1!(r0, dm,x0,φ)
    
    N = length(dm.grid)

    function fun(n)
        (i,j) = Dolo.from_linear(dm.grid, n)
        s_ = dm.grid[i,j]
        s = Dolo.QP((i,j), s_)
        x = x0.data[n]
        r = Dolo.F(dm, s, x, φ)
        r
    #     r0[i,j] = r
    end

    vmap!(fun,r0.data,1:N)

end

@time residual1!(r0, dm, x0, φ);

@benchmark residual1!(r0, dm, x0, φ)

using BenchmarkTools

@benchmark residual!(r1, dm, x0, φ)


res = @code_native Dolo.F!(r0, dm, x0, φ)

