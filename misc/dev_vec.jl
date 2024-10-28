using Dolo
using Dolo: transition

# can I vectorize evaluation of F?

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/consumption_savings.jl")

dm = Dolo.discretize(model)

wk = Dolo.time_iteration_workspace(dm)


x0 = wk.x0
φ = wk.φ
r0 = Dolo.F(dm, x0, φ)


t_engine = Dolo.get_backend(wk.x0)


using BenchmarkTools

@time Dolo.F!(r0,dm, x0, φ)
@time Dolo.F!(r0,dm, x0, φ, t_engine)



@benchmark Dolo.F!(r0,dm, x0, φ)
@benchmark Dolo.F!(r0,dm, x0, φ, t_engine)

# the second is two times slower and makes a few allocations
