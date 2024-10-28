using Dolo
using DoloYAML

# model = include("examples/ymodels/consumption_savings.jl")
model = include("examples/ymodels/rbc_iid.yaml")
Dolo.time_iteration(model; improve=false)

dmodel = Dolo.discretize(model)

wk = Dolo.time_iteration_workspace(dmodel; improve=true)

using FiniteDiff

function residual_1(dmodel,wk,v)
    x = Dolo.unravel(wk.x0, v)
    r = Dolo.F(dmodel,x, wk.ψ)
    Dolo.ravel(r)
end

ψ = deepcopy(wk.φ)

function residual_2(dmodel,wk,v)
    x = Dolo.unravel(wk.x0, v)
    Dolo.fit!(ψ, x)
    r = Dolo.F(dmodel,wk.x0, ψ)
    Dolo.ravel(r)
end



v = Dolo.ravel(wk.x0)

residual(dmodel, wk, v, wk.φ)

Jdiff = FiniteDiff.finite_difference_jacobian(u->residual(dmodel, wk, u, ψ), v)

J = Dolo.dF_1(dmodel, wk.x0, wk.φ)
Jmat = convert(Matrix, J)

@assert maximum( abs.(Jdiff - Jmat)./(1 .+ abs.(Jdiff)) ) <1e-6


Jdiff_2 = FiniteDiff.finite_difference_jacobian(u->residual_2(dmodel, wk, u), v)
J_2 = Dolo.dF_2(dmodel, wk.x0, wk.φ)

Jmat_2 = convert(Matrix, J_2)

using LinearMaps

L_2 = convert(LinearMap, J_2)

using Plots
plot(xlims=(-10, 10), ylims=(-2,2))

spy(Jmat)
spy(Jmat_2, xlims=(1,150), ylims=(1,150))

spy(Jmat_2)


Δ = Jdiff_2 - Jmat_2
maximum(abs.(Δ)./( 1 .+ Jmat_2))
# # model = DoloYAML.yaml_import("examples/ymodels/consumption_savings.yaml")


# wk = Dolo.time_iteration_workspace(dmodel; improve=true)

# # @time Dolo.time_iteration(dmodel, wk; verbose=false, improve=false, improve_wait=10);

# @time Dolo.time_iteration(dmodel, wk; 
#     verbose=false, improve=false, improve_wait=10, engine=nothing);


# dm.grid[1]

# [dm.grid...]


# # wk = Dolo.time_iteration_workspace(model; improve=true, interp_mode=:cubic)
# # @time Dolo.time_iteration(model, wk; verbose=true, improve=true, improve_wait=20);



# # (r0, model, x0, φ, t_engine) = Dolo.time_iteration(model; verbose=false, improve=true);

# # Dolo.F!(r0, model, x0, φ, t_engine)
