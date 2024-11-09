using Dolo
using DoloYAML

model = include("examples/ymodels/consumption_savings.jl")
# model = include("examples/ymodels/rbc_iid.jl")

# model = DoloYAML.yaml_import("examples/ymodels/rbc_iid.yaml")

dmodel = Dolo.discretize(model)

wk = Dolo.time_iteration_workspace(dmodel; interp_mode=:linear, improve=true)

wk = Dolo.time_iteration_workspace(dmodel; interp_mode=:cubic, improve=false)
@time sol = Dolo.time_iteration(dmodel, wk; improve=false, verbose=false)

wk = Dolo.time_iteration_workspace(dmodel; interp_mode=:cubic, improve=true)
@time sol = Dolo.time_iteration(dmodel, wk; improve=true, verbose=false, improve_K=1000, improve_wait=1)

wk = Dolo.time_iteration_workspace(dmodel; interp_mode=:linear, improve=true)
@time sol = Dolo.time_iteration(model; improve=true, verbose=false)
@code_warntype sol = Dolo.time_iteration(model; improve=true, verbose=false)
@profview sol = Dolo.time_iteration(model; improve=true, verbose=false)

J = Dolo.dF_1(dmodel, wk.x0, wk.φ)
dx = deepcopy(wk.x0)
r0 = wk.r0

@time dx.data .= J.data .\ r0.data


Dolo.time_iteration(dmodel, wk; improve=true, improve_wait=5, T=50)


(J_1, J_2, d) = Dolo.time_iteration(dmodel, wk; improve=true, improve_wait=5)
T = J_1 \ J_2
M = convert(Matrix, T)
abs.(eigvals(M))



L = Dolo.dF_2(dmodel, sol.dr.values, sol.dr )
M = convert(Matrix, L)

##

dmodel = Dolo.discretize(model)

wk = Dolo.newton_workspace(dmodel)

res = Dolo.newton(dmodel, wk; verbose=true, K=10)
dx, T, r0 = res

# Dolo.neumann(T, r0)
M = convert(Matrix, T)
using LinearAlgebra

eigvals(M)

wk = Dolo.time_iteration_workspace(dmodel; improve=true)

using FiniteDiff

function residual_1(dmodel,wk,v)
    x = Dolo.unravel(wk.x0, v)
    r = Dolo.F(dmodel,x, wk.φ)
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

residual_1(dmodel, wk, v)

Jdiff_1 = FiniteDiff.finite_difference_jacobian(u->residual_1(dmodel, wk, u), v)

J_1 = Dolo.dF_1(dmodel, wk.x0, wk.φ)
Jmat_1 = convert(Matrix, J_1)

@assert maximum( abs.(Jdiff_1 - Jmat_1)./(1 .+ abs.(Jdiff_1)) ) <1e-6


Jdiff_2 = FiniteDiff.finite_difference_jacobian(u->residual_2(dmodel, wk, u), v)
J_2 = Dolo.dF_2(dmodel, wk.x0, wk.φ)

Jmat_2 = convert(Matrix, J_2)

using LinearMaps

L_2 = convert(LinearMap, J_2)

using Plots
plot()

spy(1:5000, 1:5000, Jmat)


spy(Jmat, xlims=(1, 6000), ylims=(1,6000))

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
