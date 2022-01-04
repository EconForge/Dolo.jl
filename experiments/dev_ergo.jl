using Dolo


# model = yaml_import("examples/models/rbc_mc.yaml")
model = yaml_import("examples/models/rbc.yaml")
# model = yaml_import("examples/models/rbc_iid.yaml")

sol = Dolo.improved_time_iteration(model)

G = Dolo.distG(model, sol)

Π, dΠ = Dolo.transition_matrix(G; diff=true)


DΠ = Π .!= 0.0
DdΠ = [ maximum(abs,DΠ[c]) for c in CartesianIndices(dΠ)]

Π = Dolo.transition_matrix(G)
sum(Π; dims=2) # should be only ones

x0 = G.x0

# G contains all informations to compute the ergodic_distribution
m = Dolo.ergodic_distribution(G)

Dolo.ergodic_distribution(model, sol)


# we can compute distribution transitions

@time μ1 = G(G.μ0, G.x0);
@time μ1 = G(G.μ0, G.x0);


x0 = G.x0
μ0 = G.μ0

# we compute differential operators w.r.t. arguments
μ1, ∂G_∂μ, ∂G_∂x = G(μ0, x0, diff=true);

# we can pass vectors as arguments
x0_flat = cat(G.x0.data...; dims=1)
μ1, ∂G_∂μ, ∂G_∂x = G(μ0, x0_flat; diff=true)

# differential operate *like* matrices and can be composed like usual matrices
∂G_∂μ(μ0*0.01)
∂G_∂μ*(μ0*0.01) # equivalent

∂G_∂x*(x0_flat*0.01)


[∂G_∂μ ∂G_∂μ]
[∂G_∂μ ∂G_∂x]([μ0; x0_flat]) # yes, it is kind of beautiful!



using FiniteDiff

Jμ_num = FiniteDiff.finite_difference_jacobian(mu->G(mu, x0_flat), μ0)
Jμ_exact = convert(Matrix, ∂G_∂μ)
maximum(abs, Jμ_num-Jμ_exact)
@assert maximum(abs, Jμ_num-Jμ_exact)<1e-9  # This works

@time Jx_num = FiniteDiff.finite_difference_jacobian(x->G(μ0, x), x0_flat)
@time Jx_exact = convert(Matrix,∂G_∂x)
maximum(abs, Jx_num-Jx_exact)

@assert maximum(abs, Jx_num-Jx_exact)<1e-9


using Plots


pl1 = spy(abs.(Jx_num).>1e-8, title="Numerical")
pl2 = spy(abs.(Jx_exact).>1e-8, title="Exact")
pl3 = spy(abs.(Jx_exact - Jx_num).>1e-8, title="Diff")
plot(pl1,pl2,pl3)





# M = reshape(Jx_exact, 5, 10, 5, 20)
# MM = reshape(permutedims(M, [1, 2, 4, 3]), 50, 100)

pl1 = spy(abs.(Jx_num).>1e-8, title="Numerical")
pl2 = spy(abs.(Jx_exact).>1e-8, title="Exact")
pl3 = spy(abs.(Jx_exact - Jx_num).>1e-8, title="Diff")
plot(pl1,pl2,pl3)
