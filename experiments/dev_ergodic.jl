using Dolo 
using BenchmarkTools
using StaticArrays
using FiniteDiff


model_rbc_mc = yaml_import("examples/models/rbc_mc.yaml")
model = model_rbc_mc
sol = Dolo.improved_time_iteration(model)

G = Dolo.distG(model, sol)
z10 = SVector{2,Float64}( model.exogenous.values[1,:])
z20 = SVector{2,Float64}(model.exogenous.values[2,:])
x0 = G.x0
μ0 = G.μ0

μ1, ∂G_∂μ, ∂G_∂x, ∂G_∂z1, ∂G_∂z2 = G(μ0, x0, exo = [z10,z20], diff = true)

# Jμ_exact = convert(Matrix, ∂G_∂μ)
# Jμ_num = FiniteDiff.finite_difference_jacobian(mu -> G(mu, x0, exo = [z10,z20]), μ0)

# Jx_exact = convert(Matrix, ∂G_∂x)
# Jx_num = FiniteDiff.finite_difference_jacobian(x -> G(μ0, x, exo = [z10,z20]), x0)

Jz1_exact = convert(Matrix, ∂G_∂z1)
Jz1_num = FiniteDiff.finite_difference_jacobian(z1 -> G(μ0, x0; exo = [z1,z20]), z10)

Jz2_exact = convert(Matrix, ∂G_∂z2)
Jz2_num = FiniteDiff.finite_difference_jacobian(z2 -> G(μ0, x0, exo = [z10,z2]), z20)

# print( maximum(abs, Jμ_num - Jμ_exact) < 1e-8)

# print(maximum(abs, Jx_num - Jx_exact) < 1e-8)

print(maximum(abs, Jz1_num - Jz1_exact) < 1e-8)

print(maximum(abs, Jz2_num - Jz2_exact) < 1e-8)



maximum(Jz1_num*100000000)

