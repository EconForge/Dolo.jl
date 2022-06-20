using FiniteDiff
using Dolo

@testset "Test the G derivative w.r.t. x, μ, z1 and z2" begin

    model_rbc_mc = yaml_import("examples/models/rbc_mc.yaml")
    model_rbc = yaml_import("examples/models/rbc.yaml")
    model_rbc_iid = yaml_import("examples/models/rbc_iid.yaml")

    for model in [model_rbc_mc, model_rbc, model_rbc_iid]

        sol = Dolo.improved_time_iteration(model)
        G = Dolo.distG(model, sol)
        z10 = SVector(model.calibration[:exogenous]...) 
        z20 = z10 
        x0 = G.x0
        x0_flat = cat(G.x0.data...; dims=1)
        μ0 = G.μ0



        μ1, ∂G_∂μ, ∂G_∂x, ∂G_∂z1, ∂G_∂z2 = G(μ0, x0; exo = [z10,z20], diff = true)

        Jμ_exact = convert(Matrix, ∂G_∂μ)
        Jμ_num = FiniteDiff.finite_difference_jacobian(mu -> G(mu, x0, exo = [z10,z20]), μ0)

        Jx_exact = convert(Matrix, ∂G_∂x)
        Jx_num = FiniteDiff.finite_difference_jacobian(x -> G(μ0, x; exo = [z10,z20]), x0)

        Jz1_exact = convert(Matrix, ∂G_∂z1)
        Jz1_num = FiniteDiff.finite_difference_jacobian(z1 -> G(μ0, x0; exo = [z1,z20]), z10)

        Jz2_exact = convert(Matrix, ∂G_∂z2)
        Jz2_num = FiniteDiff.finite_difference_jacobian(z2 -> G(μ0, x0; exo = [z10,z2]), z20)



        @assert maximum(abs, Jμ_num - Jμ_exact) < 1e-8

        @assert maximum(abs, Jx_num - Jx_exact) < 1e-8

        @assert maximum(abs, Jz1_num - Jz1_exact) < 1e-8

        @assert maximum(abs, Jz2_num - Jz2_exact) < 1e-5


    end

end