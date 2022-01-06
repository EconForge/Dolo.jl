using FiniteDiff
using Dolo

@testset "Test the G derivative w.r.t. x" begin

    model_rbc_mc = yaml_import("examples/models/rbc_mc.yaml")
    model_rbc = yaml_import("examples/models/rbc.yaml")
    model_rbc_iid = yaml_import("examples/models/rbc_iid.yaml")
    for model in [model_rbc_mc, model_rbc, model_rbc_iid]

        sol = Dolo.improved_time_iteration(model)

        G = Dolo.distG(model, sol)

        x0_not_flat = G.x0
        x0_flat = cat(G.x0.data...; dims = 1)

        μ0 = G.μ0

        for x0 in [x0_not_flat, x0_flat]

            μ1, ∂G_∂μ, ∂G_∂x = G(μ0, x0, diff = true)

            Jμ_exact = convert(Matrix, ∂G_∂μ)
            Jμ_num = FiniteDiff.finite_difference_jacobian(mu -> G(mu, x0), μ0)

            Jx_exact = convert(Matrix, ∂G_∂x)
            Jx_num = FiniteDiff.finite_difference_jacobian(x -> G(μ0, x), x0)

            @assert maximum(abs, Jμ_num - Jμ_exact) < 1e-8

            @assert maximum(abs, Jx_num - Jx_exact) < 1e-8
        end
        
    end

end