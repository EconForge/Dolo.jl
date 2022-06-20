using FiniteDiff
using Dolo
using StaticArrays

@testset "Test F derivatives w.r.t. x1, x2 and p" begin

    model_rbc_mc = yaml_import("examples/models/rbc_mc.yaml")
    model_rbc = yaml_import("examples/models/rbc.yaml")
    model_rbc_iid = yaml_import("examples/models/rbc_iid.yaml")

    for model in [model_rbc_mc, model_rbc, model_rbc_iid]
        F = Dolo.Euler(model)

        x0_not_flat = F.x0
        z0 = SVector(model.calibration[:exogenous]...)
        z1 = z0
        x0_flat = cat(G.x0.data...; dims = 1)

        n_z = length(model.symbols[:exogenous])
        n_x = length(x0_flat)


        dF_A = Dolo.df_A(F, x0_not_flat, x0_not_flat; exo=(z0,z1))
        dF_B = Dolo.df_B(F, x0_not_flat, x0_not_flat; exo=(z0,z1))
        dF_e = sum(Dolo.df_e(F, x0_not_flat, x0_not_flat, z0, z1))


        J_A_num = FiniteDiff.finite_difference_jacobian(x -> cat(F(x, x0_flat; exo=(z0,z1)).data...; dims = 1), x0_flat)
        J_B_num = FiniteDiff.finite_difference_jacobian(x -> cat(F(x0_flat, x, exo=(z0,z1)).data...; dims = 1), x0_flat)
        J_e_num = FiniteDiff.finite_difference_jacobian(e -> cat(F(x0_not_flat, x0_not_flat, e, e).data...; dims = 1), z0)

            
        dx0 =  rand(n_x) ./ 1000
        dz0 = rand(n_z) ./ 1000

        @assert maximum(abs.(J_A_num * dx0 - dF_A * dx0)) < 1e-7
        @assert maximum(abs.(J_B_num * dx0 - dF_B * dx0)) < 1e-7
        @assert maximum(abs.(J_e_num * dz0 - dF_e * dz0)) < 1e-7


    end

end
