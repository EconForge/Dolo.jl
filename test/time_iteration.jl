root_dir = pkgdir(Dolo)

# maximum absolute residual of the arbitrage equations on the grid
function max_residual(dmodel, dr)
    x = Dolo.GArray(
        dmodel.grid,
        [(dmodel.grid isa Dolo.CGrid ? dr(s.val) : dr(s)) for s in Dolo.enum(dmodel.grid)]
    )
    r = Dolo.F(dmodel, x, dr)
    return maximum(u->maximum(abs.(u)), r.data)
end

@testset "Benchmark Models" verbose=true begin

    model_files = ["rbc_iid", "rbc_mc", "rbc_ar1", "consumption_savings_iid"]

    for fn in model_files

        @testset "$fn" verbose=true begin

            @testset "Julia" verbose=true begin

                model = include("$(root_dir)/examples/ymodels/$fn.jl")
                @test isbits(model)
                sol = Dolo.time_iteration(model; verbose=false)
                @test Dolo.converged(sol)
                dmodel = Dolo.discretize(model)
                @test max_residual(dmodel, sol.dr) < 1e-5

            end

            @testset "YAML" verbose=true begin
                # YAML import is disabled until DoloYAML produces YModel objects
                @test_broken false
            end

        end

    end

end

@testset "Steady state (rbc_iid)" begin

    model = include("$(root_dir)/examples/ymodels/rbc_iid.jl")
    sol = Dolo.time_iteration(model; verbose=false)
    c = model.calibration
    x = sol.dr(Dolo.SVector(c.z, c.k))
    # stochastic and deterministic steady states are close
    @test isapprox(x[1], c.i; rtol=1e-2)
    @test isapprox(x[2], c.n; rtol=1e-2)

end
