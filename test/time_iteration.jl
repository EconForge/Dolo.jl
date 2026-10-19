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

        @testset "$fn" begin

            model = include("$(root_dir)/examples/ymodels/$fn.jl")
            @test isbits(model)
            sol = Dolo.time_iteration(model; verbose=false)
            @test Dolo.converged(sol)
            dmodel = Dolo.discretize(model)
            @test max_residual(dmodel, sol.dr) < 1e-5

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

@testset "Solver variants" verbose=true begin

    # all variants should find the same solution as plain time iteration
    variants = [
        ("accelerate (neumann)", m->Dolo.time_iteration(m; verbose=false, accelerate=true)),
        ("accelerate (gmres)", m->Dolo.time_iteration(m; verbose=false, accelerate=true, linsolve=:gmres)),
        ("accelerate (safeguarded)", m->Dolo.time_iteration(m; verbose=false, accelerate=true, accelerate_wait=0, safeguarded=true)),
        ("accelerate (threaded)", m->Dolo.time_iteration(m; verbose=false, accelerate=true, engine=:threads)),
        ("newton (neumann)", m->Dolo.newton(m; verbose=false)),
        ("newton (gmres)", m->Dolo.newton(m; verbose=false, linsolve=:gmres)),
        ("newton (gmres, no preconditioning)", m->Dolo.newton(m; verbose=false, linsolve=:gmres, precondition=false, krylov_dim=300)),
        ("newton (backtracking, threaded)", m->Dolo.newton(m; verbose=false, globalization=:backtracking, engine=:threads)),
        ("newton (linesearch)", m->Dolo.newton(m; verbose=false, globalization=:linesearch)),
        ("accelerate (KernelAbstractions)", m->Dolo.time_iteration(m; verbose=false, accelerate=true, engine=:ka)),
    ]

    for fn in ["rbc_mc", "rbc_ar1", "consumption_savings_iid"]

        @testset "$fn" begin

            model = include("$(root_dir)/examples/ymodels/$fn.jl")
            dmodel = Dolo.discretize(model)
            on_grid(dr) = [(dmodel.grid isa Dolo.CGrid ? dr(s.val) : dr(s)) for s in Dolo.enum(dmodel.grid)]

            ref = Dolo.time_iteration(model; verbose=false, tol_η=1e-10)
            x_ref = on_grid(ref.dr)

            for (name, solve) in variants
                sol = solve(model)
                @test Dolo.converged(sol)
                @test sol.ε < sol.tol_ε
                @test maximum(maximum(abs.(a-b)) for (a,b) in zip(on_grid(sol.dr), x_ref)) < 1e-6
                if sol isa Dolo.TimeIterationResult
                    @test sol.iterations < ref.iterations/3
                end
            end

            # timings are only recorded on demand
            @test Dolo.time_iteration(model; verbose=false).timer === nothing
            @test Dolo.time_iteration(model; verbose=false, accelerate=true, timer=true).timer isa Dolo.TimerOutput

            # deprecated keywords
            sol_old = Dolo.time_iteration(model; verbose=false, improve=true, improve_K=1000)
            sol_new = Dolo.time_iteration(model; verbose=false, accelerate=true, accelerate_K=1000)
            @test sol_old.iterations == sol_new.iterations

        end

    end

    @testset "Invalid options" begin
        model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")
        dmodel = Dolo.discretize(model)
        @test_throws ErrorException Dolo.newton(dmodel; verbose=false, precondition=false)
        @test_throws ErrorException Dolo.time_iteration(dmodel; verbose=false, accelerate=true, linsolve=:lu)
        wk = Dolo.time_iteration_workspace(dmodel)
        @test_throws ErrorException Dolo.time_iteration(dmodel, wk; verbose=false, accelerate=true)
    end

end

@testset "Robustness" begin

    # interpolation at NaN gives NaN (and does not read out of bounds)
    model = include("$(root_dir)/examples/ymodels/consumption_savings_iid.jl")
    dm = Dolo.discretize(model)
    for interp_mode in (:linear, :cubic)
        (; x0, φ) = Dolo.time_iteration_workspace(dm; interp_mode)
        Dolo.fit!(φ, x0)
        @test all(isnan, φ(Dolo.SVector(NaN)))
    end

    # a diverging Newton stops at its last point, without error
    dm = Dolo.discretize(model, Dict(:endo => [10000]))
    sol = Dolo.newton(dm, Dolo.newton_workspace(dm; interp_mode=:cubic); verbose=false, globalization=:backtracking)
    @test !Dolo.converged(sol)
    @test isfinite(sol.ε)

    # `in` is only extended for Dolo's spaces
    @test which(in, Tuple{Dolo.SVector{1,Float64}, Vector{Any}}).module === Base

end
