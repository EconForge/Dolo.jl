# Tests of the GPU code paths, for any array type `AT` with a KernelAbstractions backend:
# CuArray (test/cuda), oneArray (test/oneapi), or JLArray (a reference implementation of a
# GPU array, on the CPU, run with the main test suite).
#
# Results on `AT` are compared with the CPU, in each of the given precisions. Scalar indexing
# of `AT` must be disallowed by the caller (e.g. `CUDA.allowscalar(false)`): it would be a bug.

using Dolo
using Test

# largest difference between two vectors of SVector/SMatrix, relative to the size of the second
function maxdiff(a, b)
    a, b = Array(a), Array(b)
    d = maximum(maximum(abs.(u-v)) for (u,v) in zip(a, b))
    return d / max(1, maximum(maximum(abs.(v)) for v in b))
end
maxdiff(a::Dolo.LL, b::Dolo.LL) = maxdiff(
    [u.F_x for A in Array(a.D) for u in A],
    [v.F_x for B in Array(b.D) for v in B]
)

# tolerances, in each precision
tol(::Type{Float64}, x) = x
tol(::Type{Float32}, x) = max(x, 1e-3)

model_in(T, model) = T == Float64 ? model : Dolo.convert_precision(T, model)

"""
    gpu_tests(AT; precisions=(Float64,), solvers=true)

Run the GPU tests with arrays of type `AT` (e.g. `CuArray`), in each of the given
precisions. `solvers=false` only tests the operations, not the complete solution.
"""
function gpu_tests(AT; precisions=(Float64,), solvers=true)

    @testset "Spline prefilter ($T): size $n" for T in precisions, n in ((10,), (10_000,), (300, 200), (40, 30, 50))
        θ = zeros(Dolo.SVector{2,T}, (k+2 for k in n)...)
        θ[(2:k+1 for k in n)...] .= rand(Dolo.SVector{2,T}, n...)
        cpu = Dolo.splines.prefilter!(copy(θ))       # sequential
        gpu = Array(Dolo.splines.prefilter!(AT(θ)))  # chunked
        @test maxdiff(gpu, cpu) < tol(T, 1e-13)
    end

    for T in precisions, fn in ["rbc_mc", "rbc_ar1", "consumption_savings_iid"]
        model = include("$(pkgdir(Dolo))/examples/ymodels/$fn.jl")
        # the model's methods were just defined: `invokelatest` makes them visible
        Base.invokelatest(model_tests, AT, T, fn, model_in(T, model); solvers)
    end

    for T in precisions
        model = include("$(pkgdir(Dolo))/examples/ymodels/integration_B.jl")
        Base.invokelatest(integration_tests, AT, T, model_in(T, model); solvers)
    end

end


"Tests with model `model` (from example `fn`), in precision `T`."
function model_tests(AT, T, fn, model; solvers=true)

    @testset "$fn ($T)" verbose=true begin

        dm = Dolo.discretize(model)
        @test isbits(dm)

        @testset "Operations" begin

            cpu = Dolo.newton_workspace(dm)
            gpu = Dolo.newton_workspace(dm; dest=AT)
            @test gpu.x0.data isa AT
            @test eltype(eltype(gpu.x0.data)) == T

            engine = Dolo.select_engine(:ka, gpu.x0)
            @test engine isa Dolo.KA

            # non-constant controls (any interpolation is exact on a constant)
            cpu.x0.data .= [x .* T(1 + 0.1*sin(n)) for (n,x) in enumerate(cpu.x0.data)]
            copyto!(gpu.x0.data, cpu.x0.data)
            for wk in (cpu, gpu)
                Dolo.fit!(wk.φ, wk.x0)
            end

            Dolo.F!(cpu.r0, dm, cpu.x0, cpu.φ)
            Dolo.F!(gpu.r0, dm, gpu.x0, gpu.φ, engine)
            @test maxdiff(gpu.r0.data, cpu.r0.data) < tol(T, 1e-12)

            Dolo.dF_1!(cpu.J, dm, cpu.x0, cpu.φ)
            Dolo.dF_1!(gpu.J, dm, gpu.x0, gpu.φ, engine)
            @test maxdiff(gpu.J.data, cpu.J.data) < tol(T, 1e-10)

            Dolo.dF_2!(cpu.L, dm, cpu.x0, cpu.φ)
            Dolo.dF_2!(gpu.L, dm, gpu.x0, gpu.φ, engine)
            @test maxdiff(gpu.L, cpu.L) < tol(T, 1e-10)

            Dolo.mul!(cpu.r0, cpu.L, cpu.x0)
            Dolo.mul!(gpu.r0, gpu.L, gpu.x0, engine)
            @test maxdiff(gpu.r0.data, cpu.r0.data) < tol(T, 1e-10)

            Dolo.mul!(cpu.L, -one(T))
            Dolo.mul!(gpu.L, -one(T), engine)
            Dolo.ldiv!(cpu.J, cpu.L)
            Dolo.ldiv!(gpu.J, gpu.L, engine)
            @test maxdiff(gpu.L, cpu.L) < tol(T, 1e-10)

            @test isapprox(Dolo.norm(gpu.r0), Dolo.norm(cpu.r0); rtol=tol(T, 1e-12))
            @test isapprox(Dolo.distance(gpu.r0, gpu.x0), Dolo.distance(cpu.r0, cpu.x0); rtol=tol(T, 1e-12))

        end

        solvers && @testset "Solvers ($interp)" for interp in (:linear, :cubic)

            # tolerances reachable in each precision
            tol_ε = T == Float64 ? 1e-8 : 1e-4
            tol_η = T == Float64 ? 1e-10 : 1e-5
            linsolve_tol = T == Float64 ? 1e-10 : 1e-5

            # reference solution, on the CPU
            on_grid(dr) = [(dm.grid isa Dolo.CGrid ? dr(s.val) : dr(s)) for s in Dolo.enum(dm.grid)]
            ref = Dolo.time_iteration(dm, Dolo.newton_workspace(dm; interp_mode=interp);
                verbose=false, accelerate=true, tol_ε=T == Float64 ? 1e-10 : 1e-5, linsolve_tol)
            x_ref = on_grid(ref.dr)

            kw = (; verbose=false, engine=:ka, tol_ε, linsolve_tol)
            variants = [
                "time iteration" => wk->Dolo.time_iteration(dm, wk; kw..., tol_η),
                "accelerated (neumann)" => wk->Dolo.time_iteration(dm, wk; kw..., accelerate=true),
                "newton (backtracking)" => wk->Dolo.newton(dm, wk; kw..., globalization=:backtracking),
                "accelerated (gmres)" => wk->Dolo.time_iteration(dm, wk; kw..., accelerate=true, linsolve=:gmres),
                "newton (gmres)" => wk->Dolo.newton(dm, wk; kw..., globalization=:backtracking, linsolve=:gmres),
            ]
            @testset "$name" for (name, solve) in variants
                sol = solve(Dolo.newton_workspace(dm; interp_mode=interp, dest=AT, linsolve=:gmres))
                @test Dolo.converged(sol)
                @test sol.dr.values.data isa AT
                x = on_grid(Dolo.adapt(Array, sol.dr))
                # time iteration stops when steps are smaller than tol_η: the error is about
                # tol_η/(1-λ), which is ~1e-3 in Float32 for slowly converging models
                slow = name == "time iteration" && T == Float32
                @test maxdiff(x, x_ref) < tol(T, 1e-6) * (slow ? 10 : 1)
            end

        end

    end

end

"Tests with integration_B (8 controls, 9 exogenous nodes), in precision `T`."
function integration_tests(AT, T, model; solvers=true)

    @testset "integration_B ($T)" begin

        # 8 controls, 9 exogenous nodes, 3 endogenous states
        dm = Dolo.discretize(model, Dict(:endo => [5, 5, 5]))
        cpu = Dolo.newton_workspace(dm; interp_mode=:cubic)
        cpu.x0.data .= [x .* (1 .+ T(0.01)*sin.(T.(1:8) .* n)) for (n,x) in enumerate(cpu.x0.data)]
        gpu = Dolo.adapt(AT, cpu)
        @test gpu.L.D isa Dolo.BlockVector
        engine = Dolo.select_engine(:ka, gpu.x0)
        for wk in (cpu, gpu)
            Dolo.fit!(wk.φ, wk.x0)
        end
        @test maxdiff(Dolo.F!(gpu.r0, dm, gpu.x0, gpu.φ, engine).data, Dolo.F!(cpu.r0, dm, cpu.x0, cpu.φ).data) < tol(T, 1e-12)
        @test maxdiff(Dolo.dF_1!(gpu.J, dm, gpu.x0, gpu.φ, engine).data, Dolo.dF_1!(cpu.J, dm, cpu.x0, cpu.φ).data) < tol(T, 1e-10)
        Dolo.dF_2!(cpu.L, dm, cpu.x0, cpu.φ)
        Dolo.dF_2!(gpu.L, dm, gpu.x0, gpu.φ, engine)
        @test maxdiff(gpu.L, cpu.L) < tol(T, 1e-10)
        @test maxdiff(Dolo.mul!(gpu.r0, gpu.L, gpu.x0, engine).data, Dolo.mul!(cpu.r0, cpu.L, cpu.x0).data) < tol(T, 1e-10)
        @test maxdiff(Dolo.ldiv!(gpu.J, gpu.L, engine), Dolo.ldiv!(cpu.J, cpu.L)) < tol(T, 1e-10)

        # a few time iteration steps, then Newton: same solution as on the CPU
        if solvers
            tol_ε = T == Float64 ? 1e-8 : 1e-4
            linsolve_tol = T == Float64 ? 1e-10 : 1e-5
            function solve(dest, engine)
                wk = Dolo.newton_workspace(dm; interp_mode=:cubic, dest=dest, linsolve=:gmres)
                Dolo.time_iteration(dm, wk; T=25, verbose=false, engine=engine, tol_ε)
                Dolo.newton(dm, wk; verbose=false, engine=engine, linsolve=:gmres, tol_ε, linsolve_tol)
            end
            ref = solve(Array, :threads)
            sol = solve(AT, :ka)
            @test Dolo.converged(sol)
            @test maxdiff(sol.dr.values.data, ref.dr.values.data) < tol(T, 1e-8)
        end

    end

end
