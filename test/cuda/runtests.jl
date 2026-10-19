# GPU tests, run in their own environment (test/cuda/Project.toml), which loads CUDA
# and therefore the DoloCUDAExt extension:
#
#     julia --project=test/cuda -e 'using Pkg; Pkg.instantiate()'
#     julia --project=test/cuda test/cuda/runtests.jl
#
# or from the main test suite with `Pkg.test("Dolo"; test_args=["cuda"])`.

using CUDA
using Dolo
using Test

@testset verbose=true "CUDA" begin

    @test Base.get_extension(Dolo, :DoloCUDAExt) !== nothing

    if !CUDA.functional()
        @info "CUDA is not functional on this machine: skipping GPU tests."
    else

        CUDA.allowscalar(false)   # any scalar indexing on the GPU is a bug

        # largest difference between a CPU and a GPU vector of SVector/SMatrix
        maxdiff(a, b) = maximum(maximum(abs.(u-v)) for (u,v) in zip(Array(a), Array(b)))
        maxdiff(a::Dolo.LL, b::Dolo.LL) = maximum(
            maximum(abs.(u.F_x-v.F_x)) for (A,B) in zip(Array(a.D), Array(b.D)) for (u,v) in zip(A,B)
        )

        @testset "Spline prefilter: size $n" for n in ((10,), (10_000,), (300, 200), (40, 30, 50))
            θ = zeros(Dolo.SVector{2,Float64}, (k+2 for k in n)...)
            θ[(2:k+1 for k in n)...] .= rand(Dolo.SVector{2,Float64}, n...)
            cpu = Dolo.splines.prefilter!(copy(θ))           # sequential
            gpu = Array(Dolo.splines.prefilter!(CuArray(θ))) # chunked
            @test maximum(maximum(abs.(u-v)) for (u,v) in zip(cpu, gpu)) < 1e-13
        end

        for fn in ["rbc_mc", "rbc_ar1", "consumption_savings_iid"]

            @testset "$fn" verbose=true begin

                model = include("$(pkgdir(Dolo))/examples/ymodels/$fn.jl")
                dm = Dolo.discretize(model)
                @test isbits(dm)

                @testset "Operations" begin

                    cpu = Dolo.newton_workspace(dm)
                    gpu = Dolo.newton_workspace(dm; dest=CuArray)
                    @test gpu.x0.data isa CuArray

                    engine = Dolo.select_engine(:ka, gpu.x0)
                    @test engine isa Dolo.KA{CUDABackend}

                    # non-constant controls (any interpolation is exact on a constant)
                    cpu.x0.data .= [x .* (1 + 0.1*sin(n)) for (n,x) in enumerate(cpu.x0.data)]
                    copyto!(gpu.x0.data, cpu.x0.data)
                    for wk in (cpu, gpu)
                        Dolo.fit!(wk.φ, wk.x0)
                    end

                    Dolo.F!(cpu.r0, dm, cpu.x0, cpu.φ)
                    Dolo.F!(gpu.r0, dm, gpu.x0, gpu.φ, engine)
                    @test maxdiff(cpu.r0.data, gpu.r0.data) < 1e-12

                    Dolo.dF_1!(cpu.J, dm, cpu.x0, cpu.φ)
                    Dolo.dF_1!(gpu.J, dm, gpu.x0, gpu.φ, engine)
                    @test maxdiff(cpu.J.data, gpu.J.data) < 1e-10

                    Dolo.dF_2!(cpu.L, dm, cpu.x0, cpu.φ)
                    Dolo.dF_2!(gpu.L, dm, gpu.x0, gpu.φ, engine)
                    @test maxdiff(cpu.L, gpu.L) < 1e-10

                    Dolo.mul!(cpu.r0, cpu.L, cpu.x0)
                    Dolo.mul!(gpu.r0, gpu.L, gpu.x0, engine)
                    @test maxdiff(cpu.r0.data, gpu.r0.data) < 1e-10

                    Dolo.mul!(cpu.L, -1.0)
                    Dolo.mul!(gpu.L, -1.0, engine)
                    Dolo.ldiv!(cpu.J, cpu.L)
                    Dolo.ldiv!(gpu.J, gpu.L, engine)
                    @test maxdiff(cpu.L, gpu.L) < 1e-10

                    @test Dolo.norm(gpu.r0) ≈ Dolo.norm(cpu.r0)
                    @test Dolo.distance(gpu.r0, gpu.x0) ≈ Dolo.distance(cpu.r0, cpu.x0)

                end

                @testset "Solvers ($interp)" for interp in (:linear, :cubic)

                    # reference solution, on the CPU
                    on_grid(dr) = [(dm.grid isa Dolo.CGrid ? dr(s.val) : dr(s)) for s in Dolo.enum(dm.grid)]
                    ref = Dolo.time_iteration(dm, Dolo.newton_workspace(dm; interp_mode=interp);
                        verbose=false, accelerate=true, tol_ε=1e-10)
                    x_ref = on_grid(ref.dr)

                    variants = [
                        "time iteration" => wk->Dolo.time_iteration(dm, wk; verbose=false, engine=:ka, tol_η=1e-10),
                        "accelerated (neumann)" => wk->Dolo.time_iteration(dm, wk; verbose=false, engine=:ka, accelerate=true),
                        "newton (backtracking)" => wk->Dolo.newton(dm, wk; verbose=false, engine=:ka, globalization=:backtracking),
                        "accelerated (gmres)" => wk->Dolo.time_iteration(dm, wk; verbose=false, engine=:ka, accelerate=true, linsolve=:gmres),
                        "newton (gmres)" => wk->Dolo.newton(dm, wk; verbose=false, engine=:ka, globalization=:backtracking, linsolve=:gmres),
                    ]
                    @testset "$name" for (name, solve) in variants
                        sol = solve(Dolo.newton_workspace(dm; interp_mode=interp, dest=CuArray, linsolve=:gmres))
                        @test Dolo.converged(sol)
                        @test sol.dr.values.data isa CuArray
                        x = on_grid(Dolo.adapt(Array, sol.dr))
                        @test maximum(maximum(abs.(a-b)) for (a,b) in zip(x, x_ref)) < 1e-6
                    end

                end

            end

        end

        @testset "integration_B" begin

            # 8 controls, 9 exogenous nodes, 3 endogenous states
            model = include("$(pkgdir(Dolo))/examples/ymodels/integration_B.jl")
            dm = Dolo.discretize(model, Dict(:endo => [5, 5, 5]))
            cpu = Dolo.newton_workspace(dm; interp_mode=:cubic)
            cpu.x0.data .= [x .* (1 .+ 0.01*sin.((1:8) .* n)) for (n,x) in enumerate(cpu.x0.data)]
            gpu = Dolo.adapt(CuArray, cpu)
            @test gpu.L.D isa Dolo.BlockVector
            engine = Dolo.select_engine(:ka, gpu.x0)
            for wk in (cpu, gpu)
                Dolo.fit!(wk.φ, wk.x0)
            end
            @test maxdiff(Dolo.F!(cpu.r0, dm, cpu.x0, cpu.φ).data, Dolo.F!(gpu.r0, dm, gpu.x0, gpu.φ, engine).data) < 1e-12
            @test maxdiff(Dolo.dF_1!(cpu.J, dm, cpu.x0, cpu.φ).data, Dolo.dF_1!(gpu.J, dm, gpu.x0, gpu.φ, engine).data) < 1e-10
            Dolo.dF_2!(cpu.L, dm, cpu.x0, cpu.φ)
            Dolo.dF_2!(gpu.L, dm, gpu.x0, gpu.φ, engine)
            @test maxdiff(cpu.L, gpu.L) < 1e-10
            @test maxdiff(Dolo.mul!(cpu.r0, cpu.L, cpu.x0).data, Dolo.mul!(gpu.r0, gpu.L, gpu.x0, engine).data) < 1e-10
            @test maxdiff(Dolo.ldiv!(cpu.J, cpu.L), Dolo.ldiv!(gpu.J, gpu.L, engine)) < 1e-10

            # a few time iteration steps, then Newton: same solution as on the CPU
            function solve(dest, engine)
                wk = Dolo.newton_workspace(dm; interp_mode=:cubic, dest=dest, linsolve=:gmres)
                Dolo.time_iteration(dm, wk; T=25, verbose=false, engine=engine)
                Dolo.newton(dm, wk; verbose=false, engine=engine, linsolve=:gmres)
            end
            ref = solve(Array, :threads)
            sol = solve(CuArray, :ka)
            @test Dolo.converged(sol)
            @test maxdiff(sol.dr.values.data, ref.dr.values.data) < 1e-8

        end

    end

end
