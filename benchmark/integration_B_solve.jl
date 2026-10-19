# Complete solution of the integration_B model: a few time iteration steps from the initial
# guess, then Newton, on CPU threads and on the GPU. Prints the timings of each part.
#
#     julia -t auto --project=benchmark benchmark/integration_B_solve.jl [--sizes=20,40] [--precisions=Float64,Float32]
#                                                    [--engines=threads,gpu] [--linsolve=gmres|neumann]
# On Intel GPUs:
#     julia -t auto --project=benchmark/oneapi benchmark/integration_B_solve.jl --backend=oneapi
#
# Tolerances depend on the precision (`tol_ε = 1e-8` cannot be reached in Float32).

using Dolo
using Printf

include("gpu_backend.jl")

option(name, default) = let a = findfirst(startswith("--$name="), ARGS)
    a === nothing ? default : split(split(ARGS[a], "=")[2], ",")
end

const SIZES = parse.(Int, option("sizes", ["20"]))
const PRECISIONS = [Dict("Float64" => Float64, "Float32" => Float32)[p] for p in option("precisions", ["Float64", "Float32"])]
const ENGINES = filter(e -> e == "threads" || GPU.functional, option("engines", ["threads", "gpu"]))
const LINSOLVE = Symbol(only(option("linsolve", ["gmres"])))
const WARMUP = 5     # time iteration steps before Newton

model64 = include(joinpath(pkgdir(Dolo), "examples", "ymodels", "integration_B.jl"))

function solve(dm, ename, T; timer=false)
    engine, dest = ename == "gpu" ? (:ka, GPU.AT) : (:threads, Array)
    tol_ε = T == Float64 ? 1e-8 : 1e-4
    linsolve_tol = T == Float64 ? 1e-10 : 1e-5
    wk = Dolo.newton_workspace(dm; dest=dest, interp_mode=:cubic, linsolve=LINSOLVE)
    t_ti = @elapsed Dolo.time_iteration(dm, wk; T=WARMUP, verbose=false, engine=engine, tol_ε=tol_ε)
    t_nw = @elapsed (sol = Dolo.newton(dm, wk; verbose=false, engine=engine, tol_ε=tol_ε,
        linsolve_tol=linsolve_tol, linsolve=LINSOLVE, timer=timer))
    return sol, t_ti, t_nw
end

println("Julia ", VERSION, ", ", Threads.nthreads(), " threads, ", Sys.cpu_info()[1].model)
GPU.functional && println("GPU ($(GPU.name)): ", GPU.description)
println("Linear solver: ", LINSOLVE)

for T in PRECISIONS
    local model = T == Float64 ? model64 : Dolo.convert_precision(T, model64)
    engines = filter(e -> e == "threads" || T == Float32 || GPU.float64, ENGINES)
    # compilation, on a small grid
    dms = Dolo.discretize(model, Dict(:endo => [4, 4, 4]))
    for e in engines
        Base.invokelatest(solve, dms, e, T)
    end
    for n in SIZES, e in engines
        dm = Dolo.discretize(model, Dict(:endo => [n, n, n]))
        try
            sol, t_ti, t_nw = Base.invokelatest(solve, dm, e, T; timer=true)
            @printf("\n%s, %d³ × 9 = %d points, %s: %d time iteration steps %.2f s, Newton %.2f s (%d iterations, ε=%.1e, converged=%s)\n",
                T, n, length(dm.grid), e, WARMUP, t_ti, t_nw, sol.iterations, sol.ε, Dolo.converged(sol))
            show(IOContext(stdout, :compact => true), sol.timer; allocations=false, sortby=:firstexec)
            println()
        catch err
            @printf("\n%s, %d³, %s: failed: %s\n", T, n, e, sprint(showerror, err)[1:min(end, 300)])
        end
        GC.gc(true); GPU.functional && GPU.reclaim()
    end
end
