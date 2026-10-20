# Complete solution of the integration_B model from the initial guess, on CPU threads and on the GPU,
# with the three algorithms of run.jl: time iteration, accelerated time iteration and Newton
# (with backtracking).
#
#     julia -t auto --project=benchmark benchmark/integration_B_solve.jl [--sizes=20,40] [--precisions=Float64,Float32]
#                    [--engines=threads,gpu] [--algos=time_iteration,accelerated,newton] [--linsolve=gmres|neumann]
# On Intel GPUs:
#     julia -t auto --project=benchmark/oneapi benchmark/integration_B_solve.jl --backend=oneapi
#
# Tolerances depend on the precision (`tol_ε = 1e-8` cannot be reached in Float32).
# Results are saved in benchmark/history/ with the platform, the software and the commit (see record.jl).

using Dolo
using Printf

include("gpu_backend.jl")
include("record.jl")

option(name, default) = let a = findfirst(startswith("--$name="), ARGS)
    a === nothing ? default : split(split(ARGS[a], "=")[2], ",")
end

const SIZES = parse.(Int, option("sizes", ["20"]))
const PRECISIONS = [Dict("Float64" => Float64, "Float32" => Float32)[p] for p in option("precisions", ["Float64", "Float32"])]
const ENGINES = filter(e -> e == "threads" || GPU.functional,
    replace(option("engines", ["threads", "gpu"]), "gpu" => GPU.name))
const LINSOLVE = Symbol(only(option("linsolve", ["gmres"])))
const ALGOS = let names = Dict("time_iteration" => "time iteration", "accelerated" => "accelerated TI", "newton" => "newton")
    [names[a] for a in option("algos", ["time_iteration", "accelerated", "newton"])]
end

model64 = include(joinpath(pkgdir(Dolo), "examples", "ymodels", "integration_B.jl"))

"""
Solve `dm` with `algo` from the initial guess; returns the solution and the time of the solver
(the allocation of the workspace is not timed). `short=true` stops after a few iterations: used to
compile the methods on a small grid.
"""
function solve(dm, algo, ename, T; short=false)
    engine, dest = ename == "threads" ? (:threads, Array) : (:ka, GPU.AT)
    tol_ε = T == Float64 ? 1e-8 : 1e-4
    linsolve_tol = T == Float64 ? 1e-10 : 1e-5
    if algo == "time iteration"
        wk = Dolo.time_iteration_workspace(dm; dest=dest, interp_mode=:cubic)
        t = @elapsed sol = Dolo.time_iteration(dm, wk; verbose=false, engine=engine, tol_ε=tol_ε,
            T=short ? 3 : 500, timer=true)
    elseif algo == "accelerated TI"
        wk = Dolo.time_iteration_workspace(dm; dest=dest, interp_mode=:cubic, accelerate=true, linsolve=LINSOLVE)
        # acceleration starts after `accelerate_wait=10` steps
        t = @elapsed sol = Dolo.time_iteration(dm, wk; verbose=false, engine=engine, tol_ε=tol_ε,
            accelerate=true, linsolve=LINSOLVE, linsolve_tol=linsolve_tol, T=short ? 12 : 500, timer=true)
    elseif algo == "newton"
        wk = Dolo.newton_workspace(dm; dest=dest, interp_mode=:cubic, linsolve=LINSOLVE)
        t = @elapsed sol = Dolo.newton(dm, wk; verbose=false, engine=engine, tol_ε=tol_ε,
            linsolve=LINSOLVE, linsolve_tol=linsolve_tol, globalization=:backtracking, T=short ? 2 : 20, timer=true)
    end
    return sol, t
end

this_run = run_metadata("integration_B_solve"; gpu=GPU, options=Dict{String,Any}("args" => join(ARGS, " "),
    "linsolve" => string(LINSOLVE)))

println("Julia ", VERSION, ", ", Threads.nthreads(), " threads, ", Sys.cpu_info()[1].model)
GPU.functional && println("GPU ($(GPU.name)): ", GPU.description)
println("Linear solver: ", LINSOLVE, "; algorithms: ", join(ALGOS, ", "), "; engines: ", join(ENGINES, ", "))
println("Results: ", relpath(this_run.path))

rows = NamedTuple[]

for T in PRECISIONS
    local model = T == Float64 ? model64 : Dolo.convert_precision(T, model64)
    engines = filter(e -> e == "threads" || T == Float32 || GPU.float64, ENGINES)
    dms = Dolo.discretize(model, Dict(:endo => [6, 6, 6]))
    for algo in ALGOS, e in engines
        # compilation, on a small grid (failures are reported on the larger grids)
        try Base.invokelatest(solve, dms, algo, e, T; short=true) catch end
        for n in SIZES
            dm = Dolo.discretize(model, Dict(:endo => [n, n, n]))
            N = length(dm.grid)
            row = (; benchmark="solve", model="integration_B", precision=string(T), N, grid="$(n)³×9", engine=e, operation=algo)
            try
                sol, t = Base.invokelatest(solve, dm, algo, e, T)
                ok = Dolo.converged(sol)
                push!(rows, merge(row, (; time=t, status=ok ? "ok" : "no convergence", iterations=sol.iterations,
                    residual=sol.ε, sections=timer_sections(sol.timer))))
                @printf("\n%s, %d³ × 9 = %d points, %s, %s: %.2f s (%d iterations, ε=%.1e, converged=%s)\n",
                    T, n, N, e, algo, t, sol.iterations, sol.ε, ok)
                show(IOContext(stdout, :compact => true), sol.timer; allocations=false, sortby=:firstexec)
                println()
            catch err
                push!(rows, merge(row, (; time=NaN, status="error", error=short_error(err))))
                @printf("\n%s, %d³, %s, %s: failed: %s\n", T, n, e, algo, short_error(err))
            end
            save(this_run, rows)
            GC.gc(true); GPU.functional && GPU.reclaim()
        end
    end
end

println()
print_table(as_history(this_run, rows))
println("\nSaved in ", relpath(this_run.path))
