# Performance benchmarks of Dolo.jl: models × grid sizes × engines.
#
#     julia -t auto --project=benchmark -e 'using Pkg; Pkg.instantiate()'
#     julia -t auto --project=benchmark benchmark/run.jl [options]
#
# Options:
#     --quick                      smaller grids, for a fast check
#     --models=rbc_mc,rbc_ar1      models (from examples/ymodels)
#     --engines=serial,threads,gpu engines (gpu only if CUDA is functional)
#     --only=kernels|solvers       one category of benchmarks
#
# Two categories of benchmarks:
# - kernels: one call of each vectorized operation on the whole grid
# - solvers: complete solution of the model, from the initial guess
#
# Results are printed as tables, and saved (markdown and CSV) in benchmark/results/.

using Dolo
using BenchmarkTools
using CUDA
using Printf, Dates

## Configuration

option(name, default) = let a = findfirst(startswith("--$name="), ARGS)
    a === nothing ? default : split(split(ARGS[a], "=")[2], ",")
end

const QUICK = "--quick" in ARGS
const MODELS = option("models", ["rbc_mc", "rbc_ar1", "consumption_savings_iid", "rbc_iid"])
const ONLY = option("only", ["kernels", "solvers"])
const KERNEL_SIZES = QUICK ? [10^3, 10^4] : [10^3, 10^4, 10^5, 10^6]
const SOLVER_SIZES = QUICK ? [10^3] : [10^3, 10^4]

# engine name => (engine, array type where the data lives)
const ENGINES = let
    all = ["serial" => (Dolo.Serial(), Array), "threads" => (Dolo.Threaded(), Array)]
    CUDA.functional() && push!(all, "gpu" => (Dolo.KA(CUDABackend()), CuArray))
    filter(e -> e[1] in option("engines", first.(all)), all)
end

const KERNELS = [
    "fit!"  => (dm, wk, e) -> Dolo.fit!(wk.φ, wk.x0),
    "F!"    => (dm, wk, e) -> Dolo.F!(wk.r0, dm, wk.x0, wk.φ, e),
    "dF_1!" => (dm, wk, e) -> Dolo.dF_1!(wk.J, dm, wk.x0, wk.φ, e),
    "dF_2!" => (dm, wk, e) -> Dolo.dF_2!(wk.L, dm, wk.x0, wk.φ, e),
    "L.v"   => (dm, wk, e) -> Dolo.mul!(wk.r0, wk.L, wk.x0, e),
    "J\\L"  => (dm, wk, e) -> Dolo.ldiv!(wk.J, wk.L, e),
]

const SOLVERS = [
    "time iteration" => (dm, wk, e) -> Dolo.time_iteration(dm, wk; verbose=false, engine=e),
    "accelerated TI" => (dm, wk, e) -> Dolo.time_iteration(dm, wk; verbose=false, engine=e, accelerate=true),
    "newton"         => (dm, wk, e) -> Dolo.newton(dm, wk; verbose=false, engine=e, globalization=:backtracking),
]


## Setup

load_model(name) = include(joinpath(pkgdir(Dolo), "examples", "ymodels", "$name.jl"))

"Discretize `model` with about `N` grid points (exogenous nodes × endogenous points)."
function discretize(model, N)
    g = Dolo.discretize(model).grid
    n_exo, cgrid = g isa Dolo.CGrid ? (1, g) : (length(g.grids[1]), g.grids[2])
    d = length(cgrid.ranges)
    n = round(Int, (N / n_exo)^(1/d))
    return Dolo.discretize(model, Dict(:endo => fill(n, d)))
end

"Workspace on the device of `dest`, ready for the kernels (φ fitted, J and L computed)."
function workspace(dm, engine, dest)
    wk = Dolo.newton_workspace(dm; dest=dest, interp_mode=:cubic)
    Dolo.fit!(wk.φ, wk.x0)
    Dolo.dF_1!(wk.J, dm, wk.x0, wk.φ, engine)
    Dolo.dF_2!(wk.L, dm, wk.x0, wk.φ, engine)
    return wk
end

function free_memory()
    GC.gc(true)
    CUDA.functional() && CUDA.reclaim()
end


## Measurements

const DIR = mkpath(joinpath(@__DIR__, "results"))
const STAMP = Dates.format(now(), "yyyy-mm-dd_HHMM")
const CSV = joinpath(DIR, "$STAMP.csv")

results = NamedTuple[]

# `status` is "ok", "no convergence" or "error".
# Results are saved as soon as they are measured, so that nothing is lost if a run is interrupted.
function record!(category, model, N, engine, name; time=NaN, allocs=0, memory=0, status="ok")
    r = (; category, model, N, engine, name, time, allocs, memory, status)
    isempty(results) && open(io -> println(io, join(keys(r), ",")), CSV, "w")
    push!(results, r)
    open(io -> println(io, join(values(r), ",")), CSV, "a")
    @printf("  %-8s %-24s %-8s %-16s %12s\n", category, model, engine, name,
        status == "ok" ? BenchmarkTools.prettytime(time*1e9) : status)
end

function record!(category, model, N, engine, name, trial::BenchmarkTools.Trial)
    t = minimum(trial)
    record!(category, model, N, engine, name; time=t.time/1e9, allocs=t.allocs, memory=t.memory)
end

# runs `f()`; if it fails, records the error and continues
function attempt(f, category, model, N, engine, name)
    try
        f()
    catch e
        @error "$category / $model / $N / $engine / $name failed" exception=(e, catch_backtrace())
        record!(category, model, N, engine, name; status="error")
    end
end

function bench_kernels(model_name, model, N)
    dm = discretize(model, N)
    for (ename, (engine, dest)) in ENGINES
        wk = workspace(dm, engine, dest)
        for (kname, f) in KERNELS
            attempt("kernels", model_name, length(dm.grid), ename, kname) do
                # `J\L` modifies L, so it starts from a copy at each evaluation
                trial = kname == "J\\L" ?
                    (@benchmark $f($dm, w, $engine) setup=(w = (; J=$wk.J, L=deepcopy($wk.L))) evals=1 seconds=2) :
                    (@benchmark $f($dm, $wk, $engine) seconds=2)
                record!("kernels", model_name, length(dm.grid), ename, kname, trial)
            end
        end
        wk = nothing; free_memory()
    end
end

function bench_solvers(model_name, model, N)
    dm = discretize(model, N)
    for (ename, (engine, dest)) in ENGINES
        wk0 = Dolo.newton_workspace(dm; dest=dest, interp_mode=:cubic)
        for (sname, f) in SOLVERS
            attempt("solvers", model_name, length(dm.grid), ename, sname) do
                sol = f(dm, deepcopy(wk0), engine)      # compilation, and check
                if Dolo.converged(sol)
                    # each run starts from a fresh copy of the initial workspace
                    trial = @benchmark $f($dm, wk, $engine) setup=(wk = deepcopy($wk0)) evals=1 samples=3 seconds=1
                    record!("solvers", model_name, length(dm.grid), ename, sname, trial)
                else
                    record!("solvers", model_name, length(dm.grid), ename, sname; status="no convergence")
                end
            end
        end
        wk0 = nothing; free_memory()
    end
end


## Report

speedup(r) = r >= 1 ? @sprintf("%.1f×", r) : @sprintf("%.2f×", r)

function table(io, category)
    rs = filter(r -> r.category == category, results)
    isempty(rs) && return
    engines = unique(r.engine for r in rs)
    println(io, "\n### ", category, "\n")
    header = ["model", "points", "operation", ("$e" for e in engines)..., ("$e/serial" for e in engines if e != "serial")...]
    println(io, "| ", join(header, " | "), " |")
    println(io, "|", repeat("---|", length(header)))
    for key in unique((r.model, r.N, r.name) for r in rs)
        row = Dict(r.engine => r for r in rs if (r.model, r.N, r.name) == key)
        ok(e) = haskey(row, e) && row[e].status == "ok"
        times = [ok(e) ? BenchmarkTools.prettytime(row[e].time*1e9) : haskey(row, e) ? row[e].status : ""
                 for e in engines]
        speedups = [ok(e) && ok("serial") ? speedup(row["serial"].time/row[e].time) : ""
                    for e in engines if e != "serial"]
        println(io, "| ", join([key[1], key[2], key[3], times..., speedups...], " | "), " |")
    end
end

function report(io)
    println(io, "## Dolo.jl benchmarks, ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
    println(io)
    println(io, "- Julia ", VERSION, ", ", Threads.nthreads(), " threads, ", Sys.cpu_info()[1].model)
    CUDA.functional() && println(io, "- GPU: ", CUDA.name(CUDA.device()), " (Float64)")
    println(io, "- times are minima over samples; speed-ups are relative to the serial engine")
    table(io, "kernels")
    table(io, "solvers")
end


## Run

Threads.nthreads() == 1 && @warn "Julia runs with one thread: start it with `julia -t auto` to benchmark the threaded engine."
println("Engines: ", join(first.(ENGINES), ", "), (CUDA.functional() ? "" : "  (CUDA not functional: no gpu)"))

for name in MODELS
    model = load_model(name)
    # the model's methods were just defined: `invokelatest` makes them visible
    "kernels" in ONLY && for N in KERNEL_SIZES
        Base.invokelatest(bench_kernels, name, model, N)
    end
    "solvers" in ONLY && for N in SOLVER_SIZES
        Base.invokelatest(bench_solvers, name, model, N)
    end
end

report(stdout)

open(io -> report(io), joinpath(DIR, "$STAMP.md"), "w")
println("\nSaved in ", joinpath(DIR, "$STAMP.{md,csv}"))
