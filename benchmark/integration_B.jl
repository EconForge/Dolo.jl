# Kernels of the integration_B model (8 controls, 5 states, 9 exogenous nodes) on CPU and GPU:
# residuals, their derivatives, and application of the linear operator L.
#
#     julia -t auto --project=benchmark benchmark/integration_B.jl [--sizes=10,20,30] [--precisions=Float64,Float32]
#                                                                   [--engines=serial,threads,gpu]
# On Intel GPUs:
#     julia -t auto --project=benchmark/oneapi benchmark/integration_B.jl --backend=oneapi
#
# Results are saved in benchmark/history/ with the platform, the software and the commit (see record.jl).

using Dolo
using Printf

include("gpu_backend.jl")
include("record.jl")

option(name, default) = let a = findfirst(startswith("--$name="), ARGS)
    a === nothing ? default : split(split(ARGS[a], "=")[2], ",")
end

const SIZES = parse.(Int, option("sizes", ["10", "20", "30", "40"]))
const PRECISIONS = [Dict("Float64" => Float64, "Float32" => Float32)[p] for p in option("precisions", ["Float64", "Float32"])]
gpu_ok(T) = GPU.functional && (T == Float32 || GPU.float64)   # can the GPU run in precision T?
const SERIAL_MAX = 20    # the serial engine is only timed up to this size

# engine name (as recorded) => (engine, array type where the data lives)
const ENGINES = let
    all = ["serial" => (Dolo.Serial(), Array), "threads" => (Dolo.Threaded(), Array)]
    GPU.functional && push!(all, GPU.name => (Dolo.KA(GPU.backend), GPU.AT))
    selected = replace(option("engines", ["serial", "threads", "gpu"]), "gpu" => GPU.name)
    filter(e -> e[1] in selected, all)
end

const KERNELS = [
    "fit!"     => (dm, wk, e) -> Dolo.fit!(wk.φ, wk.x0),
    "F!"       => (dm, wk, e) -> Dolo.F!(wk.r0, dm, wk.x0, wk.φ, e),
    "dF_1!"    => (dm, wk, e) -> Dolo.dF_1!(wk.J, dm, wk.x0, wk.φ, e),
    "F+dF_1!"  => (dm, wk, e) -> Dolo.F_and_dF_1!(wk.r0, wk.J, dm, wk.x0, wk.φ, e),
    "dF_2!"    => (dm, wk, e) -> Dolo.dF_2!(wk.L, dm, wk.x0, wk.φ, e),
    "L.v"      => (dm, wk, e) -> Dolo.mul!(wk.r0, wk.L, wk.x0, e),
    "J\\L"     => (dm, wk, e) -> Dolo.ldiv!(wk.J, wk.L, e),
]

model64 = include(joinpath(pkgdir(Dolo), "examples", "ymodels", "integration_B.jl"))

sync(dest) = dest === GPU.AT && GPU.sync()

"Minimum time of `f()` over repeated runs (about one second in total, at least 3 runs)."
function timeit(f, dest)
    f(); sync(dest)
    ts = Float64[]
    while length(ts) < 3 || (sum(ts) < 1.0 && length(ts) < 20)
        push!(ts, @elapsed (f(); sync(dest)))
    end
    return minimum(ts)
end

"""
Workspace with φ fitted, J and L computed on the CPU (from a slightly perturbed guess), then moved
to the device of `dest`: a kernel that fails on the GPU does not prevent the others from being timed.
"""
function workspace(dm, dest)
    wk = Dolo.newton_workspace(dm; interp_mode=:cubic)
    T = eltype(eltype(wk.x0))
    wk.x0.data .= [x .* (1 .+ T(0.01)*sin.(T.(1:length(x)) .* n)) for (n, x) in enumerate(wk.x0.data)]
    Dolo.fit!(wk.φ, wk.x0)
    Dolo.dF_1!(wk.J, dm, wk.x0, wk.φ, Dolo.Threaded())
    Dolo.dF_2!(wk.L, dm, wk.x0, wk.φ, Dolo.Threaded())
    return Dolo.adapt(dest, wk)
end

maxdiff(a, b) = maximum(maximum(abs.(u - v)) for (u, v) in zip(a, b))
host(x) = Array(x)

"Largest difference between the results of each kernel on the GPU and on the CPU (threads); NaN if it fails."
function check_gpu(dm)
    e = Dolo.KA(GPU.backend)
    wc = workspace(dm, Array)
    wg = workspace(dm, GPU.AT)
    checks = [
        "F!" => () -> maxdiff(host(Dolo.F!(wg.r0, dm, wg.x0, wg.φ, e).data), Dolo.F!(wc.r0, dm, wc.x0, wc.φ, Dolo.Threaded()).data),
        "dF_1!" => () -> (Dolo.dF_1!(wg.J, dm, wg.x0, wg.φ, e); maxdiff(host(wg.J.data), wc.J.data)),
        "F+dF_1!" => () -> (Dolo.F_and_dF_1!(wg.r0, wg.J, dm, wg.x0, wg.φ, e); max(maxdiff(host(wg.r0.data), wc.r0.data), maxdiff(host(wg.J.data), wc.J.data))),
        "dF_2!" => () -> (Dolo.dF_2!(wg.L, dm, wg.x0, wg.φ, e);
            maximum(maxdiff(map(b -> b.F_x, u), map(b -> b.F_x, v)) for (u, v) in zip(host(wg.L.D), wc.L.D))),
        "L.v" => () -> maxdiff(host(Dolo.mul!(wg.r0, wg.L, wg.x0, e).data), Dolo.mul!(wc.r0, wc.L, wc.x0, Dolo.Threaded()).data),
    ]
    return Dict(k => (try f() catch; NaN end) for (k, f) in checks)
end

this_run = run_metadata("integration_B"; gpu=GPU, options=Dict{String,Any}("args" => join(ARGS, " ")))

println("Julia ", VERSION, ", ", Threads.nthreads(), " threads, ", Sys.cpu_info()[1].model)
GPU.functional && println("GPU ($(GPU.name)): ", GPU.description)
println("Engines: ", join(first.(ENGINES), ", "))
println("Results: ", relpath(this_run.path))

rows = NamedTuple[]

for T in PRECISIONS
    local model = T == Float64 ? model64 : Dolo.convert_precision(T, model64)
    for n in SIZES
        dm = Dolo.discretize(model, Dict(:endo => [n, n, n]))
        N = length(dm.grid)
        wk = Dolo.newton_workspace(Dolo.discretize(model, Dict(:endo => [2, 2, 2])))
        bytes_L = sizeof(eltype(wk.L.D))*N
        @printf("\n%s, %d³ endogenous points × 9 = %d points (L: %.2f GiB)\n", T, n, N, bytes_L/2^30)
        diffs = Dict{String,Float64}()
        if gpu_ok(T) && n == first(SIZES) && any(e -> e[1] == GPU.name, ENGINES)
            diffs = check_gpu(dm)
            for (k, v) in sort(collect(diffs))
                @printf("  %s vs cpu: %-6s max diff %.2e\n", GPU.name, k, v)
            end
        end
        for (ename, (engine, dest)) in ENGINES
            ename == "serial" && n > SERIAL_MAX && continue
            ename == GPU.name && !gpu_ok(T) && continue
            local wk
            try
                wk = workspace(dm, dest)
            catch e
                println("  ", ename, ": workspace failed: ", short_error(e))
                for (kname, _) in KERNELS
                    push!(rows, (; benchmark="kernels", model="integration_B", precision=string(T), N, grid="$(n)³×9",
                        engine=ename, operation=kname, time=NaN, status="error", error=short_error(e)))
                end
                continue
            end
            for (kname, f) in KERNELS
                t, status, err = try
                    t = if kname == "J\\L"
                        # J\L modifies L: start from the same L each time (the copy is not timed)
                        L0 = deepcopy(wk.L)
                        minimum(begin
                            copyto!(wk.L.D, L0.D); sync(dest)
                            @elapsed (f(dm, wk, engine); sync(dest))
                        end for _ in 1:5)
                    else
                        timeit(() -> f(dm, wk, engine), dest)
                    end
                    t, "ok", nothing
                catch e
                    NaN, "error", short_error(e)
                end
                maxdiff_ = ename == GPU.name ? get(diffs, kname, nothing) : nothing
                push!(rows, (; benchmark="kernels", model="integration_B", precision=string(T), N, grid="$(n)³×9",
                    engine=ename, operation=kname, time=t, status, error=err, maxdiff=maxdiff_))
                status == "ok" ?
                    @printf("  %-8s %-6s %10.3f ms  %8.1f ns/point\n", ename, kname, t*1e3, t/N*1e9) :
                    @printf("  %-8s %-6s failed: %s\n", ename, kname, err)
            end
            wk = nothing; GC.gc(true); GPU.functional && GPU.reclaim()
        end
        save(this_run, rows)
    end
end

println()
print_table(as_history(this_run, rows))
println("\nSaved in ", relpath(this_run.path))
