# Kernels of the integration_B model (8 controls, 5 states, 9 exogenous nodes) on CPU and GPU:
# residuals, their derivatives, and application of the linear operator L.
#
#     julia -t auto --project=benchmark benchmark/integration_B.jl [--sizes=10,20,30] [--precisions=Float64,Float32]
#                                                                   [--engines=serial,threads,gpu]

using Dolo
using CUDA
using Printf

option(name, default) = let a = findfirst(startswith("--$name="), ARGS)
    a === nothing ? default : split(split(ARGS[a], "=")[2], ",")
end

const SIZES = parse.(Int, option("sizes", ["10", "20", "30", "40"]))
const PRECISIONS = [Dict("Float64" => Float64, "Float32" => Float32)[p] for p in option("precisions", ["Float64", "Float32"])]
const SERIAL_MAX = 20    # the serial engine is only timed up to this size

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

model64 = include(joinpath(pkgdir(Dolo), "examples", "ymodels", "integration_B.jl"))

sync(dest) = dest === CuArray && CUDA.synchronize()

"Minimum time of `f()` over repeated runs (about one second in total, at least 3 runs)."
function timeit(f, dest)
    f(); sync(dest)
    ts = Float64[]
    while length(ts) < 3 || (sum(ts) < 1.0 && length(ts) < 20)
        push!(ts, @elapsed (f(); sync(dest)))
    end
    return minimum(ts)
end

"Workspace on the device of `dest`, with φ fitted, J and L computed (from a slightly perturbed guess)."
function workspace(dm, engine, dest)
    wk = Dolo.newton_workspace(dm; interp_mode=:cubic)
    T = eltype(eltype(wk.x0))
    wk.x0.data .= [x .* (1 .+ T(0.01)*sin.(T.(1:length(x)) .* n)) for (n, x) in enumerate(wk.x0.data)]
    wk = Dolo.adapt(dest, wk)
    Dolo.fit!(wk.φ, wk.x0)
    Dolo.dF_1!(wk.J, dm, wk.x0, wk.φ, engine)
    Dolo.dF_2!(wk.L, dm, wk.x0, wk.φ, engine)
    return wk
end

maxdiff(a, b) = maximum(maximum(abs.(u - v)) for (u, v) in zip(a, b))
host(x) = Array(x)

"Largest difference between the results of each kernel on the GPU and on the CPU (threads)."
function check_gpu(dm)
    wc = workspace(dm, Dolo.Threaded(), Array)
    wg = workspace(dm, Dolo.KA(CUDABackend()), CuArray)
    d = Dict{String,Float64}()
    d["F!"] = maxdiff(host(Dolo.F!(wg.r0, dm, wg.x0, wg.φ, Dolo.KA(CUDABackend())).data), Dolo.F!(wc.r0, dm, wc.x0, wc.φ, Dolo.Threaded()).data)
    d["dF_1!"] = maxdiff(host(wg.J.data), wc.J.data)
    d["dF_2!"] = maximum(maxdiff(map(b -> b.F_x, u), map(b -> b.F_x, v)) for (u, v) in zip(host(wg.L.D), wc.L.D))
    d["L.v"] = maxdiff(host(Dolo.mul!(wg.r0, wg.L, wg.x0, Dolo.KA(CUDABackend())).data), Dolo.mul!(wc.r0, wc.L, wc.x0, Dolo.Threaded()).data)
    return d
end

println("Julia ", VERSION, ", ", Threads.nthreads(), " threads, ", Sys.cpu_info()[1].model)
CUDA.functional() && println("GPU: ", CUDA.name(CUDA.device()), ", ", round(CUDA.totalmem(CUDA.device())/2^30, digits=1), " GiB")
println("Engines: ", join(first.(ENGINES), ", "))

rows = []

for T in PRECISIONS
    local model = T == Float64 ? model64 : Dolo.convert_precision(T, model64)
    for n in SIZES
        dm = Dolo.discretize(model, Dict(:endo => [n, n, n]))
        N = length(dm.grid)
        wk = Dolo.newton_workspace(Dolo.discretize(model, Dict(:endo => [2, 2, 2])))
        bytes_L = sizeof(eltype(wk.L.D))*N
        @printf("\n%s, %d³ endogenous points × 9 = %d points (L: %.2f GiB)\n", T, n, N, bytes_L/2^30)
        if CUDA.functional() && n == first(SIZES)
            for (k, v) in sort(collect(check_gpu(dm)))
                @printf("  gpu vs cpu: %-6s max diff %.2e\n", k, v)
            end
        end
        for (ename, (engine, dest)) in ENGINES
            ename == "serial" && n > SERIAL_MAX && continue
            local wk
            try
                wk = workspace(dm, engine, dest)
            catch e
                println("  ", ename, ": workspace failed: ", sprint(showerror, e)[1:min(end, 200)])
                continue
            end
            for (kname, f) in KERNELS
                t = try
                    if kname == "J\\L"
                        # J\L modifies L: start from the same L each time (the copy is not timed)
                        L0 = deepcopy(wk.L)
                        minimum(begin
                            copyto!(wk.L.D, L0.D); sync(dest)
                            @elapsed (f(dm, wk, engine); sync(dest))
                        end for _ in 1:5)
                    else
                        timeit(() -> f(dm, wk, engine), dest)
                    end
                catch e
                    println("  ", ename, " ", kname, " failed: ", sprint(showerror, e)[1:min(end, 200)])
                    NaN
                end
                push!(rows, (; T, n, N, engine=ename, kernel=kname, t))
                @printf("  %-8s %-6s %10.3f ms  %8.1f ns/point\n", ename, kname, t*1e3, t/N*1e9)
            end
            wk = nothing; GC.gc(true); CUDA.functional() && CUDA.reclaim()
        end
    end
end

# summary: speed-ups of the GPU relative to threads
println("\n| precision | points | operation | ", join(first.(ENGINES), " | "), " | gpu/threads |")
println("|---|---|---|", repeat("---|", length(ENGINES) + 1))
for key in unique((r.T, r.N, r.kernel) for r in rows)
    rs = Dict(r.engine => r.t for r in rows if (r.T, r.N, r.kernel) == key)
    ts = [haskey(rs, e) ? @sprintf("%.3f ms", rs[e]*1e3) : "" for e in first.(ENGINES)]
    su = haskey(rs, "gpu") && haskey(rs, "threads") ? @sprintf("%.1f×", rs["threads"]/rs["gpu"]) : ""
    println("| ", key[1], " | ", key[2], " | ", key[3], " | ", join(ts, " | "), " | ", su, " |")
end
