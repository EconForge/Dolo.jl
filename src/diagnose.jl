# Diagnostic of the computing resources (CPU, memory, GPUs), and of what a model needs from them:
# size of the grid, memory used by the solvers in each precision, size of the kernel arguments,
# compatibility with GPUs and with Float32.
#
#     Dolo.diagnose()                                                  # this machine
#     Dolo.diagnose(model; discretization=Dict(:endo => [20, 20, 20])) # and a model
#     Dolo.diagnose(dmodel)                                            # and a discretized model

import Adapt
import Printf: @sprintf


## GPUs, as described by the extension of each GPU package

# package name => function returning a description of its current device:
# (; package, functional, reason, device, memory_total, memory_free, integrated, float64,
#    float64_ratio, details, param_limit)
const GPU_PACKAGES = Dict{String,Function}()

"Called by the extension of a GPU package (e.g. DoloCUDAExt) to describe its device."
register_gpu!(name::String, info::Function) = (GPU_PACKAGES[name] = info; nothing)

const KNOWN_GPU_PACKAGES = ("CUDA", "oneAPI")

# hardware present, as far as can be told without the GPU package (Linux only)
function gpu_hardware(package)
    Sys.islinux() || return false
    if package == "CUDA"
        return isfile("/proc/driver/nvidia/version") || isfile("/dev/nvidia0")
    elseif package == "oneAPI"
        cards = try filter(c -> occursin(r"^card\d+$", c), readdir("/sys/class/drm")) catch; String[] end
        vendor(c) = try strip(read("/sys/class/drm/$c/device/vendor", String)) catch; "" end
        return any(c -> vendor(c) == "0x8086", cards)
    end
    return false
end

function gpu_infos()
    map(KNOWN_GPU_PACKAGES) do name
        if haskey(GPU_PACKAGES, name)
            info = try
                GPU_PACKAGES[name]()
            catch e
                (; functional=false, reason="error while querying the device: " * sprint(showerror, e))
            end
            merge((; package=name, loaded=true), info)
        else
            (; package=name, loaded=false, functional=false, hardware=gpu_hardware(name))
        end
    end
end


## Sizes of the workspaces, measured on small grids

# bytes of the arrays contained in x (on the CPU)
hostbytes(x::Array) = sizeof(x)
hostbytes(x::Union{Tuple,NamedTuple}) = sum(hostbytes, values(x); init=0)
hostbytes(x) = isbits(x) ? 0 : sum((hostbytes(getfield(x, f)) for f in fieldnames(typeof(x))); init=0)

# the same discretized model, with k points in each endogenous dimension
function with_endo_points(dmodel::DYModel, k)
    endo(g::CGrid) = CGrid(map(r -> (r[1], r[2], k), g.ranges))
    grid = dmodel.grid isa CGrid ? endo(dmodel.grid) : PGrid(dmodel.grid.grids[1], endo(dmodel.grid.grids[2]))
    return DYModel(dmodel.model, grid, dmodel.dproc)
end

endo_grid(dmodel) = dmodel.grid isa CGrid ? dmodel.grid : dmodel.grid.grids[2]
n_exo_nodes(dmodel) = dmodel.grid isa CGrid ? 1 : length(dmodel.grid.grids[1])

"""
Bytes used by each part of the Newton workspace (with GMRES) of `dmodel`. All are linear in
the number of grid points, except the coefficients of the splines (`coefficients`), which
have padding in each dimension.
"""
function workspace_parts(dmodel; interp_mode=:cubic, krylov_dim=50)
    wk = newton_workspace(dmodel; interp_mode, linsolve=:gmres, krylov_dim)
    return (;
        vectors = sum(hostbytes, (wk.x0, wk.x1, wk.x2, wk.r0, wk.dx)),
        J = hostbytes(wk.J),
        values = hostbytes(wk.φ.values) + hostbytes(wk.L.φ.values),  # (φ.values: a copy on GPUs)
        L = hostbytes(wk.L.D),
        neumann = hostbytes(wk.memn),
        gmres = hostbytes(wk.gmres),
        coefficients = hostbytes(wk.φ.itp) + hostbytes(wk.L.φ.itp),   # φ and its copy in L
    ), wk
end

# memory of each solver, from the parts of the workspace (half of `values` and `coefficients` is in L)
solver_bytes(p) = [
    "time iteration" => p.vectors + p.J + (p.values + p.coefficients) ÷ 2,
    "accelerated TI / Newton (Neumann)" => p.vectors + p.J + p.values + p.coefficients + p.L + p.neumann,
    "accelerated TI / Newton (GMRES)" => p.vectors + p.J + p.values + p.coefficients + p.L + p.neumann + p.gmres,
]

# number of spline coefficients per exogenous node, for the endogenous grid g
n_coefficients(g::CGrid, pad) = prod(r[3] + pad for r in g.ranges)

# Layout of a GPU array once passed to a kernel (CUDA, oneAPI): pointer, length, dimensions.
struct ArgLayout{T,N}
    ptr::UInt64
    maxsize::Int
    dims::NTuple{N,Int}
end
struct ArgSizeAdaptor end
Adapt.adapt_storage(::ArgSizeAdaptor, a::Array{T,N}) where {T,N} = ArgLayout{T,N}(0, length(a), size(a))

"Bytes of the arguments of the main kernels, once the data is on a GPU."
function kernel_arguments(dmodel, wk)
    argsize(args) = (a = Adapt.adapt(ArgSizeAdaptor(), args); isbits(a) ? sizeof(a) : missing)
    return (;
        residuals = argsize((wk.r0.data, dmodel, wk.x0, wk.φ)),   # F!, dF_1!, dF_2!
        operator = argsize((wk.r0.data, wk.L.D, wk.L.φ)),         # L.v
    )
end

"Measures for `dmodel` in its own precision: memory per solver, kernel arguments, checks."
function measure(dmodel; interp_mode=:cubic, krylov_dim=50)
    d = length(endo_grid(dmodel).ranges)
    dm1, dm2 = with_endo_points(dmodel, 4), with_endo_points(dmodel, 5)
    (p1, wk), (p2, _) = workspace_parts(dm1; interp_mode, krylov_dim), workspace_parts(dm2; interp_mode, krylov_dim)
    N1, N2, N = length(dm1.grid), length(dm2.grid), length(dmodel.grid)
    # linear parts are extrapolated, the coefficients of the splines are computed for the grid
    extrapolate(b1, b2) = round(Int, b1 + (b2 - b1)/(N2 - N1)*(N - N1))
    pad = interp_mode == :cubic ? 2 : 0
    coefficients = p1.coefficients ÷ n_coefficients(endo_grid(dm1), pad) * n_coefficients(endo_grid(dmodel), pad)
    parts = merge(map(extrapolate, p1, p2), (; coefficients))
    solvers = solver_bytes(parts)
    # residuals at the initial guess
    fit!(wk.φ, wk.x0)
    ε = try norm(F!(wk.r0, dm1, wk.x0, wk.φ)) catch e; e end
    T = eltype(eltype(wk.x0.data))
    return (; T, solvers, kernel_args=kernel_arguments(dm1, wk), ε, future_nodes=nblocks(wk.L.D), dm_small=dm1, wk)
end

# precision of the results of the model's functions at a grid point, each called with arguments
# in the precision of the model (a Float64 literal promotes the result to Float64)
function result_precisions(dm, wk)
    s = point(dm.grid, 1)
    T = eltype(s.val)
    x0 = wk.x0[1]
    x = SVector{length(x0),T}(x0)
    w, S = first(τ(dm, s, x))
    r = arbitrage(dm, s, x, S, x)
    return (;
        initial_guess = eltype(x0),
        transition = eltype(S.val),
        arbitrage = eltype(r),
        complementarities = eltype(complementarities(dm.model, s, x, SVector{length(r),T}(r))),
    )
end


## Report

struct Diagnosis
    system
    gpus
    model         # nothing, or a description of the discretized model
    precisions    # [(T, measures or error)]
    advice        # [(level, text)], level in :ok, :info, :warn, :fail
end

"""
    diagnose()
    diagnose(model; discretization=Dict(), interp_mode=:cubic, krylov_dim=50)
    diagnose(dmodel; interp_mode=:cubic, krylov_dim=50)

Report on the computing resources available to Dolo (threads, memory, GPUs) and, given a
model or a discretized model, on what solving it needs: grid size, memory of each solver
in the model's precision and in Float32, size of the GPU kernel arguments, compatibility with
GPUs and Float32, with recommendations.

GPUs are described when their package is loaded (`using CUDA`, `using oneAPI`).
"""
function diagnose()
    sys, gpus = system_info(), gpu_infos()
    return Diagnosis(sys, gpus, nothing, [], nothing_to_report(general_advice(sys, gpus)))
end

diagnose(model::YModel; discretization=Dict(), kwargs...) = diagnose(discretize(model, discretization); kwargs...)

function diagnose(dmodel::DYModel; interp_mode=:cubic, krylov_dim=50)
    sys = system_info()
    gpus = gpu_infos()
    eg = endo_grid(dmodel)
    model = (;
        name = name(dmodel),
        states = variables(dmodel.model.states),
        exogenous = variables(dmodel.model.exogenous),
        controls = variables(dmodel.model.controls),
        process = string(typeof(dmodel.model.exogenous).name.name),
        n_exo = n_exo_nodes(dmodel),
        endo = map(r -> r[3], eg.ranges),
        N = length(dmodel.grid),
        isbits = isbits(dmodel),
    )
    precisions = []
    m = try measure(dmodel; interp_mode, krylov_dim) catch e; e end
    push!(precisions, (m isa Exception ? eltype(dmodel.model) : m.T, m))
    if !(m isa Exception) && m.T == Float64
        m32 = try
            d32 = convert_precision(Float32, dmodel)
            r = measure(d32; interp_mode, krylov_dim)
            merge(r, (; results = result_precisions(r.dm_small, r.wk)))
        catch e
            e
        end
        push!(precisions, (Float32, m32))
    end
    advice = nothing_to_report(vcat(general_advice(sys, gpus), model_advice(sys, gpus, model, precisions)))
    return Diagnosis(sys, gpus, model, precisions, advice)
end

system_info() = (;
    julia = string(VERSION),
    threads = Threads.nthreads(),
    cpu_threads = Sys.CPU_THREADS,
    cpu = strip(Sys.cpu_info()[1].model),
    memory_total = Int(Sys.total_memory()),
    memory_free = Int(Sys.free_memory()),
)

usable(g) = g.functional

const LEVELS = (:ok, :info, :warn, :fail)
max_level(a, b) = findfirst(==(a), LEVELS) >= findfirst(==(b), LEVELS) ? a : b


## Advice

function general_advice(sys, gpus)
    advice = Tuple{Symbol,String}[]
    if sys.threads == 1 && sys.cpu_threads > 1
        push!(advice, (:warn, "Julia runs on 1 thread out of $(sys.cpu_threads): start it with `julia -t auto` to use the `:threads` engine."))
    end
    for g in gpus
        if !g.loaded && g.hardware
            hw = g.package == "CUDA" ? "An NVIDIA GPU" : "An Intel GPU"
            push!(advice, (:info, "$hw seems present: load $(g.package) (`using $(g.package)`) to use it."))
        elseif g.loaded && !g.functional
            push!(advice, (:warn, "$(g.package) is loaded but not functional: $(g.reason)"))
        end
    end
    return advice
end

fits(bytes, available) = bytes <= 0.9*available ? :ok : bytes <= available ? :warn : :fail

ratio(r::Rational) = "$(numerator(r))/$(denominator(r))"

# memory available on each device: (label, bytes, gpu or nothing)
function devices(sys, gpus)
    ds = Any[("CPU RAM", sys.memory_free, nothing)]
    for g in filter(usable, gpus)
        push!(ds, (g.package, something(g.memory_free, g.memory_total), g))
    end
    return ds
end

function model_advice(sys, gpus, model, precisions)
    advice = Tuple{Symbol,String}[]
    T0, m0 = first(precisions)
    if m0 isa Exception
        push!(advice, (:fail, "The workspace of the solvers could not be built: " * short_error(m0)))
        return advice
    end
    if !(m0.ε isa Number) || !isfinite(m0.ε)
        push!(advice, (:fail, "The residuals at the initial guess are not finite (or fail): " * (m0.ε isa Exception ? short_error(m0.ε) : string(m0.ε))))
    end
    ok_gpus = filter(usable, gpus)
    if !isempty(ok_gpus)
        model.isbits || push!(advice, (:fail, "The discretized model is not `isbits`: it cannot be passed to GPU kernels."))
        if model.N < 10_000
            push!(advice, (:info, "With $(model.N) grid points, a GPU is unlikely to be faster than CPU threads (break-even measured around 10⁴ points)."))
        end
    end
    # Float32
    i32 = findfirst(p -> p[1] == Float32, precisions)
    f32 = i32 === nothing ? nothing : precisions[i32][2]
    if f32 isa Exception
        push!(advice, (:warn, "No Float32 copy of the model (`Dolo.convert_precision(Float32, dmodel)` fails): define its methods for `hypeof(model)` (see examples/ymodels/rbc_mc.jl). Error: " * short_error(f32)))
    elseif f32 !== nothing
        for (fun, T) in pairs(f32.results)
            T == Float32 || push!(advice, (:warn, "In the Float32 copy, `$fun` returns $T values: a Float64 literal promotes the computation (write e.g. `T(0.95)` with the precision `T` of the model). GPUs without Float64 cannot run it."))
        end
    end
    for g in ok_gpus
        dev = g.device
        if T0 == Float64 && !g.float64
            push!(advice, (:warn, "$dev has no Float64 support: solve a Float32 copy of the model on it (`Dolo.convert_precision(Float32, dmodel)`, with `tol_ε ≈ 1e-4`)."))
        elseif T0 == Float64 && g.float64_ratio !== nothing && g.float64_ratio <= 1//8 && f32 isa NamedTuple
            push!(advice, (:info, "On $dev, Float64 runs at about $(ratio(g.float64_ratio)) of the speed of Float32: a Float32 copy of the model will be much faster on it (`Dolo.convert_precision(Float32, dmodel)`, with `tol_ε ≈ 1e-4`)."))
        end
        # memory, in each precision the GPU can run
        available = something(g.memory_free, g.memory_total)
        runnable = [(T, r) for (T, r) in precisions if r isa NamedTuple && (T == Float32 || g.float64)]
        summaries = String[]
        worst = :ok
        for (T, r) in runnable
            st = [fits(b, available) for (_, b) in r.solvers]
            ok = [name for ((name, _), s) in zip(r.solvers, st) if s != :fail]
            if all(==(:ok), st)
                continue
            elseif isempty(ok)
                push!(summaries, "in $T no solver fits ($(human(r.solvers[1][2])) for time iteration)"); worst = :fail
            elseif st[3] != :fail
                push!(summaries, "in $T GMRES is tight ($(human(r.solvers[3][2])))"); worst = max_level(worst, :warn)
            elseif st[2] != :fail
                push!(summaries, "in $T GMRES does not fit ($(human(r.solvers[3][2]))): use `linsolve=:neumann` ($(human(r.solvers[2][2])))"); worst = max_level(worst, :warn)
            else
                push!(summaries, "in $T only plain time iteration fits (Newton needs $(human(r.solvers[2][2])))"); worst = max_level(worst, :warn)
            end
        end
        if !isempty(summaries)
            push!(advice, (worst, "Memory on $dev ($(human(available)) available): " * join(summaries, "; ") * "."))
        end
        # kernel arguments
        # in the precision the GPU would run (Float32 if it has no Float64)
        mk = (T0 == Float64 && !g.float64 && f32 isa NamedTuple) ? f32 : m0
        ka = maximum(skipmissing(values(mk.kernel_args)); init=0)
        if g.param_limit !== nothing && ka > g.param_limit
            push!(advice, (:fail, "The kernel arguments ($ka bytes) exceed the limit of $dev ($(g.param_limit) bytes)."))
        elseif g.param_limit === nothing && ka > 2048
            push!(advice, (:warn, "The kernel arguments ($ka bytes) may exceed the limit of $dev (often 2 KB on Intel GPUs): if kernels fail to launch, this is a likely cause."))
        end
    end
    return advice
end

# when there is no recommendation
nothing_to_report(advice) = isempty(advice) ? [(:ok, "Nothing to report.")] : advice


## Printing

# an error in one line, without the (long) types of the arguments
function short_error(e)
    e isa MethodError && return "no method `$(nameof(e.f))` for these argument types"
    msg = first(split(sprint(showerror, e), '\n'))
    return length(msg) > 160 ? first(msg, 160) * "…" : msg
end

function human(b)
    b === missing && return "?"
    b < 1e3 && return "$(round(Int, b)) B"
    b < 1e6 && return @sprintf("%.1f KB", b/1e3)
    b < 1e9 && return @sprintf("%.1f MB", b/1e6)
    return @sprintf("%.2f GB", b/1e9)
end

const MARKS = Dict(:ok => ("✓", :green), :warn => ("!", :yellow), :fail => ("✗", :red), :info => ("•", :cyan))

mark(io, level) = (m = MARKS[level]; printstyled(io, m[1]; color=m[2], bold=true))

function header(io, title)
    println(io)
    printstyled(io, "━━ ", title, " ", repeat("━", max(3, 66 - length(title))); bold=true, color=:blue)
    println(io)
end

entry(io, key, value) = (printstyled(io, "  ", rpad(key, 20); bold=true); println(io, value))

function Base.show(io::IO, ::MIME"text/plain", d::Diagnosis)

    s = d.system
    header(io, "System")
    entry(io, "Julia", "$(s.julia), $(s.threads) thread$(s.threads > 1 ? "s" : "") ($(s.cpu_threads) available)")
    entry(io, "CPU", s.cpu)
    entry(io, "Memory", "$(human(s.memory_total)), $(human(s.memory_free)) available")
    engines = [":serial", ":threads ($(s.threads))"]
    for g in filter(usable, d.gpus)
        push!(engines, ":ka ($(g.package))")
    end
    entry(io, "Engines", join(engines, ", "))

    header(io, "GPUs")
    for g in d.gpus
        print(io, "  ")
        if !g.loaded
            mark(io, g.hardware ? :info : :fail)
            print(io, " ", rpad(g.package, 8))
            println(io, "not loaded", g.hardware ? " (a device seems present: `using $(g.package)`)" : "")
        elseif !g.functional
            mark(io, :fail); print(io, " ", rpad(g.package, 8)); println(io, "not functional: ", g.reason)
        else
            mark(io, :ok); print(io, " ", rpad(g.package, 8)); printstyled(io, g.device; bold=true); println(io, " (", g.details, ")")
            mem = g.memory_free === nothing ? human(g.memory_total) : "$(human(g.memory_total)), $(human(g.memory_free)) free"
            println(io, " "^13, "memory ", mem, g.integrated ? " (shared with the system)" : "")
            f64 = !g.float64 ? "no Float64" :
                g.float64_ratio === nothing ? "Float64 supported" :
                "Float64 at ≈$(ratio(g.float64_ratio)) of Float32 throughput"
            println(io, " "^13, f64, g.param_limit === nothing ? "" : ", kernel arguments ≤ $(g.param_limit) bytes")
        end
    end

    if d.model !== nothing
        m = d.model
        header(io, "Model: $(m.name)")
        exo = filter(v -> v in m.exogenous, m.states)
        entry(io, "States", join(m.states, ", ") * (isempty(exo) ? "" : "  (exogenous: $(join(exo, ", ")))"))
        entry(io, "Controls", "$(length(m.controls)): " * join(m.controls, ", "))
        entry(io, "Exogenous process", m.process)
        T0, m0 = first(d.precisions)
        nodes = m0 isa Exception ? "" : ", $(m0.future_nodes) future states per point"
        entry(io, "Grid", "$(m.n_exo) exogenous node$(m.n_exo > 1 ? "s" : "") × $(join(m.endo, "×")) = $(m.N) points$nodes")
        entry(io, "Precision", string(T0))

        header(io, "Memory of the solvers")
        devs = devices(d.system, d.gpus)
        for (T, r) in d.precisions
            if r isa Exception
                print(io, "  "); mark(io, :fail); println(io, " $T: ", short_error(r))
                continue
            end
            printstyled(io, "  ", rpad("$T", 38), lpad("needed", 10); bold=true)
            for (label, avail, _) in devs
                printstyled(io, lpad(label, 10); bold=true)
            end
            println(io)
            for (name, b) in r.solvers
                print(io, "  ", rpad(name, 38), lpad(human(b), 10))
                for (label, avail, g) in devs
                    if g !== nothing && T == Float64 && !g.float64
                        printstyled(io, lpad("no F64", 10); color=:red)
                    else
                        st = fits(b, avail)
                        print(io, " "^9); mark(io, st)
                    end
                end
                println(io)
            end
            ka = r.kernel_args
            println(io, "  ", rpad("kernel arguments", 38), lpad(human(ka.residuals), 10), "  (residuals; operator L: $(human(ka.operator)))")
            println(io)
        end
        avail = join(("$label $(human(a))" for (label, a, _) in devs), ", ")
        println(io, "  available: ", avail, "  (", MARKS[:ok][1], " < 90%, ", MARKS[:warn][1], " < 100%)")

        header(io, "Checks")
        T0, m0 = first(d.precisions)
        check(level, text) = (print(io, "  "); mark(io, level); println(io, " ", text))
        check(m.isbits ? :ok : :fail, "discretized model is isbits (required by GPU kernels)")
        if !(m0 isa Exception)
            ε = m0.ε
            check(ε isa Number && isfinite(ε) ? :ok : :fail, "residuals at the initial guess: " * (ε isa Number ? @sprintf("max |F| = %.2e", ε) : short_error(ε)))
        end
        i32 = findfirst(p -> p[1] == Float32, d.precisions)
        if i32 !== nothing
            r32 = d.precisions[i32][2]
            if r32 isa Exception
                check(:warn, "Float32 copy of the model: not available ($(short_error(r32)))")
            else
                bad = [string(k) for (k, T) in pairs(r32.results) if T != Float32]
                check(isempty(bad) ? :ok : :warn, "Float32 copy of the model" * (isempty(bad) ? ": all computations in Float32" : ": Float64 in $(join(bad, ", "))"))
            end
        end
    end

    header(io, "Recommendations")
    for (level, text) in d.advice
        print(io, "  "); mark(io, level); println(io, " ", text)
    end
end

Base.show(io::IO, d::Diagnosis) = print(io, "Dolo.Diagnosis (display it for the report)")
