using LinearMaps
using Adapt
using TimerOutputs: TimerOutput, @timeit, disable_timer!

"""
    time_iteration_workspace(dmodel; interp_mode=:linear, accelerate=false, linsolve=:neumann, krylov_dim=50, dest=Array)

Preallocate the memory used by `time_iteration` and `newton`.
The linear operator `L` (and Neumann buffers) are allocated when `accelerate=true`,
the GMRES workspace (of dimension `krylov_dim`) when `linsolve=:gmres`.
"""
function time_iteration_workspace(dmodel; interp_mode=:linear, accelerate=false, linsolve=:neumann, krylov_dim=50, dest=Array)

    x0 = (Dolo.initial_guess(dmodel))
    T = eltype(eltype(x0))

    x1 = deepcopy(x0)
    x2 = deepcopy(x0)
    r0 = deepcopy(x0)
    dx = deepcopy(x0)
    N = length(dx)
    n = length(dx.data[1])
    J = GArray(
        dmodel.grid,
        zeros(SMatrix{n,n,T,n*n}, N)
    )
    vars = variables(dmodel.model.controls)
    φ = DFun(dmodel.model.states, x0, vars; interp_mode=interp_mode)

    tt = (;x0, x1, x2, r0, dx, J, φ)

    if accelerate
        L = Dolo.dF_2(dmodel, x1, φ) # TODO : speed up (keep only the allocation, not the computation)
        memn = (;du=deepcopy(x0), dv=deepcopy(x0))
        tt = merge(tt, (;L, memn))
    end

    if linsolve == :gmres
        gmres = gmres_workspace(T, N*n, krylov_dim)
        tt = merge(tt, (;gmres))
    end

    return adapt(dest, tt)

end

newton_workspace(dmodel::ADModel; kwargs...) = time_iteration_workspace(dmodel; accelerate=true, kwargs...)


function check_workspace(workspace, accelerate, linsolve)
    linsolve in (:neumann, :gmres) || error("Unknown linear solver $linsolve. Use :neumann or :gmres.")
    if accelerate && !haskey(workspace, :L)
        error("Workspace must be created with `accelerate=true`.")
    end
    if accelerate && linsolve == :gmres && !haskey(workspace, :gmres)
        error("Workspace must be created with `linsolve=:gmres`.")
    end
end

is_domain_error(e) = e isa DomainError
is_domain_error(e::TaskFailedException) = is_domain_error(e.task.exception)
is_domain_error(e::CompositeException) = all(is_domain_error, e.exceptions)

# fits φ on x then computes the residuals F(x,φ) in r0
# returns their norm (Inf if the model can't be evaluated)
function residual!(r0, model, x, φ, t_engine, tmr)
    @timeit tmr "fit φ" fit!(φ, x)
    try
        @timeit tmr "eval F" F!(r0, model, x, φ, t_engine)
    catch e
        is_domain_error(e) || rethrow()
        return convert(getprecision(x), Inf)
    end
    return norm(r0)
end

"""
Solve `(I-L) dx = r` with a Neumann series or GMRES. Returns the number of iterations.
"""
function solve_I_minus_L!(dx, L, r, workspace, t_engine; linsolve=:neumann, K=1000, tol=1e-10, tmr=NO_TIMER)
    if linsolve == :neumann
        return neumann!(dx, L, r, t_engine; mem=workspace.memn, K=K, τ_η=tol, tmr=tmr)
    else
        x = ravel(dx)
        b = ravel(r)
        x .= b      # first order approximation
        A! = (w,v) -> (mul!(w, L, v, t_engine); w .= v .- w)
        k, _ = gmres!(x, A!, b, workspace.gmres; tol=tol, tmr=tmr)
        return k
    end
end

const DEPRECATED_KWARGS = (
    improve = :accelerate,
    improve_K = :accelerate_K,
    improve_wait = :accelerate_wait,
    mbsteps = :max_bsteps,
)

function rename_deprecated(kwargs)
    kw = Dict{Symbol,Any}(kwargs)
    for (old, new) in pairs(DEPRECATED_KWARGS)
        if haskey(kw, old)
            Base.depwarn("`$old` is deprecated, use `$new` instead.", :time_iteration)
            kw[new] = pop!(kw, old)
        end
    end
    return kw
end


"""
    time_iteration(model; kwargs...)
    time_iteration(dmodel, workspace=nothing; kwargs...)

Solve the model by time iteration. With `accelerate=true`, after `accelerate_wait`
iterations each time iteration step `x1=T(x0)` is followed by the approximate Newton step
`x0 ← x0 + (I-T')^{-1}(x1-x0)`, where the linear system is solved by a Neumann series
(`linsolve=:neumann`, at most `accelerate_K` terms) or GMRES (`linsolve=:gmres`, Krylov space of
dimension `krylov_dim`). With `safeguarded=true`, the accelerated step is damped until the residual
decreases and replaced by the time iteration step if it doesn't.

Other options:
- `T`: maximum number of iterations
- `tol_ε`: tolerance on the residuals `|F(x,x)|`
- `tol_η`: tolerance on successive time iteration steps `|x1-x0|`
- `newton_steps`, `max_bsteps`: Newton and backtracking steps to compute `T(x0)`
- `engine`: `:serial`, `:threads` or `:ka` (KernelAbstractions, e.g. on GPU), or an `Engine`
- `timer`: record timings (returned in the result)
- `trace`: keep all successive decision rules
"""
function time_iteration(model::YModel; discretization=Dict(), interpolation=:cubic, kwargs...)
    kw = rename_deprecated(kwargs)
    dmodel = discretize(model, discretization)
    wksp = time_iteration_workspace(dmodel;
        interp_mode=interpolation,
        accelerate=get(kw, :accelerate, false),
        linsolve=get(kw, :linsolve, :neumann),
        krylov_dim=get(kw, :krylov_dim, 50)
    )
    time_iteration(dmodel, wksp; kw...)
end

function time_iteration(model::DYModel, workspace=nothing; kwargs...)
    _time_iteration(model, workspace; rename_deprecated(kwargs)...)
end

function _time_iteration(model::DYModel, workspace;
    T=500,
    tol_ε=1e-8,
    tol_η=1e-6,
    newton_steps=10,
    max_bsteps=5,
    accelerate=false,
    accelerate_K=1000,
    accelerate_wait=10,
    linsolve=:neumann,
    linsolve_tol=1e-10,
    krylov_dim=50,
    safeguarded=false,
    verbose=true,
    trace=false,
    engine=:serial,
    timer=false
)

    tmr = TimerOutput()
    timer || disable_timer!(tmr)

    if workspace === nothing
        @timeit tmr "allocation" workspace = time_iteration_workspace(model;
            accelerate=accelerate, linsolve=linsolve, krylov_dim=krylov_dim)
    end
    check_workspace(workspace, accelerate, linsolve)

    log = IterationLog(
            it = ("n", Int),
            err =  ("ϵₙ=|F(xₙ,xₙ)|", Float64),
            sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
            lam = ("λₙ=ηₙ/ηₙ₋₁",Float64),
            elapsed = ("Time", Float64),
            k_newton = ("Newton", Int64),
            k = ("Backsteps", Int64)
        )
    initialize(log, verbose=verbose; message="Time Iteration")

    (;x0, x1, x2, r0, dx, J, φ) = workspace

    Tf = eltype(eltype(x0))
    lam = convert(Tf, 0.5)

    t_engine = select_engine(engine, x0)

    ti_trace = trace ? IterationTrace(typeof(φ)[]) : nothing

    local η_0 = NaN
    local η = NaN
    local k_newton_ = 0
    local k_ = 0
    iterations = 0

    ε = residual!(r0, model, x0, φ, t_engine, tmr)

    for t=1:T

        t1 = time_ns()

        trace && push!(ti_trace.data, deepcopy(φ))

        if ε<tol_ε
            break
        end

        iterations = t

        # solve u->F(u,φ)=0 where φ interpolates x0
        # result in x1

        x1.data .= x0.data

        @timeit tmr "time iteration step" for k_newton=1:newton_steps

            @timeit tmr "eval F" F!(r0, model, x1,  φ, t_engine)
            ε_n = norm(r0)
            if ε_n<tol_ε
                break
            end

            @timeit tmr "compute F_1" dF_1!(J, model, x1,  φ, t_engine)

            # TODO: fix regression or revert to Julia 1.10
            @inline dx.data .= J.data .\ r0.data

            for k=0:max_bsteps
                # fails on GPU
                x2.data .= x1.data .- dx.data .* convert(Tf,lam^k)
                @timeit tmr "eval F" F!(r0, model, x2,  φ, t_engine)
                ε_b = norm(r0)
                if ε_b<ε_n   # false if ε_b is NaN
                    k_newton_ = k_newton
                    k_ = k
                    break
                end
            end

            x1.data .= x2.data

        end

        η = distance(x0, x1)
        gain = η/η_0
        η_0 = η

        if accelerate && t>accelerate_wait

            # x = T(x)
            # x1 = T(x0)
            # x - x1 = T'(x) (x - x0)
            # x = x0 + (I-T')^{-1} (x1 - x0)
            # with T' = L = -F_1^{-1}F_2

            (;L) = workspace

            @timeit tmr "compute L=-F_1\\F_2" begin
                dF_1!(J, model, x1, φ, t_engine)
                dF_2!(L, model, x1, φ, t_engine)
                mul!(L, convert(Tf,-1.0), t_engine)
                ldiv!(J, L, t_engine)
            end

            r0.data .= x1.data .- x0.data

            @timeit tmr "compute (I-L)^(-1).dr" solve_I_minus_L!(dx, L, r0, workspace, t_engine;
                linsolve=linsolve, K=accelerate_K, tol=linsolve_tol, tmr=tmr)

            if safeguarded

                # Armijo-like condition
                c_1 = 1e-4
                α = convert(Tf, 0.5)
                λ = sqrt(1-2*c_1*α)

                accepted = false
                @timeit tmr "safeguard" for i=0:max_bsteps
                    x2.data .= x0.data .+ dx.data .* α^i
                    ε1 = residual!(r0, model, x2, φ, t_engine, tmr)
                    if ε1 < λ*ε
                        x0.data .= x2.data
                        ε = ε1
                        accepted = true
                        break
                    end
                end

                if !accepted
                    # fall back to the time iteration step
                    x0.data .= x1.data
                    ε = residual!(r0, model, x0, φ, t_engine, tmr)
                end

            else
                x0.data .+= dx.data
                ε = residual!(r0, model, x0, φ, t_engine, tmr)
            end

        else
            x0.data .= x1.data
            ε = residual!(r0, model, x0, φ, t_engine, tmr)
        end

        elapsed = time_ns() - t1
        elapsed /= 1000000000

        verbose ? append!(log; verbose=verbose, it=t, err=ε, sa=η, lam=gain, elapsed=elapsed, k_newton=k_newton_, k=k_) : nothing

        if η < tol_η
            break
        end

    end

    verbose ?  finalize(log, verbose=verbose) : nothing

    TimeIterationResult(
        φ,
        iterations,
        tol_ε,
        ε,
        tol_η,
        η,
        log,
        ti_trace,
        timer ? tmr : nothing
    )
end


"""
    newton(model; kwargs...)
    newton(dmodel, workspace=nothing; kwargs...)

Solve the model with a Newton-Krylov method on `x -> F(x, φ(x))`, where `φ` interpolates `x`.

- `linsolve`: `:neumann` or `:gmres` (Krylov space of dimension `krylov_dim`)
- `precondition`: precondition the linear system by `F_1^{-1}` so that it takes the form
  `(I-L) dx = F_1^{-1} r`. Required by `:neumann`. Without preconditioning, GMRES
  typically needs a much larger `krylov_dim`.
- `K`: maximum number of terms in the Neumann series
- `globalization`: `:none` (full Newton steps), `:backtracking` (halve the step until the residual
  decreases, at most `max_bsteps` times) or `:linesearch` (golden section search on the step size,
  `linesearch_maxiter` iterations)
- `T`: maximum number of Newton iterations
- `engine`: `:serial`, `:threads` or `:ka` (KernelAbstractions, e.g. on GPU), or an `Engine`
- `tol_ε`: tolerance on the residuals
"""
function newton(model::YModel; discretization=Dict(), interpolation=:cubic, kwargs...)
    dmodel = discretize(model, discretization)
    wksp = newton_workspace(dmodel;
        interp_mode=interpolation,
        linsolve=get(kwargs, :linsolve, :neumann),
        krylov_dim=get(kwargs, :krylov_dim, 50)
    )
    newton(dmodel, wksp; kwargs...)
end

function newton(model::DYModel, workspace=nothing;
    T=20,
    K=1000,
    tol_ε=1e-8,
    linsolve=:neumann,
    linsolve_tol=1e-10,
    krylov_dim=50,
    precondition=true,
    globalization=:none,
    max_bsteps=10,
    linesearch_maxiter=5,
    verbose=true,
    engine=:serial,
    timer=false
)

    tmr = TimerOutput()
    timer || disable_timer!(tmr)

    if !precondition && linsolve != :gmres
        error("`precondition=false` requires `linsolve=:gmres`.")
    end
    globalization in (:none, :backtracking, :linesearch) || error("Unknown globalization $globalization. Use :none, :backtracking or :linesearch.")

    if workspace === nothing
        @timeit tmr "allocation" workspace = newton_workspace(model; linsolve=linsolve, krylov_dim=krylov_dim)
    end
    check_workspace(workspace, true, linsolve)

    log = IterationLog(
        it = ("n", Int),
        err =  ("ϵₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        elapsed = ("Time", Float64),
        k_lin = ("Linear it.", Int64),
        k = ("Backsteps", Int64)
    )
    initialize(log, verbose=verbose; message="Newton")

    (;x0, x1, dx, r0, J, φ, L) = workspace

    Tf = eltype(eltype(x0))

    t_engine = select_engine(engine, x0)

    iterations = 0
    ε = residual!(r0, model, x0, φ, t_engine, tmr)

    for t=1:T

        if ε<tol_ε
            break
        end

        iterations = t
        t1 = time_ns()

        x1.data .= x0.data

        @timeit tmr "compute Jacobian" begin
            @timeit tmr "compute F_1" dF_1!(J, model, x0, φ, t_engine)
            @timeit tmr "compute F_2" dF_2!(L, model, x0, φ, t_engine)
        end

        if precondition
            # (F_1 + F_2) dx = r  <=>  (I - L) dx = F_1^{-1} r  with L = -F_1^{-1} F_2
            @timeit tmr "compute L=-F_1\\F_2" begin
                mul!(L, convert(Tf,-1.0), t_engine)
                ldiv!(J, L, t_engine)
            end
            r0.data .= J.data .\ r0.data
            k_lin = @timeit tmr "compute (I-L)^(-1).dr" solve_I_minus_L!(dx, L, r0, workspace, t_engine;
                linsolve=linsolve, K=K, tol=linsolve_tol, tmr=tmr)
        else
            wk = workspace.gmres
            x = ravel(dx)
            x .= zero(Tf)
            A! = (w,v) -> (mul!(w, J, v); mul!(wk.aux, L, v, t_engine); w .+= wk.aux)
            k_lin, _ = @timeit tmr "solve (F_1+F_2).dx=r" gmres!(x, A!, ravel(r0), wk; tol=linsolve_tol, tmr=tmr)
        end

        k_ = 0
        if globalization == :none
            x0.data .= x1.data .- dx.data
            ε = residual!(r0, model, x0, φ, t_engine, tmr)
        elseif globalization == :backtracking
            @timeit tmr "backtracking" for i=0:max_bsteps
                k_ = i
                x0.data .= x1.data .- dx.data .* convert(Tf, 0.5^i)
                ε1 = residual!(r0, model, x0, φ, t_engine, tmr)
                if ε1 < ε || i == max_bsteps
                    ε = ε1
                    break
                end
            end
        else
            ε = @timeit tmr "linesearch" golden_section!(x0, x1, dx, r0, model, φ, t_engine, tmr, ε; maxiter=linesearch_maxiter)
        end

        η = distance(x0, x1)

        elapsed = (time_ns() - t1)/1000000000

        verbose ? append!(log; verbose=verbose, it=t, err=ε, sa=η, elapsed=elapsed, k_lin=k_lin, k=k_) : nothing

    end

    verbose ?  finalize(log, verbose=verbose) : nothing

    NewtonResult(
        φ,
        iterations,
        tol_ε,
        ε,
        log,
        timer ? tmr : nothing
    )

end

# Golden section search for the step size α∈[0,1] in x0 = x1 - α dx, minimizing the residual.
# The full Newton step is used if it decreases the residual.
# On exit x0 contains the best point and φ interpolates it. Returns the residual.
function golden_section!(x0, x1, dx, r0, model, φ, t_engine, tmr, ε; maxiter=5)

    Tf = eltype(eltype(x0))

    f = α -> begin
        x0.data .= x1.data .- dx.data .* convert(Tf, α)
        residual!(r0, model, x0, φ, t_engine, tmr)
    end

    ε1 = f(1.0)
    if ε1 < ε
        return ε1
    end

    gr = (sqrt(5.0) - 1.0)/2.0
    a, b = 0.0, 1.0
    c = b - gr*(b-a)
    d = a + gr*(b-a)
    fc = f(c)
    fd = f(d)
    best_α, best_f = fc < fd ? (c, fc) : (d, fd)

    for i=1:maxiter
        if fc < fd
            # minimum is in [a,d]
            b, d, fd = d, c, fc
            c = b - gr*(b-a)
            fc = f(c)
        else
            # minimum is in [c,b]
            a, c, fc = c, d, fd
            d = a + gr*(b-a)
            fd = f(d)
        end
        fc < best_f && ((best_α, best_f) = (c, fc))
        fd < best_f && ((best_α, best_f) = (d, fd))
    end

    return f(best_α)

end
