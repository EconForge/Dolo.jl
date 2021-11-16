

function new_time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    maxit=500
)
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([x0]) : nothing


    z0 = deepcopy(F.x0)
    z1 = deepcopy(z0)

    log = TimeIterationLog()
    initialize(log, verbose=verbose)

    local err_η

    err_η_0 = NaN

    it = 0
    while it<=maxit

        it += 1

        t1 = time_ns()

        res = F(z0, z0, true)

        err_ε= norm(res)
        if err_ε<tol_ε
            break
        end

        sol = newton(u-> F(u,z0,false), z0)

        z1 = sol.solution

        trace && push!(ti_trace.trace, z1)


        δ = z1 - z0

        err_η = norm(δ)
        gain = err_η_0 / err_η
        err_η_0 = err_η

        if err_η<tol_η
            break 
        end
        z0 = z1

        elapsed = time_ns() - t1

        append!(log; verbose=verbose, it=it, epsilon=err_ε, err=err_η, gain=gain, time=elapsed, nit=sol.iterations)

    end

    finalize(log, verbose=verbose)


    res = TimeIterationResult(F.dr.dr, it, complementarities, F.dprocess, err_η<tol_η, tol_η, err_η, log, ti_trace)

    return res

end

function loop_ti(model::AbstractModel{V}; kwargs...) where V
    F = Euler(model)
    z0 = F.x0
    return loop_ti(F, z0; kwargs...)
end


function loop_iti(model::AbstractModel{V}; kwargs...) where V
    F = Euler(model)
    z0 = F.x0
    return loop_iti(F, z0; kwargs...)
end


function loop_ti(F::Euler, z0::MSM{V}; tol_η=1e-7, tol_ε=1e-8, T=100) where V

    z1 = deepcopy(z0)

    local err
    for t=1:T

        res = F(z0, z0, true)
        err_1 = norm(res)
        if err_1<tol_ε
            break
        end

        sol = newton(u-> F(u,z0,false), z0)

        z1 = sol.solution


        δ = z1 - z0

        err = norm(δ)

        if err<tol_η
            break 
        end
        z0 = z1

    end

    return err
end

function improved_time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    tol_ν = 1e-10,
    maxit=500,
    smaxit=500
)
# function loop_iti(F::Euler, z0::MSM{V}; verbose=true, tol_η=1e-7, tol_ε=1e-8, tol_κ=1e-8, T=500, K=500, switch=5, mode=:iti) where V
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([x0]) : nothing


    z0 = deepcopy(F.x0)

    local err_η, err_ε

    log = TimeIterationLog(
        ["It", "ϵₙ", "ηₙ=|xₙ-xₙ₋₁|", "Time"],
        [:it, :err, :sa, :time],
        []
    )
    initialize(log; verbose=verbose)

 
    it = 0
    while it<=maxit

        t1 = time_ns()

        it += 1

        r = F(z0, z0, true)

        err_ε = norm(r)



        J = df_A(F, z0, z0 ; set_future=false)

        L = df_B(F, z0, z0 ; set_future=false)

        mult!(L, -1.0)
        prediv!(L, J) # J\L

        π = -J\r

        count = 0

        u = π
        δ = π
        for i=1:smaxit
            count +=1
            u = L*u
            δ += u
            if norm(u)<tol_ν
                break
            end
        end
        

        z0 = z0 + δ

        err_η = norm(δ)

        elapsed = time_ns() - t1

        append!(log; verbose=verbose,
            it = it, 
            err = err_ε,
            sa = err_η,
            time = elapsed
        )

        if err_η<tol_η
            break 
        end

    end

    finalize(log, verbose=verbose)

    res = ImprovedTimeIterationResult(
        F.dr.dr,
        it,
        err_ε,
        err_η,
        err_η<tol_η,
        complementarities,
        F.dprocess,
        tol_ν,
        NaN,
        0,
        0,
        0,
        ti_trace
    )
        
    res
end
