

"""
Computes the residuals of the arbitrage equations. The general form of the arbitrage equation is

    `0 = E_t [f(m, s, x, M, S, X; p)]`

where `m` are current exogenous variables, `s` are current states,
`x` are current controls, `M` are next period's exogenous variables, `S` are next period's states, `X` are next period's controls, and `p` are the model parameters. This function evaluates the right hand side of the arbitrage equation for the given inputs.

If the list of current controls `x` is provided as a two-dimensional array (`ListOfPoints`), it is transformed to a one-dimensional array (`ListOfListOfPoints`).


# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dprocess`: Discretized exogenous process.
* `s::ListOfPoints`: List of state variable values.
* `x::ListOfListOfPoints`: List of control variable values associated with each exogenous shock.
* `p::Vector{Float64}`: Model parameters.
* `dr`: Current guess for the decision rule.
# Returns
* `res`: Residuals of the arbitrage equation associated with each exogenous shock.
"""
function euler_residuals_ti(model, dprocess::AbstractDiscretizedProcess,s::ListOfPoints{d}, x::Vector{ListOfPoints{n_x}}, p::SVector, dr) where n_x where d

    N = length(s)
    # TODO: allocate properly...
    res = deepcopy(x)
    for i_m=1:length(res)
        res[i_m][:] *= 0.0
    end
    for i in 1:size(res, 1)
        m = node(Point, dprocess, i)
        for (w, M, j) in get_integration_nodes(Point, dprocess,i)
            # Update the states
            S = Dolo.transition(model, m, s, x[i], M, p)::ListOfPoints{d}
            X = dr(i, j, S)::ListOfPoints{n_x}
            res[i] += w*Dolo.arbitrage(model, m, s, x[i], M, S, X, p)::ListOfPoints{n_x}
        end
    end
    return res::Vector{ListOfPoints{n_x}}
end



immutable TimeIterationLog
    header::Array{String, 1}
    keywords::Array{Symbol, 1}
    entries::Array{Any, 1}
end

function TimeIterationLog()
    header = [ "It", "ϵₙ", "ηₙ=|xₙ-xₙ₋₁|", "λₙ=ηₙ/ηₙ₋₁", "Time", "Newton steps"]
    keywords = [:it, :err, :gain, :time, :nit]
    TimeIterationLog(header, keywords, [])
end

function append!(log::TimeIterationLog; verbose=true, entry...)
    d = Dict(entry)
    push!(log.entries, d)
    verbose && show_entry(log, d)
end

function initialize(log::TimeIterationLog; verbose=true)
    verbose && show_start(log)
end

function finalize(log::TimeIterationLog; verbose=true)
    verbose && show_end(log)
end

function show(log::TimeIterationLog)
    show_start(log)
    for entry in log.entries
        show_entry(log,entry)
    end
    show_end(log)
end

function show_start(log::TimeIterationLog)
    println(repeat("-", 66))
    @printf "%-6s%-16s%-16s%-16s%-16s%-5s\n" "It" "ϵₙ" "ηₙ=|xₙ-xₙ₋₁|" "λₙ=ηₙ/ηₙ₋₁" "Time" "Newton steps"
    println(repeat("-", 66))
end


function show_entry(log::TimeIterationLog, entry)
    it = entry[:it]
    epsilon = entry[:epsilon]
    err = entry[:err]
    gain = entry[:gain]
    time = entry[:time]
    nit = entry[:nit]
    @printf "%-6i%-16.2e%-16.2e%-16.2e%-16.2e%-5i\n" it epsilon err gain time nit
end

function show_end(log::TimeIterationLog)
    println(repeat("-", 66))
end

###


type IterationTrace
    trace::Array{Any,1}
end


type TimeIterationResult
    dr::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    dprocess::AbstractDiscretizedProcess
    x_converged::Bool
    x_tol::Float64
    err::Float64
    log::TimeIterationLog
    trace::Union{Void,IterationTrace}
end

converged(r::TimeIterationResult) = r.x_converged
function Base.show(io::IO, r::TimeIterationResult)
    @printf io "Results of Time Iteration Algorithm\n"
    @printf io " * Complementarities: %s\n" string(r.complementarities)
    @printf io " * Discretized Process type: %s\n" string(typeof(r.dprocess))
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * Convergence: %s\n" converged(r)
    @printf io "   * |x - x'| < %.1e: %s\n" r.x_tol r.x_converged
end



"""
Computes a global solution for a model via backward time iteration. The time iteration is applied to the residuals of the arbitrage equations.

If the initial guess for the decision rule is not explicitly provided, the initial guess is provided by `ConstantDecisionRule`.
If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, `model.exogenous`

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `process`: The stochastic process associated with the exogenous variables in the model.
* `init_dr`: Initial guess for the decision rule.
# Returns
* `dr`: Solved decision rule.
"""
function time_iteration(model::Model, dprocess::AbstractDiscretizedProcess,
                        grid, init_dr;
                        verbose::Bool=true,
                        maxit::Int=500, tol_η::Float64=1e-7, trace::Bool=false, maxbsteps=5, dampen=1.0,
                        solver=Dict(), complementarities=true)

    if get(solver, :type, :__missing__) == :direct
        return time_iteration_direct(
            model, dprocess, grid, init_dr; verbose=verbose,
            maxit=maxit, tol_η=tol_η
        )
    end



    endo_nodes = nodes(ListOfPoints, grid)
    N = n_nodes(grid)
    n_s_endo = size(endo_nodes, 2)
    n_s_exo = n_nodes(dprocess)

    # initial guess
    # number of smooth decision rules
    n_m = max(n_s_exo, 1)
    p = SVector(model.calibration[:parameters]...)

    x0 = [init_dr(i, endo_nodes) for i=1:n_m]

    ti_trace = trace ? IterationTrace([x0]) : nothing

    n_x = length(model.calibration[:controls])
    if complementarities == true
        x_lb = [controls_lb(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m]
        x_ub = [controls_ub(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m]
        BIG = 100000
        for i=1:n_m
          for n=1:N
            x_lb[i][n] = max.(x_lb[i][n],-BIG)
            x_ub[i][n] = min.(x_ub[i][n], BIG)
          end
        end
    end

    # create decision rule (which interpolates x0)
    dr = CachedDecisionRule(dprocess, grid, x0)

    steps = 0.5.^collect(0:maxbsteps)


    err_0 = NaN
    err = 1.0

    log = TimeIterationLog()
    initialize(log, verbose=verbose)

    it = 0
    while it<maxit && err>tol_η

        it += 1
        tic()
        set_values!(dr, x0)
        fobj(u) = euler_residuals_ti(model, dprocess, endo_nodes, u, p, dr)

        tt = euler_residuals_ti(model, dprocess, endo_nodes, x0, p, dr)

        if complementarities
            res = newton(fobj, x0, x_lb, x_ub)
        else
            res = newton(fobj, x0)
        end

        x1 = res.solution
        nit = res.iterations
        epsil = res.errors[1]


        trace && push!(ti_trace.trace, x1)

        err = maxabs(x1-x0)
        x0=x1

        gain = err / err_0
        err_0 = err

        elapsed = toq()

        append!(log; verbose=verbose, it=it, epsilon=epsil, err=err, gain=gain, time=elapsed, nit=nit)
    end

    finalize(log, verbose=verbose)

    # TODO: somehow after defining `fobj` the `dr` object gets `Core.Box`ed
    #       making the return type right here non-inferrable.

    converged = err < tol_η
    res = TimeIterationResult(dr.dr, it, true, dprocess, converged, tol_η, err, log, ti_trace)

end

# get grid for endogenous
function time_iteration(model, dprocess, dr; grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    return time_iteration(model, dprocess, grid, dr;  kwargs...)
end

# get stupid initial rule
function time_iteration(model, dprocess::AbstractDiscretizedProcess; grid=Dict(), kwargs...)

    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration(model, dprocess, init_dr;  grid=grid, kwargs...)
end


function time_iteration(model, init_dr; grid=Dict(), kwargs...)
    dprocess = discretize( model.exogenous )
    return time_iteration(model, dprocess, init_dr; grid=grid, kwargs...)
end


function time_iteration(model; grid=Dict(), kwargs...)
    dprocess = discretize( model.exogenous )
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration(model, dprocess, init_dr; grid=grid, kwargs...)
end
