using Optim


struct ValueIterationLog
    header::Array{String, 1}
    keywords::Array{Symbol, 1}
    entries::Array{Any, 1}
end

function ValueIterationLog()
    header = [ "It", "ηₙ=|xₙ-xₙ₋₁|", "νₙ=|vₙ-vₙ₋₁|", "Time", "Eval steps"]
    keywords = [:it, :sa_x, :sa_v, :time, :nit]
    ValueIterationLog(header, keywords, [])
end


mutable struct ValueIterationResult
    dr::AbstractDecisionRule
    drv::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    dprocess::AbstractDiscretizedProcess
    x_converged::Bool
    x_tol::Float64
    x_err::Float64
    v_converged::Bool
    v_tol::Float64
    v_err::Float64
    log::ValueIterationLog
    trace::Union{Nothing,IterationTrace}
end

converged(r::ValueIterationResult) = r.x_converged && r.v_converged
function Base.show(io::IO, r::ValueIterationResult)
    @printf io "Results of Value Iteration Algorithm\n"
    @printf io " * Complementarities: %s\n" string(r.complementarities)
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Discretized Process type: %s\n" string(typeof(r.dprocess))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * Convergence: %s\n" converged(r)
    @printf io "   * |x - x'| < %.1e: %s\n" r.x_tol r.x_converged
    @printf io "   * |v - v'| < %.1e: %s\n" r.v_tol r.v_converged

end

function append!(log::ValueIterationLog; verbose=true, entry...)
    d = Dict(entry)
    push!(log.entries, d)
    verbose && show_entry(log, d)
end

function initialize(log::ValueIterationLog; verbose=true)
    verbose && show_start(log)
end

function finalize(log::ValueIterationLog; verbose=true)
    verbose && show_end(log)
end

function show(log::ValueIterationLog)
    show_start(log)
    for entry in log.entries
        show_entry(log,entry)
    end
    show_end(log)
end

function show_start(log::ValueIterationLog)
    println(repeat("-", 66))
    @printf "%-6s%-16s%-16s%-16s%-5s\n" "It" "|xₙ-xₙ₋₁|" "|vₙ-vₙ₋₁|" "Time" "Eval steps"
    println(repeat("-", 66))
end


function show_entry(log::ValueIterationLog, entry)
    it = entry[:it]
    sa_x = entry[:sa_x]
    sa_v = entry[:sa_v]
    time = entry[:time]
    nit = entry[:nit]
    @printf "%-6i%-16.2e%-16.2e%-16.2e%-5i\n" it sa_x sa_v time nit
end

function show_end(log::ValueIterationLog)
    println(repeat("-", 66))
end


"""
Evaluate the value function under the given decision rule.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dr`: Current guess for the decision rule.
# Returns
* `drv`: Value function.
"""
function evaluate_policy(model, dprocess::AbstractDiscretizedProcess, grid, dr;
                            verbose::Bool=true, maxit::Int=5000, tol=1e-6)


    # extract parameters
    p = SVector(model.calibration[:parameters]...)
    β = model.calibration.flat[:beta]

    # states today are the grid
    s = nodes(ListOfPoints, grid)

    # Number of endogenous nodes
    N = length(s)

    # number of smooth decision rules
    nsd = max(n_nodes(dprocess), 1)

    n_u = length(model.calibration[:rewards])

    res = [zeros(Point{n_u},N) for i=1:nsd]

    # Value function : v_t = u + β*E[V_{t+1}]
    u = deepcopy(res)
    v = deepcopy(res)
    v0 = deepcopy(res)

    # Expected utility
    E_V = deepcopy(res)

    # evaluate u at time t. This is constant because dr is constant within this
    # function
    x = [dr(i, s) for i=1:nsd]
    for i = 1:size(res, 1)
        m = node(Point, dprocess, i)
        u[i][:] = Dolo.felicity(model, m, s, x[i], p)
    end

    # Initial guess for the value function
    v0 = deepcopy(u)

    #Preparation for a loop
    err = 10.0
    it = 0

    drv = CachedDecisionRule(dprocess, grid, v0)

    verbose && @printf "%-6s%-12s\n" "It" "SA"
    verbose && println(repeat("-", 14))

    while err > tol && it < maxit
        it += 1
        # Interpolate v0 on the grid
        # Compute value function for each smooth decision rule
        for i = 1:nsd
            m = node(Point, dprocess, i)
            for (w, MM, j) in get_integration_nodes(dprocess,i)
                    M = SVector(MM...)
                     # Update the states
                     S = Dolo.transition(model, m, s, x[i], M, p)
                     # Update value function
                      E_V[i] += w*drv(i, j, S)
            end
        end


        err = 0.0
        for i in 1:length(v)
            v[i][:] = u[i] + β*E_V[i]
            err = max(err, maxabs(v[i] - v0[i]))
            copyto!(v0[i], v[i])
            E_V[i] *= 0.0
        end

        verbose && @printf "%-6i%-12.2e\n" it err

        set_values!(drv, v0)
    end

    return (drv, it)
end

function evaluate_policy(model, dr; grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    dprocess = discretize(model.exogenous)
    return evaluate_policy(model, dprocess, grid, dr; kwargs...)[1]

end

"""
Evaluate the right hand side of the value function at given values of states, controls, and exogenous variables.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `β::Float64`: Value of the discount factor.
* `dprocess::`: Discretized exogenous process.
* `drv`: Current guess for the value function object.
* `i::Int64`: Index of node for exogeous variable(s).
* `s::ListOfPoints`: List of state variables.
* `x0::ListOfListOfPoints`: List of control variables.
* `p::Vector{Float64}`: Model parameters.
# Returns
* `E_V::`: Right hand side of the value function.
"""
function update_value(model, β::Float64, dprocess, drv, i, s::Point{d},
                      x0::Point{n_x}, p::Point{n_p}) where d where n_x where n_p
    m = node(Point,dprocess, i)
    N = length(s)
    E_V = Point{1}(0.0)
    for (w, M, j) in get_integration_nodes(Point,dprocess,i)
        # Update the states
        S = Dolo.transition(model, m, s, x0, M, p)
        E_V += w*drv(i, j, S)[1]
    end
    u = Dolo.felicity(model, m, s, x0, p)[1]
    V = u + β*E_V
    return V
end



# function update_value(model, β::Float64, dprocess, drv, i, s::Point{d},
#                       x0::Float64, p::Point{n_p}) where d where n_p
#     update_value(model, β, dprocess, drv, i, s, SVector{1,Float64}(x0), p)
# end
#


function call_optim(fobj, initial_x, lower, upper, optim_opts)
    if length(initial_x) == 1
        return optimize(fobj, lower[1], upper[1])
    end
    return optimize(fobj, lower, upper, initial_x)
end

"""
Solve for the value function and associated decision rule using value function iteration.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `pdr`: Initial guess for the decision rule.
# Returns
* `dr`: Solved decision rule object.
* `drv`: Solved value function object.
"""
function value_iteration(
        model, dprocess::AbstractDiscretizedProcess, grid, init_dr;
        discount_symbol=:beta,
        maxit::Int=1000, tol_x::Float64=1e-8, tol_v::Float64=1e-8,
        optim_options=Dict(), eval_options=Dict(),
        verbose::Bool=true, trace::Bool=false
    )


    β = model.calibration.flat[discount_symbol]

    # compute the value function

    p = SVector(model.calibration[:parameters]...)

    endo_nodes = nodes(ListOfPoints,grid)

    # Number of endogenous nodes
    N = length(endo_nodes)

    # number of smooth decision rules
    n_m = nsd = max(n_nodes(dprocess), 1)

    res = [zeros(Point{1},N) for i=1:nsd]

    # States at time t+1
    S = copy(endo_nodes)

    # Controls at time t
    x = [init_dr(i, endo_nodes) for i=1:nsd]
    x0 = deepcopy(x)

    dr = CachedDecisionRule(dprocess, grid, x0)

    drv = dr # never used. Just to escape loop scope
    v0 = [zeros(Point{1},N) for i=1:nsd]
    v = deepcopy(v0)

    ti_trace = trace ? IterationTrace([x0, v0]) : nothing

    # fourth time I repeat that one...
    n_x = length(model.calibration[:controls])
    x_lb = [controls_lb(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m]
    x_ub = [controls_ub(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m]
    BIG = 100000
    for i=1:n_m
      for n=1:N
        x_lb[i][n] = max.(x_lb[i][n],-BIG)
        x_ub[i][n] = min.(x_ub[i][n], BIG)
      end
    end

    #Preparation for a loop
    err_v = 10.0
    err_x = 10.0
    err_eval = 10.0

    it = 0
    it_eval = 0

    maxit_eval = get(eval_options, :maxit, 1000)
    tol_eval = get(eval_options, :tol, 1e-8)

    optim_opts = Optim.Options(optim_options...)

    log = ValueIterationLog()
    initialize(log, verbose=verbose)

    mode = :improve
    done = false

    while !done

        tic()
        it += 1

        drv, it_eval = evaluate_policy(model, dprocess, grid, dr; verbose=false, eval_options...)
        mode = :improve

        # else
        for i = 1:size(res, 1)
            m = node(dprocess, i)
            for n = 1:N
                s = endo_nodes[n]
                # optimize vals
                fobj(u) = -update_value(model, β, dprocess, drv, i, s, SVector(u...), p)[1]*1000
                lower = Float64[x_lb[i][n]...]
                upper = Float64[x_ub[i][n]...]
                upper = clamp.(upper, -Inf, 1000000)
                lower = clamp.(lower, -1000000, Inf)
                initial_x = Float64[x0[i][n]...] #, :]
                test = fobj(initial_x)
                results = call_optim(fobj, initial_x, lower, upper, optim_opts)
                xn = Optim.minimizer(results)
                nv = -Optim.minimum(results)/1000.0
                x[i][n] = SVector(xn...)
                v[i][n] = SVector(nv)
            end
        end
        # compute diff in values
        err_v = maxabs(v-v0)
        err_x = maxabs(x-x0)
        for i in 1:nsd
            copyto!(v0[i], v[i])
            copyto!(x0[i], x[i])
        end

        # update values and policies
        set_values!(drv, v0)
        set_values!(dr, x0)

        # terminate only if policy didn't move
        done = (err_x<tol_x) || (err_v<tol_v) || (it>=maxit)

        trace && push!(ti_trace.trace, [x0, v0]) # this is ridiculous since only one of them is updated !

        elapsed = toq()

        append!(log; verbose=verbose, it=it, sa_x=err_x, sa_v=err_v, time=elapsed, nit=it_eval)

    end

    finalize(log, verbose=verbose)

    converged_x = err_x<tol_x
    converged_v = err_v<tol_v

    ValueIterationResult(dr.dr, drv.dr, it, true, dprocess, converged_x, tol_x, err_x, converged_v, tol_v, err_v, log, ti_trace)

end

# get stupid initial rule
function value_iteration(model, dprocess::AbstractDiscretizedProcess, init_dr; grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    return value_iteration(model, dprocess, grid, init_dr;  kwargs...)
end

# get stupid initial rule
function value_iteration(model, dprocess::AbstractDiscretizedProcess; grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return value_iteration(model, dprocess, grid, init_dr;  kwargs...)
end


function value_iteration(model, init_dr; grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    dprocess = discretize( model.exogenous )
    return value_iteration(model, dprocess, grid, init_dr; kwargs...)
end


function value_iteration(model; grid=Dict(), kwargs...)
    grid = get_grid(model; options=grid)
    dprocess = discretize( model.exogenous )
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return value_iteration(model, dprocess, grid, init_dr; kwargs...)
end

# compatibility
const solve_policy = value_iteration
