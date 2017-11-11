using Optim

type ValueIterationResult
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
    trace::Union{Void,IterationTrace}
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
    s = nodes(Point, grid)

    # Number of endogenous nodes
    N = length(s)

    # number of smooth decision rules
    nsd = max(n_nodes(dprocess), 1)

    n_u = length(model.calibration[:rewards])

    res = [to_LOP(zeros(N, n_u)) for i=1:nsd]

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
            copy!(v0[i], v[i])
            E_V[i] *= 0.0
        end

        verbose && @printf "%-6i%-12.2e\n" it err

        set_values!(drv, v0)
    end

    return drv
end

function evaluate_policy(model, dr; grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    dprocess = discretize(model.exogenous)
    return evaluate_policy(model, dprocess, grid, dr; kwargs...)

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
function update_value(model, β::Float64, dprocess, drv, i, s::Vector{Float64},
                      x0::Vector{Float64}, p::Vector{Float64})
    m = node(dprocess, i)
    E_V = 0.0
    for (w, M, j) in get_integration_nodes(dprocess,i)
        # Update the states
        S = Dolo.transition(model, m, s, x0, M, p)
        E_V += w*drv(i, j, S)[1]
    end
    u = Dolo.felicity(model, m, s, x0, p)[1]
    E_V = u + β.*E_V
    return E_V
end
function update_value(model, β::Float64, dprocess, drv, i, s::Vector{Float64},
                      x0::Float64, p::Vector{Float64})
    update_value(model, β, dprocess, drv, i, s, [x0], p)
end


@static if Pkg.installed("Optim") < v"0.9-"
    function call_optim(fobj, initial_x, lower, upper, optim_opts)
        if length(initial_x) == 1
            return optimize(fobj, lower[1], upper[1])
        end
        results = optimize(
            Optim.OnceDifferentiable(fobj), initial_x, lower, upper,
            Fminbox(), optimizer=NelderMead,
            x_tol=1e-10, f_tol=1e-10
        )
    end
else
    function call_optim(fobj, initial_x, lower, upper, optim_opts)
        if length(initial_x) == 1
            return optimize(fobj, lower[1], upper[1])
        end
        results = optimize(
            Optim.OnceDifferentiable(fobj, initial_x), initial_x, lower, upper,
            Fminbox{NelderMead}(), optimizer_o=optim_opts
        )
    end
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
        model, dprocess::AbstractDiscretizedProcess, grid, pdr;
        discount_symbol=:beta,
        maxit::Int=1000, tol_x::Float64=1e-8, tol_v::Float64=1e-8,
        optim_options=Dict(), eval_options=Dict(),
        verbose::Bool=true, trace::Bool=false
    )


    β = model.calibration.flat[discount_symbol]

    dr = CachedDecisionRule(pdr, dprocess)
    # compute the value function
    p = model.calibration[:parameters]

    endo_nodes = nodes(grid)

    # Number of endogenous nodes
    N = size(endo_nodes, 1)

    # number of smooth decision rules
    nsd = max(n_nodes(dprocess), 1)
    res = [zeros(N, 1) for i=1:nsd]

    # bounds on controls
    x_lb = Array{Float64,2}[cat(1, [Dolo.controls_lb(model, node(dprocess, i), endo_nodes[n, :], p)' for n=1:N]...) for i=1:nsd]
    x_ub = Array{Float64,2}[cat(1, [Dolo.controls_ub(model, node(dprocess, i), endo_nodes[n, :], p)' for n=1:N]...) for i=1:nsd]

    # States at time t+1
    S = copy(endo_nodes)

    # Controls at time t
    x = [dr(i, endo_nodes) for i=1:nsd]
    x0 = deepcopy(x)

    dr = CachedDecisionRule(dprocess, grid, x0)

    # this could be integrated in the main loop.
    verbose && print("Evaluating initial policy")

    drv = evaluate_policy(model, dr; verbose=false, eval_options...)

    verbose && print(" (done)\n")

    v0 = [drv(i, endo_nodes) for i=1:nsd]
    v = deepcopy(v0)

    ti_trace = trace ? IterationTrace([x0, v0]) : nothing


    #Preparation for a loop
    err_v = 10.0
    err_x = 10.0
    err_eval = 10.0

    it = 0
    it_eval = 0

    maxit_eval = get(eval_options, :maxit, 1000)
    tol_eval = get(eval_options, :tol, 1e-8)

    optim_opts = Optim.Options(optim_options...)

    mode = :improve
    converged = false

    while !converged

        # it += 1
        if (mode == :eval)

            it_eval = 0
            converged_eval = false
            while !converged_eval
                it_eval += 1
                for i = 1:size(res, 1)
                    m = node(dprocess, i)
                    for n = 1:N
                        s = endo_nodes[n, :]
                        # update vals
                        nv = update_value(model, β, dprocess, drv, i, s, x0[i][n, :], p)
                        v[i][n, 1] = nv
                    end
                end
                # compute diff in values
                err_eval = 0.0
                for i in 1:nsd
                    err_eval = max(err_eval, maximum(abs, v[i] - v0[i]))
                    copy!(v0[i], v[i])
                end
                converged_eval = (it_eval>=maxit_eval) || (err_eval<tol_eval)
                set_values!(drv, v0)
                # if verbose
                #     println("    It: ", it_eval, " ; SA: ", err_eval)
                # end
            end

            mode = :improve

        else

            it += 1
            for i = 1:size(res, 1)
                m = node(dprocess, i)
                for n = 1:N
                    s = endo_nodes[n, :]
                    # optimize vals
                    fobj(u) = -update_value(model, β, dprocess, drv, i, s, u, p)*1000
                    lower = x_lb[i][n, :]
                    upper = x_ub[i][n, :]
                    upper = clamp!(upper, -Inf, 1000000)
                    lower = clamp!(lower, -1000000, Inf)
                    initial_x = x0[i][n, :]
                    results = call_optim(fobj, initial_x, lower, upper, optim_opts)
                    xn = Optim.minimizer(results)
                    nv = -Optim.minimum(results)/1000.0
                    # ii = (fobj(initial_x))
                    # xn, nv = goldensearch(fobj, lower, upper; maxit=1000, tol=1e-10)
                    # jj = (fobj(xn))
                    #
                    # xvec =
                    # println([m,s, lower, upper, initial_x, ii, xn, jj])


                    x[i][n, :] = xn
                    v[i][n, 1] = nv
                end
            end

            # compute diff in values
            err_v = 0.0
            for i in 1:nsd
                err_v = max(err_v, maximum(abs, v[i] - v0[i]))
                copy!(v0[i], v[i])
            end
            # compute diff in policy
            err_x = 0.0
            for i in 1:nsd
                err_x = max(err_x, maximum(abs, x[i] - x0[i]))
                copy!(x0[i], x[i])
            end
            # update values and policies
            set_values!(drv, v0)
            set_values!(dr, x0)

            if verbose
                println("It: ", it, " ; SA: ", err_v, " ; SA_x: ", err_x, " ; (nit) ", it_eval)
            end

            # terminate only if policy didn't move
            converged = ((err_x<tol_x) && (err_v<tol_v)) || (it>=maxit)

            mode = :eval

        end

        trace && push!(ti_trace.trace, [x0, v0]) # this is ridiculous since only one of them is updated !

    end

    converged_x = err_x<tol_x
    converged_v = err_v<tol_v
    ValueIterationResult(dr.dr, drv.dr, it, true, dprocess, converged_x, tol_x, err_x, converged_v, tol_v, err_v, ti_trace)

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
