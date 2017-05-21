type ValueIterationResult
    dr::AbstractDecisionRule
    drv::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    x_converged::Bool
    x_tol::Float64
    x_err::Float64
    v_converged::Bool
    v_tol::Float64
    v_err::Float64
end

converged(r::ValueIterationResult) = r.x_converged && r.v_converged
function Base.show(io::IO, r::ValueIterationResult)
    @printf io "Results of Value Iteration Algorithm\n"
    @printf io " * Complementarities: %s\n" string(r.complementarities)
    @printf io " * Decision Rule type: %s\n" string(typeof(r))
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
function evaluate_policy(model, dr; verbose::Bool=true, maxit::Int=5000)

    # get grid for endogenous
    grid = model.grid

    # obtain discrete exogenous process
    process = model.exogenous
    dprocess = discretize(process)

    # extract parameters
    p = model.calibration[:parameters]
    β=model.calibration.flat[:beta]

    # states today are the grid
    s = nodes(grid)

    # Number of endogenous nodes
    N = size(s, 1)

    # number of smooth decision rules
    nsd = max(n_nodes(dprocess), 1)
    n_u = length(model.calibration[:rewards])
    res = [zeros(N, n_u) for i=1:nsd]

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
        m = node(dprocess, i)
        for n = 1:N
            u[i][n, :] = Dolo.felicity(model, m, s[n, :], x[i][n, :], p)
        end
    end

    # Initial guess for the value function
    v0 = deepcopy(u)

    #Preparation for a loop
    tol = 1e-6
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
            m = node(dprocess, i)
            for j = 1:n_inodes(dprocess, i)
                M = inode(dprocess, i, j)
                w = iweight(dprocess, i, j)
                 for n=1:N
                     # Update the states
                     S = Dolo.transition(model, m, s[n, :], x[i][n, :], M, p)

                     # Update value function
                      E_V[i][n, :] += w*drv(i, j, S)
                 end
            end
        end


        err = 0.0
        for i in 1:length(v)
            broadcast!((_u, _E_v) -> _u + β * _E_v, v[i], u[i], E_V[i])
            err = max(err, maxabs(v[i] - v0[i]))
            copy!(v0[i], v[i])
            fill!(E_V[i], 0.0)
        end

        verbose && @printf "%-6i%-12.2e\n" it err

        set_values!(drv, v0)
    end

    return drv
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
    for j=1:n_inodes(dprocess, i)
        M = inode(dprocess, i, j)
        w = iweight(dprocess, i, j)
        S = Dolo.transition(model, m, s, x0, M, p)
        E_V += w*drv(i, j, S)[1]
    end
    u = Dolo.felicity(model, m, s, x0, p)[1]
    E_V = u + β.*E_V
    return E_V
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
function solve_policy(model, pdr; maxit::Int=1000, verbose::Bool=true, details::Bool=false)

    # get grid for endogenous
    grid = model.grid

    β=model.calibration.flat[:beta]

    # process = dr.process
    process = model.exogenous
    dprocess = discretize(process)

    dr = CachedDecisionRule(pdr, dprocess)
    # compute the value function
    absmax(x) = max( [maximum(abs(x[i])) for i=1:length(x)]... )
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
    verbose && println("Evaluating initial policy")

    drv = evaluate_policy(model, dr; maxit=1000, verbose=verbose)

    verbose && println("Evaluating initial policy (done)")

    v0 = [drv(i, endo_nodes) for i=1:nsd]
    v = deepcopy(v0)

    #Preparation for a loop
    tol_x = 1e-8
    tol_v = 1e-8
    tol_eval = 1e-8
    maxit_eval = 1000
    err_v = 10.0
    err_x = 10.0
    err_eval = 10.0

    it = 0
    it_eval = 0


    optim_opts = Optim.OptimizationOptions(x_tol=1e-9, f_tol=1e-9)

    mode = :improve

    converged = false

    while !converged

        # it += 1
        if (mode == :eval)
    #
            it_eval = 0
            converged_eval = false
            while !converged_eval
                it_eval += 1
                for i = 1:size(res, 1)
                    m = node(dprocess, i)
                    for n = 1:N
                        s = endo_nodes[n, :]
    #                     # update vals
                        nv = update_value(model, β, dprocess, drv, i, s, x0[i][n, :], p)
                        v[i][n, 1] = nv
                    end
                end
                # compute diff in values
                err_eval = 0.0
                for i in 1:nsd
                    err_eval = max(err_eval, maxabs(v[i] - v0[i]))
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
    #
            it += 1
            for i = 1:size(res, 1)
                m = node(dprocess, i)
                for n = 1:N
                    s = endo_nodes[n, :]
                    # optimize vals
                    fobj(u) = -update_value(model, β, dprocess, drv, i, s, u, p)*1000.0
                    lower = x_lb[i][n, :]
                    upper = x_ub[i][n, :]
                    upper = clamp!(upper, -Inf, 1000000)
                    lower = clamp!(lower, -1000000, Inf)
                    initial_x = x0[i][n, :]
                    # try
                    results = optimize(
                        DifferentiableFunction(fobj), initial_x, lower, upper,
                        Fminbox(), optimizer=NelderMead,
                        x_tol=1e-10, f_tol=1e-10
                    )
                    xn = Optim.minimizer(results)
                    nv = -Optim.minimum(results)/1000.0
                    x[i][n, :] = xn
                    v[i][n, 1] = nv
                end
            end

            # compute diff in values
            err_v = 0.0
            for i in 1:nsd
                err_v = max(err_v, maxabs(v[i] - v0[i]))
                copy!(v0[i], v[i])
            end
            # compute diff in policy
            err_x = 0.0
            for i in 1:nsd
                err_x = 0
                err_x = max(err_x, maxabs(x[i] - x0[i]))
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
    #
        end
    #
    end


    if !details
        dr = CachedDecisionRule(dprocess, grid, x0)
        return (dr.dr, drv.dr)
    else
        converged_x = err_x<tol_x
        converged_v = err_v<tol_v
        ValueIterationResult(dr.dr, drv.dr, it, true, converged_x, tol_x, err_x, converged_v, tol_v, err_v)
    end
end
