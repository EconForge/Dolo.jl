
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
function residual(model, dprocess::AbstractDiscretizedProcess, s, x::Array{Array{Float64,2},1}, p, dr)
    N = size(s, 1)
    res = [zeros(size(x[1])) for i=1:length(x)]
    S = zeros(size(s))
    # X = zeros(size(x[1]))
    for i=1:size(res, 1)
        m = node(dprocess, i)
        for j=1:n_inodes(dprocess, i)
            M = inode(dprocess, i, j)
            w = iweight(dprocess, i, j)
            # Update the states
            for n=1:N
                S[n, :] = Dolo.transition(model, m, s[n, :], x[i][n, :], M, p)
            end

            X = dr(i, j, S)
            for n=1:N
                res[i][n, :] += w*Dolo.arbitrage(model, m, s[n, :], x[i][n, :], M, S[n, :], X[n, :], p)
            end
        end
    end
    return res
end

using NLsolve


function residual(model, dprocess, s, x::Array{Float64,2}, p, dr)
    n_m = max(1, n_nodes(dprocess))
    xx = destack0(x, n_m)
    res = residual(model, dprocess, s, xx, p, dr)
    return stack0(res)
end

"""
#TODO
"""
function destack0(x::Array{Float64,2}, n_m::Int)
    N = div(size(x, 1), n_m)
    xx = reshape(x, N, n_m, size(x, 2))
    return Array{Float64,2}[xx[:, i, :] for i=1:n_m]
end

"""
#TODO
"""
function stack0(x::Array{Array{Float64,2},1})::Array{Float64,2}
     return cat(1, x...)
end

type TimeIterationResult
    dr::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    x_converged::Bool
    x_tol::Float64
    err::Float64
end

converged(r::TimeIterationResult) = r.x_converged
function Base.show(io::IO, r::TimeIterationResult)
    @printf io "Results of Time Iteration Algorithm\n"
    @printf io " * Complementarities: %s\n" string(r.complementarities)
    @printf io " * Decision Rule type: %s\n" string(typeof(r))
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
function time_iteration(model, dprocess::AbstractDiscretizedProcess, init_dr;
      verbose::Bool=true, maxit::Int=100, tol::Float64=1e-8, details::Bool=false)

    # get grid for endogenous
    grid = model.grid
    # grid = CartesianGrid(gg.a, gg.b, gg.orders)  # temporary compatibility

    endo_nodes = nodes(grid)
    N = size(endo_nodes, 1)
    n_s_endo = size(endo_nodes,2)

    n_s_exo = n_nodes(dprocess)

    # initial guess
    # number of smooth decision rules
    nsd = max(n_nodes(dprocess), 1)

    p = model.calibration[:parameters]

    x0 = [init_dr(i, endo_nodes) for i=1:nsd]

    n_x = length(model.calibration[:controls])
    lb = Array(Float64, N*nsd, n_x)
    ub = Array(Float64, N*nsd, n_x)
    ix = 0
    for i in 1:nsd
        node_i = node(dprocess, i)
        for n in 1:N
            ix += 1
            endo_n = endo_nodes[n, :]
            lb[ix, :] = Dolo.controls_lb(model, node_i, endo_n, p)
            ub[ix, :] = Dolo.controls_ub(model, node_i, endo_n, p)
        end
    end

    # create decision rule (which interpolates x0)
    dr = CachedDecisionRule(dprocess, grid, x0)

    # loop option
    init_res = residual(model, dprocess, endo_nodes, x0, p, dr)
    it = 0
    err = maximum(abs, stack0(init_res))
    err_0 = err

    verbose && @printf "%-6s%-12s%-12s%-5s\n" "It" "SA" "gain" "nit"
    verbose && println(repeat("-", 35))
    verbose && @printf "%-6i%-12.2e%-12.2e%-5i\n" 0 err NaN 0

    while it<maxit && err>tol

        it += 1

        set_values!(dr, x0)

        xx0 = stack0(x0)
        fobj(u) = residual(model, dprocess, endo_nodes, u, p, dr)
        xx1, nit = serial_solver(fobj, xx0, lb, ub, maxit=10, verbose=false)
        x1 = destack0(xx1, nsd)

        err = maximum(abs, xx1 - xx0)
        copy!(x0, x1)
        gain = err / err_0
        err_0 = err

        verbose && @printf "%-6i%-12.2e%-12.2e%-5i\n" it err gain nit
    end

    # TODO: somehow after defining `fobj` the `dr` object gets `Core.Box`ed
    #       making the return type right here non-inferrable.
    if !details
        return dr.dr
    else
        converged = err<tol
        TimeIterationResult(dr.dr, it, true, converged, tol, err)
    end

end


# get stupid initial rule
function time_iteration(model, dprocess::AbstractDiscretizedProcess; kwargs...)
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration(model, dprocess, init_dr;  kwargs...)
end


function time_iteration(model, init_dr; kwargs...)
    dprocess = discretize( model.exogenous )
    return time_iteration(model, dprocess, init_dr; kwargs...)
end


function time_iteration(model; kwargs...)
    dprocess = discretize( model.exogenous )
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration(model, dprocess, init_dr; kwargs...)
end
