
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
function residual(model, dprocess, s, x::Array{Array{Float64,2},1}, p, dr)
    N = size(s,1)
    res = [zeros(size(x[1])) for i=1:length(x)]
    S = zeros(size(s))
    # X = zeros(size(x[1]))
    for i=1:size(res,1)
        m = node(dprocess,i)  ::Vector{Float64}
        for j=1:n_inodes(dprocess,i)
            M = inodes(dprocess,i,j) ::Vector{Float64}
            w = iweights(dprocess,i,j) ::Float64
            # Update the states
            for n=1:N
                S[n,:] = Dolo.transition(model, m, s[n,:], x[i][n,:], M, p)
            end
            X = evaluate(dr, i, j, S)
            for n=1:N
                res[i][n,:] += w*Dolo.arbitrage(model, m, s[n,:], x[i][n,:], M, S[n,:], X[n,:], p)
            end
        end
    end
    return res
end

using NLsolve


function residual(model, dprocess, s, x::Array{Float64,2}, p, dr)
    n_m = max(1,n_nodes(dprocess))
    xx = destack(x,n_m)
    res = residual(model, dprocess, s, xx, p, dr)
    return stack(res)
end


"""
#TODO
"""
function destack(x::Array{Float64,2},n_m::Int)
    N = div(size(x,1),n_m)
    xx = reshape(x,N,n_m,size(x,2))
    return Array{Float64,2}[xx[:,i,:] for i=1:n_m]
end


"""
#TODO
"""
function stack(x::Array{Array{Float64,2},1})
     return cat(1,x...)
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
function time_iteration(model, process, init_dr; verbose=true, maxit=100, tol=1e-8)

    # get grid for endogenous
    gg = model.options.grid
    grid = CartesianGrid(gg.a, gg.b, gg.orders) # temporary compatibility

    endo_nodes = nodes(grid)
    N = size(endo_nodes,1)
    n_s_endo = size(endo_nodes,2)

    dprocess = discretize(process)
    n_s_exo = n_nodes(dprocess)

    # initial guess
    number_of_smooth_drs(dprocess) = max(n_nodes(dprocess),1)
    nsd = number_of_smooth_drs(dprocess)

    p = model.calibration[:parameters] :: Vector{Float64}

    x0 = [evaluate(init_dr, i, endo_nodes) for i=1:nsd]


    x_lb = Array{Float64,2}[cat(1,[Dolo.controls_lb(model,node(dprocess,i) ,endo_nodes[n,:],p)' for n=1:N]...) for i=1:nsd]
    x_ub = Array{Float64,2}[cat(1,[Dolo.controls_ub(model,node(dprocess,i),endo_nodes[n,:],p)' for n=1:N]...) for i=1:nsd]

    lb = stack(x_lb)
    ub = stack(x_ub)

    n_x = length(model.calibration[:controls])

    absmax(x) = max([maximum(abs(x[i])) for i=1:length(x)]...)

    # create decision rule (which interpolates x0)
    dr = DecisionRule(process, grid, x0)


    # loop option

    it = 0
    err = 1.0
    err_0 = err
    maxit_inner = 10


    init_res = residual(model, dprocess, endo_nodes, x0, p, dr)
    println("Initial Residual: ", maximum(abs(stack(init_res))))

    while it<maxit && err>tol

        it+=1

        # dr = DecisionRule(process, grid, x0)
        set_values(dr, x0)

        xx0 = stack(x0)
        fobj(u) = residual(model, dprocess, endo_nodes, u, p, dr)
        xx1, nit = serial_solver(fobj, xx0, lb, ub, maxit=maxit_inner, verbose=false)
        x1 = destack(xx1, nsd)

        err = maximum(abs(xx1 - xx0))
        x0 = x1
        gain = err/err_0
        if verbose
            println("It: ", it, " ; SA: ", err, " ; gain: ", gain, " ; nit: ", nit)
        end
        err_0 = err
    end

    return dr

end


# get stupid initial rule
function time_iteration(model, process::AbstractExogenous; kwargs...)
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration(model, process, init_dr;  kwargs...)
end


function time_iteration(model, init_dr::AbstractDecisionRule; kwargs...)
    process = model.exogenous
    return time_iteration(model, process, init_dr; kwargs...)
end


function time_iteration(model; kwargs...)
    process = model.exogenous
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration(model, process, init_dr; kwargs...)
end
