using Optim


"""
Evaluate the value function under the given decision rule.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dr`: Current guess for the decision rule.
# Returns
* `drv`: Value function.
"""
function evaluate_policy(model, dr; verbose=true, maxit=100, )

    β = model.calibration.flat[:beta]

    # get grid for endogenous
    gg = model.options.grid
    grid = CartesianGrid(gg.a, gg.b, gg.orders) # temporary compatibility

    process = model.exogenous
    dprocess = discretize(process)

    # compute the value function
    absmax(x) = max([maximum(abs(x[i])) for i=1:length(x)]...)
    p = model.calibration[:parameters] :: Vector{Float64}

    endo_nodes = nodes(grid)
    # Number of endogenous nodes
    N = size(endo_nodes,1)
    number_of_smooth_drs(dprocess) = max(n_nodes(dprocess),1)
    res = [zeros(N,1) for i=1:number_of_smooth_drs(dprocess)]

    # Value function : v_t = u + β*V_{t+1}
    # u is constant because the initial guess is a constant decision rule:
    # whatever the state you make the same decision.
    u = deepcopy(res)
    v = deepcopy(res)
    v0 = deepcopy(res)
    # Expected utility
    E_V = deepcopy(res)
    m = nothing
    M = nothing
    # States at time t+1
    S =  copy(endo_nodes)

    # Controls at time t
    x= [evaluate(dr, i, endo_nodes) for i=1:number_of_smooth_drs(dprocess)]
    s = deepcopy(endo_nodes)
    for i=1:size(res,1)
        m = node(dprocess,i)  ::Vector{Float64}
        for n=1:N
            u[i][n,:] = Dolo.felicity(model, m, s[n,:], x[i][n,:], p)[1]
        end
    end

    # Initial guess for the value function
    v0 = deepcopy(u)

    #Preparation for a loop
    tol = 1e-6
    err=10
    Err=zeros(maxit)
    it = 0

    drv = DecisionRule(process, grid, v0)
    # return drv
    while (err>tol && it<maxit)
        it +=1
        # Interpolate v0 on the grid
        # Compute value function
        for i=1:size(res,1)
            m = node(dprocess,i)  ::Vector{Float64}
            for j=1:n_inodes(dprocess,i)
                M = inodes(dprocess,i,j) ::Vector{Float64}
                w = iweights(dprocess,i,j) ::Float64
                 # Update the states
                 for n=1:N
                     S[n,:] = Dolo.transition(model, m, s[n,:], x[i][n,:], M, p)
                 end
                 # Update value function
                 for n=1:N
                     E_V[i][n,1] += w*evaluate(drv,i, j, S[n,:])[1]
                 end
            end
        end

        v = u + β.*E_V
        err = absmax(v-v0)
        v0 = deepcopy(v)
        E_V *= 0

        if verbose
            println("It: ", it, " ; SA: ", err, " ; nit: ", it)
        end

        set_values(drv,v0)
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
function update_value(model, β::Float64, dprocess, drv, i, s::Vector{Float64}, x0::Vector{Float64}, p::Vector{Float64})

    m = node(dprocess,i)  ::Vector{Float64}
    E_V = 0.0
    for j=1:n_inodes(dprocess,i)
        M = inodes(dprocess,i,j) ::Vector{Float64}
        w = iweights(dprocess,i,j) ::Float64
        S = Dolo.transition(model, m, s, x0, M, p)
        E_V += w*evaluate(drv, i, j, S)[1]
    end
    u = Dolo.felicity(model, m, s, x0, p)[1]
    E_V = u + β.*E_V
    return E_V
end


"""
Solve for the value function and associated decision rule using value function iteration.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dr`: Initial guess for the decision rule.
# Returns
* `dr`: Solved decision rule object.
* `drv`: Solved value function object.
"""
function solve_policy(model, dr; verbose=true, maxit=5000, )

    β = model.calibration.flat[:beta]
    # get grid for endogenous
    gg = model.options.grid
    grid = CartesianGrid(gg.a, gg.b, gg.orders) # temporary compatibility

    # process = dr.process
    process = model.exogenous
    dprocess = discretize(process)

    # compute the value function
    absmax(x) = max([maximum(abs(x[i])) for i=1:length(x)]...)
    p = model.calibration[:parameters] :: Vector{Float64}

    endo_nodes = nodes(grid)
    # Number of endogenous nodes
    N = size(endo_nodes,1)
    number_of_smooth_drs(dprocess) = max(n_nodes(dprocess),1)
    nsd = number_of_smooth_drs(dprocess)
    res = [zeros(N,1) for i=1:number_of_smooth_drs(dprocess)]

    x_lb = Array{Float64,2}[cat(1,[Dolo.controls_lb(model,node(dprocess,i) ,endo_nodes[n,:],p)' for n=1:N]...) for i=1:nsd]
    x_ub = Array{Float64,2}[cat(1,[Dolo.controls_ub(model,node(dprocess,i),endo_nodes[n,:],p)' for n=1:N]...) for i=1:nsd]


    # Value function : v_t = u + β*V_{t+1}
    # u is constant because the initial guess is a constant decision rule:
    # whatever the state you make the same decision.
    u = deepcopy(res)

    # Expected utility
    m = nothing
    M = nothing
    # States at time t+1
    S =  copy(endo_nodes)

    # Controls at time t
    x= [evaluate(dr, i, endo_nodes) for i=1:number_of_smooth_drs(dprocess)]
    x0 = deepcopy(x)

    if verbose
      println("Evaluating initial policy")
    end
    drv = evaluate_policy(model, dr, maxit=500, verbose=false)
    if verbose
      println("Evaluating initial policy (done)")
    end

    v0 = [evaluate(drv, i, endo_nodes) for i=1:number_of_smooth_drs(dprocess)]
    v = deepcopy(v0)
    #Preparation for a loop
    tol = 1e-6
    err=10
    err_x = 10
    Err=zeros(maxit)
    it = 0

    n_eval = 50

    while ( (err>tol) || (err_x>tol) ) && (it<maxit)
        it +=1
        if (n_eval*div(it,n_eval) == it) || (err<tol)
            optim = true
        else
            optim = false
        end
        # Interpolate v0 on the grid
        # Compute value function
        v = deepcopy(v0)
        x = deepcopy(x0)
        for i=1:size(res,1)
            for n=1:N
                m = node(dprocess,i)  ::Vector{Float64}
                s = endo_nodes[n,:]
                # update vals
                if !optim
                    nv = update_value(model, β, dprocess, drv, i, s, x0[i][n,:], p)
                    v[i][n,1] = nv
                else
                # optimize vals
                    fobj(u) = -update_value(model, β, dprocess, drv, i, s, u, p)*1000.0
                    lower = x_lb[i][n,:]
                    upper = x_ub[i][n,:]
                    upper[upper.>1000000] = 1000000.0
                    lower[lower.<-1000000] = -1000000.0
                    initial_x = x0[i][n,:]
                    # try
                    results = optimize(DifferentiableFunction(fobj), initial_x, lower, upper, Fminbox(), optimizer = NelderMead)
                    # results = optimize(DifferentiableFunction(fobj), initial_x, NelderMead())
                    # println(results)
                    xn = Optim.minimizer(results)
                    nv = -Optim.minimum(results)/1000.0
                    x[i][n,:] = xn
                    v[i][n,1] = nv

                end
                # println(x1)
                # x[i,1][n,:] = x1
            end
        end


        err = absmax(v-v0)
        if optim
          err_x = absmax(x-x0)
        end

        v0 = deepcopy(v)
        if optim
            x0 = deepcopy(x)

        end
        if verbose
            if optim
                println("It: ", it, " ; SA: ", err, " ; SA_x: ", err_x, " ; nit: ", it)
            else
                println("It: ", it, " ; SA: ", err, " ; nit: ", it)
            end
        end

        set_values(drv,v0)

    end

    dr = DecisionRule(process, grid, x0)
    return (dr, drv)
end
