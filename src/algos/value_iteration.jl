function evaluate_policy(model, dr; verbose=true, maxit=100, )

    β = model.calibration.flat[:beta]

    # get grid for endogenous
    gg = model.options.grid
    grid = CartesianGrid(gg.a, gg.b, gg.orders) # temporary compatibility

    process = dr.process
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
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    x= [evaluate(init_dr, i, endo_nodes) for i=1:number_of_smooth_drs(dprocess)]
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
    maxit = 700
    err=10
    Err=zeros(maxit)
    it = 0

    drv = DecisionRule(process, grid, v0)

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
                     E_V[i][n,:] += w*evaluate(drv,i, j, S[n,:])[1]
                 end
            end
        end
        # Instead of
        # V_t <= U(s_t, x_t) + \\beta sum w_M  E[ V(g(s_t, x_t, M_t) ]
        # I compute w_M  E[ V(g(s_t, x_t, M_t)
        # and then since v0 = v, then  v = v0 + β.*V = u + β*(w*Vi) + β*(w*Vi) + β*(w*Vi) + ...
        s = deepcopy(S)

        v = u + β.*E_V
        err = absmax(v-v0)
        v0 = deepcopy(v)
        E_V *= 0

        if verbose
            println("It: ", it, " ; SA: ", err, " ; nit: ", it)
        end

        set_values(drv,v0)
    end
    return (drv)
end
