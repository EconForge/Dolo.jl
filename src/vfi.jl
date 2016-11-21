module temp
    using Dolo
    using Optim
    import Optim
    import Dolo:CartesianGrid, n_nodes, discretize, nodes, ConstantDecisionRule, node, DecisionRule, n_inodes, inodes, iweights
    import Dolo:set_values
    function update_value(model, β, dprocess, drv, i, s, x0, p)
        m = node(dprocess,i)  ::Vector{Float64}
        E_V = 0.0
        for j=1:n_inodes(dprocess,i)
            M = inodes(dprocess,i,j) ::Vector{Float64}
            w = iweights(dprocess,i,j) ::Float64
            S = Dolo.transition(model, m, s, x0, M, p)
            E_V += w*evaluate(drv, i, j, S)[1]
        end
        u = Dolo.felicity(model, m, s, x0, p)[1]
        E_V = u + β*E_V
        return E_V
    end

    function solve_policy(model, dr; verbose=true, maxit=10000, )

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
        nsd = number_of_smooth_drs(dprocess)
        res = [zeros(N,1) for i=1:number_of_smooth_drs(dprocess)]

        x_lb = Array{Float64,2}[cat(1,[Dolo.controls_lb(model,node(dprocess,i) ,endo_nodes[n,:],p)' for n=1:N]...) for i=1:nsd]
        x_ub = Array{Float64,2}[cat(1,[Dolo.controls_ub(model,node(dprocess,i),endo_nodes[n,:],p)' for n=1:N]...) for i=1:nsd]


        # Value function : v_t = u + β*V_{t+1}
        # u is constant because the initial guess is a constant decision rule:
        # whatever the state you make the same decision.
        u = deepcopy(res)
        v = deepcopy(res)
        v0 = deepcopy(res)
        # Expected utility
        m = nothing
        M = nothing
        # States at time t+1
        S =  copy(endo_nodes)

        # Controls at time t
        init_dr = ConstantDecisionRule(model.calibration[:controls])
        x= [evaluate(init_dr, i, endo_nodes) for i=1:number_of_smooth_drs(dprocess)]
        x0 = deepcopy(x)

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

        n_eval = 500
        # while (err>tol && it<maxit)
        while (it<maxit)
            it +=1
            if n_eval*div(it,n_eval) == it
                optim = true
            else
                optim = false
            end
            optim = false
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
                        upper = [2, 50.0]
                        initial_x = x0[i][n,:]
                        try
                            results = optimize(DifferentiableFunction(fobj), initial_x, lower, upper, Fminbox(), optimizer = GradientDescent)
                            xn = Optim.minimizer(results)
                            nv = -Optim.minimum(results)/1000.0
                            x[i][n,:] = xn
                            v[i][n,1] = nv
                        catch
                            return  initial_x, lower, upper, drv
                        end
                    end
                    # println(x1)
                    # x[i,1][n,:] = x1
                end
            end

            s = deepcopy(S)

            err = absmax(v-v0)
            err_x = absmax(x-x0)

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
        return (drv, x0)
    end

end