function time_iteration_direct(model, process, init_dr, verbose=true, maxit=100, tol=1e-8)

    # Grid
    gg = model.options.grid
    grid = CartesianGrid(gg.a, gg.b, gg.orders) # temporary compatibility
    endo_nodes = nodes(grid)
    N = size(endo_nodes, 1)

    # Discretized exogenous process
    dprocess = discretize(process)
    number_of_smooth_drs(dprocess) = max(n_nodes(dprocess), 1)
    nsd = number_of_smooth_drs(dprocess)

    p = model.calibration[:parameters]

    # initial guess for controls
    x0 = [init_dr(i, endo_nodes) for i=1:nsd]

    # set the bound for the controls to check during the iterations not to violate them
    x_lb = Array{Float64,2}[cat(1, [Dolo.controls_lb(model, node(dprocess, i), endo_nodes[n, :], p)' for n=1:N]...) for i=1:nsd]
    x_ub = Array{Float64,2}[cat(1, [Dolo.controls_ub(model, node(dprocess, i), endo_nodes[n, :], p)' for n=1:N]...) for i=1:nsd]

    # create decision rule (which interpolates x0)
    dr = DecisionRule(process, grid, x0)

    # Define controls of tomorrow
    x1 = [zeros(N, 2) for i=1:number_of_smooth_drs(dprocess)]

    # define states of today
    s = deepcopy(endo_nodes);

    foobar = 10

    # loop option
    it = 0
    err = 1.0

    verbose && @printf "%-6s%-12s\n" "It" "SA"
    verbose && println(repeat("-", 14))

    maxabsdiff(_a, _b) = maxabs(_a - _b)

    ###############################   Iteration loop

    while it<maxit && err>tol

      it+=1
      # dr = DecisionRule(process, grid, x0)
      set_values!(dr, x0)
      # Compute expectations function E_f and states of tomorrow
      E_f = [zeros(N, 1) for i=1:number_of_smooth_drs(dprocess)]
      S = zeros(size(s))

      for i=1:size(E_f, 1)
          m = node(dprocess, i)
          for j=1:n_inodes(dprocess, i)
              M = inodes(dprocess, i, j)
              w = iweights(dprocess, i, j)
              # Update the states
              for n=1:N
                  S[n, :] = Dolo.transition(model, m, s[n, :], x0[i][n, :], M, p)
              end
              # interpolate controles conditional states of tomorrow
              X = dr(i, j, S)
              # Compute expectations as a weited average of the exo states w_j
              for n=1:N
                  E_f[i][n, :] += w*Dolo.expectation(model, M, S[n, :], X[n, :], p)
              end
          end
          # compute controles of tomorrow
          for n=1:N
             x1[i][n, :] = Dolo.direct_response(model, m, s[n, :], E_f[i][n, :], p)
          end
      end

      err = 0.0
      for i in 1:size(x1, 1)
          # apply bounds
          broadcast!(clamp, x1[i], x1[i], x_lb[i], x_ub[i])

          # update error
          err = max(err, maxabs(x1[i] - x0[i]))

          # copy controls back into x0
          copy!(x0[i], x1[i])
      end

      verbose && @printf "%-6i%-12.2e\n" it err
    end

    return dr
end


function time_iteration_direct(model, process::AbstractExogenous; kwargs...)
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration_direct(model, process, init_dr; kwargs...)
end


function time_iteration_direct(model, init_dr::AbstractDecisionRule; kwargs...)
    process = model.exogenous
    return time_iteration_direct(model, process, init_dr; kwargs...)
end


function time_iteration_direct(model; kwargs...)
    process = model.exogenous
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration_direct(model, process, init_dr; kwargs...)
end
