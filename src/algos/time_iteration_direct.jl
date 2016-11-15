import Dolo

function time_iteration_direct(model, process; verbose=true, maxit=100)

  # Grid
  gg = model.options.grid
  grid = CartesianGrid(gg.a, gg.b, gg.orders) # temporary compatibility
  endo_nodes = nodes(grid)
  N = size(endo_nodes,1)

  # Discretized exogenous process
  dprocess = discretize(process)
  number_of_smooth_drs(dprocess) = max(n_nodes(dprocess),1)
  nsd = number_of_smooth_drs(dprocess)

  p = model.calibration[:parameters] :: Vector{Float64}

  init_dr = ConstantDecisionRule(model.calibration[:controls])
  # initial guess for controls
  x0 = [evaluate(init_dr, i, endo_nodes) for i=1:nsd]

  function stack(x::Array{Array{Float64,2},1})
       return cat(1,x...)
  end

  absmax(x) = max([maximum(abs(x[i])) for i=1:length(x)]...)

  # create decision rule (which interpolates x0)
  dr = DecisionRule(process, grid, x0)
  # Define controls of tomorrow
  x1 =[zeros(N,2) for i=1:number_of_smooth_drs(dprocess)]
  # define states of today
  s=deepcopy(endo_nodes);

  # loop option
  tol = 1e-8
  it = 0
  err = 1

  ###############################   Iteration loop

  while it<maxit && err>tol

    it+=1
    # dr = DecisionRule(process, grid, x0)
    set_values(dr, x0)
    xx0 = stack(x0)
    # Compute expectations function E_f and states of tomorrow
    E_f = [zeros(N,1) for i=1:number_of_smooth_drs(dprocess)]
    S = zeros(size(s))

    for i=1:size(E_f,1)
        m = node(dprocess,i)  ::Vector{Float64}
        for j=1:n_inodes(dprocess,i)
            M = inodes(dprocess,i,j) ::Vector{Float64}
            w = iweights(dprocess,i,j) ::Float64
            # Update the states
            for n=1:N
                S[n,:] = Dolo.transition(model, m, s[n,:], x0[i][n,:], M, p)
            end
            # interpolate controles conditional states of tomorrow
            X = evaluate(dr, i, j, S)
            # Compute expectations as a weited average of the exo states w_j
            for n=1:N
                E_f[i][n,:] += w*Dolo.expectation(model, M, S[n,:], X[n,:], p)
            end
        end
        # compute controles of tomorrow
        for n=1:N
           x1[i][n,:] = Dolo.direct_response(model, m, s[n,:], E_f[i][n,:], p)
        end
    end

    xx1 = stack(x1)

    # update error
    err = maximum(abs(xx1 - xx0))

    # Update control vector
    x0 = x1

    if verbose
        println("It: ", it, " ; SA: ", err, " ; nit: ", it)
    end
  end
  return dr
end



function time_iteration_direct(model; verbose=true)
    process = model.exogenous
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration_direct(model, process; verbose=verbose, maxit = 100)
end
