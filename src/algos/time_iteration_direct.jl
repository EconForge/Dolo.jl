"""
Computes a global solution for a model via backward time iteration.
The time iteration is  applied directly to the decision rule of the model.

If the initial guess for the decision rule is not explicitly provided, the initial guess is provided by `ConstantDecisionRule`.
If the stochastic process for the model is not explicitly provided, the process is taken from the default provided by the model object, `model.exogenous`.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `process`: The stochastic process associated with the exogenous variables in the model.
* `init_dr`: Initial guess for the decision rule.
# Returns
* `dr`: Solved decision rule.
"""
function time_iteration_direct(model, dprocess::AbstractDiscretizedProcess,
                               grid, init_dr::AbstractDecisionRule;
                               verbose::Bool=true, maxit::Int=100, trace::Bool=false,
                               tol_η::Float64=1e-7)

    # Grid
    endo_nodes = nodes(ListOfPoints,grid)
    N = length(endo_nodes)

    # Discretized exogenous process
    number_of_smooth_drs(dprocess) = max(n_nodes(dprocess), 1)
    n_m = nsd = number_of_smooth_drs(dprocess)

    p = SVector(model.calibration[:parameters]...)

    # initial guess for controls
    x0 = [init_dr(i, endo_nodes) for i=1:nsd]

    ti_trace = trace ? IterationTrace([x0]) : nothing

    n_x = length(model.calibration[:controls])
    complementarities = true
    if complementarities == true
        x_lb = [controls_lb(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m]
        x_ub = [controls_ub(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m]
        BIG = 100000
        for i=1:n_m
          for n=1:N
            x_lb[i][n] = max.(x_lb[i][n],-BIG)
            x_ub[i][n] = min.(x_ub[i][n], BIG)
          end
        end
    end

    # create decision rule (which interpolates x0)
    dr = CachedDecisionRule(dprocess, grid, x0)

    # Define controls of tomorrow
    x1 = deepcopy(x0)

    # define states of today
    s = deepcopy(endo_nodes);

    # loop option
    it = 0
    err_0 = NaN
    err = 1.0

    log = TimeIterationLog()
    initialize(log, verbose=verbose)
    append!(log; verbose=verbose, it=0, err=NaN, gain=NaN, epsilon=NaN, time=0.0, nit=NaN)


    ###############################   Iteration loop
    n_h = length(model.symbols[:expectations])


    while it<maxit && err>tol_η

        it += 1

        t1 = time_ns()

        E_f = [zeros(Point{n_h},N) for i=1:number_of_smooth_drs(dprocess)]

        set_values!(dr, x0)
        # Compute expectations function E_f and states of tomorrow

        # S = zeros(size(s))

        for i in 1:size(E_f, 1)
            m = node(Point,dprocess, i)
            for (w, M, j) in get_integration_nodes(Point,dprocess,i)
                # Update the states
                # S[:,:] = Dolo.transition(model, m, s, x0[i], M, p)
                S = Dolo.transition(model, m, s, x0[i], M, p)
                # interpolate controles conditional states of tomorrow
                X = dr(i, j, S)
                # Compute expectations as a weighted average of the exo states w_j
                E_f[i] += w*Dolo.expectation(model, M, S, X, p)

            end
            # compute controles of tomorrow
            x1[i][:] = Dolo.direct_response(model, m, s, E_f[i], p)
        end

        err = 0.0
        for i in 1:size(x1, 1)
            # apply bounds # TODO: improve
            x1[i] = clamp.(x1[i], x_lb[i], x_ub[i])
            # update error
            err = max(err, maxabs(x1[i] - x0[i]))
            # copy controls back into x0
            copy!(x0[i], x1[i])
        end

        gain = err/err_0
        err_0 = err

        elapsed = time_ns()-t1

        append!(log; verbose=verbose, it=it, err=err, gain=gain, time=elapsed, epsilon=NaN, nit=NaN)

    end

    finalize(log, verbose=verbose)


    converged = err<tol_η
    TimeIterationResult(dr.dr, it, true, dprocess, converged, tol_η, err, log, ti_trace)

end

# get grid for endogenous
function time_iteration_direct(model, dprocess, init_dr::AbstractDecisionRule; grid=Dict(), kwargs...)
    grid = get_grid(model, options=grid)
    return time_iteration_direct(model, dprocess, grid, init_dr;  kwargs...)
end

# get stupid initial rule
function time_iteration_direct(model, dprocess::AbstractDiscretizedProcess; grid=Dict(), kwargs...)

    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration_direct(model, dprocess, init_dr;  grid=grid, kwargs...)
end


function time_iteration_direct(model, init_dr::AbstractDecisionRule; grid=Dict(), kwargs...)
    dprocess = discretize( model.exogenous )
    return time_iteration_direct(model, dprocess, init_dr; grid=grid, kwargs...)
end


function time_iteration_direct(model; grid=Dict(), kwargs...)
    dprocess = discretize( model.exogenous )
    init_dr = ConstantDecisionRule(model.calibration[:controls])
    return time_iteration_direct(model, dprocess, init_dr; grid=grid, kwargs...)
end
