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
    endo_nodes = nodes(grid)
    N = size(endo_nodes, 1)

    # Discretized exogenous process
    number_of_smooth_drs(dprocess) = max(n_nodes(dprocess), 1)
    nsd = number_of_smooth_drs(dprocess)

    p = model.calibration[:parameters]

    # initial guess for controls
    x0 = [init_dr(i, endo_nodes) for i=1:nsd]

    ti_trace = trace ? IterationTrace([x0]) : nothing


    # set the bound for the controls to check during the iterations not to violate them
    x_lb = Array{Float64,2}[cat(1, [Dolo.controls_lb(model, node(dprocess, i), endo_nodes[n, :], p)' for n=1:N]...) for i=1:nsd]
    x_ub = Array{Float64,2}[cat(1, [Dolo.controls_ub(model, node(dprocess, i), endo_nodes[n, :], p)' for n=1:N]...) for i=1:nsd]

    # create decision rule (which interpolates x0)
    dr = CachedDecisionRule(dprocess, grid, x0)

    # Define controls of tomorrow
    x1 = [zeros(N, 2) for i=1:number_of_smooth_drs(dprocess)]

    # define states of today
    s = deepcopy(endo_nodes);

    # loop option
    it = 0
    err_0 = NaN
    err = 1.0

    log = TimeIterationLog()
    initialize(log, verbose=verbose)
    append!(log; verbose=verbose, it=0, err=NaN, gain=NaN, time=0.0, nit=NaN)

    maxabsdiff(_a, _b) = maximum(abs, _a - _b)

    ###############################   Iteration loop

    while it<maxit && err>tol_η

        it += 1

        tic()


        set_values!(dr, x0)
        # Compute expectations function E_f and states of tomorrow
        E_f = [zeros(N, 1) for i=1:number_of_smooth_drs(dprocess)]
        S = zeros(size(s))

        for i in 1:size(E_f, 1)
            m = node(dprocess, i)
            for (w, M, j) in get_integration_nodes(dprocess,i)
                # Update the states
                S[:,:] = Dolo.transition(model, m, s, x0[i], M, p)
                # interpolate controles conditional states of tomorrow
                X = dr(i, j, S)
                # Compute expectations as a weighted average of the exo states w_j
                E_f[i][:,:] += w*Dolo.expectation(model, M, S, X, p)

            end
            # compute controles of tomorrow
            x1[i][:] = Dolo.direct_response(model, m, s, E_f[i], p)
        end

        err = 0.0
        for i in 1:size(x1, 1)
            # apply bounds
            broadcast!(clamp, x1[i], x1[i], x_lb[i], x_ub[i])
            # update error
            err = max(err, maximum(abs, x1[i] - x0[i]))
            # copy controls back into x0
            copy!(x0[i], x1[i])
        end

        gain = err/err_0
        err_0 = err

        elapsed = toq()

        append!(log; verbose=verbose, it=it, err=err, gain=gain, time=elapsed, nit=NaN)
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
