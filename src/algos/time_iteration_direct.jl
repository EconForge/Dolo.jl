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

    d = length(endo_nodes[1])
    # Discretized exogenous process
    number_of_smooth_drs(dprocess) = max(n_nodes(dprocess), 1)
    n_m = nsd = number_of_smooth_drs(dprocess)

    p = SVector(model.calibration[:parameters]...)
    n_x = length(model.calibration[:controls])


    ti_trace = trace ? IterationTrace([x0]) : nothing

    # initial guess for controls
    ds0 = DRStore(Point{n_x}, [N for i=1:n_m])
    for i=1:n_m
        ds0.data[i][:] = init_dr(i, endo_nodes)
    end

    complementarities = true
    if complementarities == true
        ds_lb = DRStore([controls_lb(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m])
        ds_ub = DRStore([controls_ub(model, node(Point,dprocess,i),endo_nodes,p) for i=1:n_m])
        BIG = 100000
        for n=1:length(ds_lb.flat)
            ds_lb.flat[n] = max.(ds_lb.flat[n],-BIG)
            ds_ub.flat[n] = min.(ds_ub.flat[n], BIG)
        end
        clamp!(ds0, ds_lb, ds_ub)
    end
    ds1 = copy(ds0)


    # create decision rule (which interpolates x0)
    dr = CachedDecisionRule(dprocess, grid, [ds0.data...])

    # Define controls of tomorrow


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

    # temporary vars
    E_f = DRStore(Point{n_h}, [length(e) for e in ds0.data])
    d_E_f = deepcopy(E_f.data[1])
    S = deepcopy(s)
    X = deepcopy(ds0.data[1])


    while it<maxit && err>tol_η

        it += 1

        tic()

        set_values!(dr, [ds0.data...])

        # Compute expectations function E_f and states of tomorrow
        E_f.flat[:] *= 0.0

        for i in 1:size(E_f, 1)

            x0 = ds0.data[i]
            x1 = ds1.data[i]
            m = node(Point, dprocess, i)
            for (w, M, j) in get_integration_nodes(Point,dprocess,i)
                # Update the states
                Dolo.transition(model, Val{(0,)}, m, s, x0, M, p, (S,))
                # interpolate controles conditional states of tomorrow
                X[:] = dr(i, j, S)
                # Compute expectations as a weighted average of the exo states w_j
                Dolo.expectation(model, Val{(0,)}, M, S, X, p, (d_E_f,))
                d_E_f[:] *= w
                E_f[i][:] +=  d_E_f
            end
            # compute tomorrow's control
            Dolo.direct_response(model, Val{(0,)}, m, s, E_f[i], p, (x1,) )
        end

        if complementarities
            clamp!(ds1, ds_lb, ds_ub)
        end
        err = distance(ds0,ds1)

        copy!(ds0, ds1)

        gain = err/err_0
        err_0 = err

        elapsed = toq()

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
