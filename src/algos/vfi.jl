
using Optim

function Q(dmodel::DYModel, s, x, φv)

    β = discount_factor(dmodel)
    r = reward(dmodel, s, x)
    contv = sum( w*φv(S)  for (w,S) in Dolo.τ(dmodel, s, x))
    return r + β*contv

end

function Bellman_update(model, x0::GVector, φv)

    nx = deepcopy(x0)
    nv = deepcopy(φv.values)

    for (n,(s,x)) in enumerate(zip(enum( model.grid ),x0))
        
        cs = bounds(model, s)
        lb, ub = cs.min, cs.max

        x_ = max.(min.(x, ub), lb)

        res = Optim.optimize(
            u->-Q(model, s, SVector(u...), φv),
            [lb...],
            [ub...],
            Vector(x_);
        )

        nx[n] = res.minimizer
        nv[n] = -res.minimum
    end

    return nx, nv

end

function Bellman_eval(model, x0::GVector, φv)
    data = [Q(model, s, x, φv)  for (s,x) in zip(enum( model.grid ),x0)]
    GVector(model.grid,data)
end


function vfi(model::YModel; kwargs...)
    discr_options = get(kwargs, :discretization, Dict())
    dmodel = discretize(model, discr_options...)
    kwargs2 = pairs(NamedTuple( k=>v for (k,v) in kwargs if !(k in (:discretization,))))
    vfi(dmodel; kwargs2...)
end


function vfi(model::DYModel; verbose=true, improve=false, improve_wait=-1, improve_K=50, trace=false, tol_η_x=1e-8, T=1000, interpolation=:cubic)

    t0 = time_ns()

    nx = initial_guess(model)

    nv = GArray(model.grid, [1.0 for i=1:length(model.grid)])
    φv = Dolo.DFun(model.grid, nv, :value; interp_mode=interpolation)

    ti_trace = trace ? IterationTrace(typeof((;x=nx,v=φv))[]) : nothing

    log = ValueIterationLog()
    initialize(log, verbose=verbose)
    
    local kmax, η_x

    for k=1:T

        t1 = time_ns()

        fit!(φv, nv)

        trace && push!( ti_trace.data, deepcopy((;x=nx,v=φv)) )


        nx1,nv = Bellman_update(model, nx, φv)

        η_x = maximum((distance(a,b) for (a,b) in zip(nx1,nx)))

        if η_x<tol_η_x
            kmax = k
            elapsed = time_ns() - t1
            append!(log; verbose=verbose, it=k, sa_x=η_x, sa_v=NaN, time=elapsed, nit=NaN)
            break
        end
        nx = nx1

        it_eval = 0
        if improve && (k > improve_wait)
            for n=1:improve_K
                it_eval +=1
                fit!(φv, nv)
                nv = Bellman_eval(model, nx, φv)
            end
        end


        elapsed = time_ns() - t1
        append!(log; verbose=verbose, it=k, sa_x=η_x, sa_v=NaN, time=elapsed, nit=it_eval)

    end

    finalize(log, verbose=verbose)
    
    φ = Dolo.DFun(model.grid, nx; interp_mode=interpolation)

    total_time = time_ns() - t0
    
    ValueIterationResult(φ, φv, kmax,  tol_η_x, η_x, NaN, NaN, log, ti_trace)

end


