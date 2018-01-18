"""
Document pf.
"""
function perfect_foresight(model, exo::AbstractMatrix{Float64}; T=200, verbose=true, complementarities=true)

    T = T+1 # compat with t=0 conventions

    T_e = size(exo, 1)
    n_e = length(model.symbols[:exogenous])
    n_s = length(model.symbols[:states])
    n_x = length(model.symbols[:controls])

    driving_process = zeros(T, n_e)
    driving_process[1:T_e,:] = exo
    for t=(T_e+1):T
        driving_process[t,:] = exo[end,:]
    end

    # find initial steady-state

    s0 = model.calibration[:states]
    x0 = model.calibration[:controls]
    p0 = model.calibration[:parameters]
    m0 = driving_process[1,:]

    function res_ss(u, m0)
        cc = Dict(:exogenous=>m0,
                  :parameters=>p0)
        cc[:states] = u[1:n_s]
        cc[:controls] = u[n_s+1:end]
        res = Dolo.residuals(model, cc)

        return [res[:transition]; res[:arbitrage]]
    end

    u0 = [s0; x0]

    sol_init = NLsolve.nlsolve(not_in_place(u->res_ss(u,m0)), u0)
    if ~NLsolve.converged(sol_init)
        error("Couldn't find initial guess.")
    end
    us = sol_init.zero
    s0 = us[1:n_s]
    x0 = us[(n_s+1):end]

    # construct dumb initial guess
    initial_guess = zeros(T, n_s+n_x)
    for t=1:T
        initial_guess[t,:] = [s0 ;x0]
    end

    function residuals(model, s0, driving_process, vv)

        n_s = length(model.symbols[:states])
        n_x = length(model.symbols[:controls])
        n_e = length(model.symbols[:exogenous])
        p = model.calibration[:parameters]
        T = size(vv,1)
        ss = vv[:,1:n_s]       # x0 ... xT
        xx = vv[:,(n_s+1):end] # s0 ... sT
        mm = driving_process

        s_bar = ss[1,:]

        ss_p = [s_bar'; ss[1:end-1,:]]      # s_bar, s0,   ..., s_{T-1}
        mm_p = [mm[1,:]'; mm[1:end-1,:]]
        xx_p = [xx[1,:]'; xx[1:end-1,:]]
        mm_f = [mm[2:end,:]; mm[end,:]']  # s1, ... s_{T-1}, s_T, s_T
        ss_f = [ss[2:end,:]; ss[end,:]']  # s1, ... s_{T-1}, s_T, s_T
        xx_f = [xx[2:end,:]; xx[end,:]']

        G = Dolo.transition(model, mm_p,ss_p,xx_p,mm,p) - ss
        F = Dolo.arbitrage(model, mm,ss,xx,mm_f,ss_f,xx_f,p)

        # special cases
        G[1,:] = ss[1,1:n_s] - s0
        F[end,:] = Dolo.arbitrage(model,
                    mm[end,:],
                    ss[end,:],
                    xx[end,:],
                    mm[end,:],
                    ss[end,:],
                    xx[end,:],
                    p)

        R = [G F]
        # res = zeros(T, 2)
        return R

    end

    sh = size(initial_guess)

    fun = u->residuals(model,s0,driving_process,reshape(u, sh...))[:]

    vv0 = initial_guess[:]

    if ~complementarities
        sol = NLsolve.nlsolve(not_in_place(fun), vv0, show_trace=verbose)
    else
        ss0 = initial_guess[:, 1:n_s]
        mm0 = driving_process
        lb = [ss0*0-Inf Dolo.controls_lb(model, mm0, ss0, p0)]
        ub = [ss0*0+Inf Dolo.controls_ub(model, mm0, ss0, p0)]

        R0 = fun(vv0)
        sol = NLsolve.mcpsolve(not_in_place(fun), lb[:], ub[:], vv0, show_trace=verbose)
    end

    if ~NLsolve.converged(sol)
        error("Stacked-time system couldn't be solved.")
    end

    headers = [model.symbols[:exogenous]; model.symbols[:states]; model.symbols[:controls]]
    resp = [driving_process reshape(sol.zero, sh...)]

    out = AxisArray(resp', Axis{:V}(headers) ,Axis{:T}(0:T-1))
    defs = evaluate_definitions(model, out)
    merge(out, defs)

end

# Constructs the matrix exo given a dictionary exo and calls the original method
function perfect_foresight(model, exo::Dict{}; kwargs... )
    convert(Dict{Symbol,Array{Float64,1}}, exo)
    n_e = length(model.symbols[:exogenous])
    T_e = maximum(length(e) for e in values(exo))
    exo_new = zeros(T_e, n_e)
    for (i, key) in enumerate(model.symbols[:exogenous])
        exo_new[:, i] = model.calibration.flat[key]
    end

    for key in keys(exo)
        ind = find(model.symbols[:exogenous] .== key)
        T_key = length(exo[key])
        exo_new[1:T_key,ind] = exo[key]

        for t=(T_key+1):T_e
            exo_new[t,ind] =  exo[key][end]
        end

    end

    exo = exo_new

    return perfect_foresight(model, exo; kwargs... )

end
