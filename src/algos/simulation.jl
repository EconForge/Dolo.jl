function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule,
                           s0::AbstractVector, e0::AbstractVector;
                           n_exp::Int=1, horizon::Int=40,  seed::Int=42, stochastic=true)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]

    # simulate exogenous shocks: size (ne.N.T)
    epsilons = simulate(model.exogenous, n_exp, horizon, e0; stochastic=stochastic)
    epsilons = permutedims(epsilons, [2,1,3]) # (N,ne,T)

    # calculate initial controls using decision rule
    println(size(epsilons))
    x0 = dr(epsilons[1,:,1],s0)

    # get number of states and controls
    ns = length(s0)
    nx = length(x0)
    nsx = nx+ns

    # TODO: talk to Pablo and make a decision about this
    s_simul = Array(Float64, n_exp, ns, horizon)
    x_simul = Array(Float64, n_exp, nx, horizon)
    for i in 1:n_exp
      s_simul[i, :, 1] = s0
      x_simul[i, :, 1] = x0
    end
    # NOTE: this will be empty if ny is zero. That's ok. Our call to `cat`  #       below will work either way  y_simul = Array(Float64, n_exp, ny, horizon)

    for t in 1:horizon
        s = copy(view(s_simul, :, :, t))
        m = copy(view(epsilons, :, :, t))
        x = dr(m,s)
        x_simul[:, :, t] = x
        if t < horizon
          M = view(epsilons, :, :, t+1)
          ss = view(s_simul, :, :, t+1)
          println([size(e) for e in [ss,m,s,x,M]])
          ss = Dolo.transition!(model, (ss), (m), (s), (x), (M), params)
          s_simul[:, :, t+1] = ss
        end
    end

    return cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}

end

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule, s0::AbstractVector; kwargs...)
    e0 = model.calibration[:exogenous]
    return simulate(model, dr, s0, e0; kwargs...)
end

function simulate(model::AbstractNumericModel,  dr::AbstractDecisionRule; kwargs...)
    s0 = model.calibration[:states]
    e0 = model.calibration[:exogenous]
    return simulate(model, dr, s0, e0; kwargs...)
end

using DataFrames

function response(model::AbstractNumericModel,  dr::AbstractDecisionRule, e0::AbstractVector; kwargs...)
    s0 = model.calibration[:states]
    sims = simulate(model, dr, s0, e0; stochastic=false, n_exp=1, kwargs...)
    sim = sims[1,:,:]
    columns = cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
    return DataFrame(Dict(columns[i]=>sim[i,:] for i=1:length(columns)))
end
