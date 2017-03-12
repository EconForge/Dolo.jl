# function QuantEcon.simulate_2(m::AbstractNumericModel, dr::AbstractDecisionRule,
#                            s0::AbstractVector=m.calibration[:states];
#                            n_exp::Int=0, horizon::Int=40,
#                            seed::Int=42,
#                            forcing_shocks::AbstractMatrix=zeros(0, 0))
#
function simulation(model::AbstractNumericModel, dr::AbstractDecisionRule,
                           s0::AbstractVector;
                           n_exp::Int=1, horizon::Int=40,  seed::Int=42, stochastic=true)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]

    # simulate exogenous shocks: size (ne.N.T)
    epsilons = simulate(model.exogenous, n_exp, horizon; stochastic=stochastic)
    epsilons = permutedims(epsilons, [2,1,3])

    # calculate initial controls using decision rule

    println(s0)
    println(epsilons[:,1,1])
    x0 = dr(epsilons[:,1,1],s0)

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
        m = copy(view(epsilons,:,:,t))
        x = dr(m,s)
        x_simul[:, :, t] = x
        if t < horizon
          M = view(epsilons,:,:,t+1)
          ss = view(s_simul, :, :, t+1)
          ss = Dolo.transition!(model, vec(ss), vec(m), vec(s), vec(x), vec(M), params)
          s_simul[:, :, t+1] = ss
        end
    end

    return cat(2, s_simul, x_simul)::Array{Float64,3}

end


function simulation(model::AbstractNumericModel,  dr::AbstractDecisionRule; kwargs...)
    s0 = model.calibration[:states]
    return simulation(model, dr, s0; kwargs...)
end
