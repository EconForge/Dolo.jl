function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule,
                  s0::AbstractVector, e0::AbstractVector, driving_process::AbstractArray=zeros(0, 0);
                  n_exp::Int=1, T::Int=40,  seed::Int=42, stochastic=true)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]
    epsilons=zeros(length(model.exogenous.mu), n_exp, T)
    # simulate exogenous shocks: size (ne.N.T)
    if !isempty(driving_process)  #   "!" means NOT
      for ii in 1:size(driving_process)[2]
          epsilons[:,1,ii] = driving_process[:,ii]
      end
    else
        epsilons = simulate(model.exogenous, n_exp, T, e0; stochastic=stochastic)
    end
    epsilons = permutedims(epsilons, [2,1,3]) # (N,ne,T)
    # calculate initial controls using decision rule
    println(size(epsilons))
    x0 = dr(epsilons[1,:,1],s0)

    # get number of states and controls
    ns = length(s0)
    nx = length(x0)
    nsx = nx+ns

    # TODO: talk to Pablo and make a decision about this
    s_simul = Array(Float64, n_exp, ns, T)
    x_simul = Array(Float64, n_exp, nx, T)
    for i in 1:n_exp
      s_simul[i, :, 1] = s0
      x_simul[i, :, 1] = x0
    end
    # NOTE: this will be empty if ny is zero. That's ok. Our call to `cat`  #       below will work either way  y_simul = Array(Float64, n_exp, ny, horizon)

    for t in 1:T
        s = copy(view(s_simul, :, :, t))
        m = copy(view(epsilons, :, :, t))
        x = dr(m,s)
        x_simul[:, :, t] = x
        if t < T
          M = view(epsilons, :, :, t+1)
          ss = view(s_simul, :, :, t+1)
          println([size(e) for e in [ss,m,s,x,M]])
          Dolo.transition!(model, (ss), (m), (s), (x), (M), params)
          # s_simul[:, :, t+1] = ss
        end
    end
    if !isempty(driving_process)
          return cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}
    else
          out = cat(2, epsilons, s_simul, x_simul)
          out = out[1,:,:]
          columns = cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
          return DataFrame(Dict(columns[i]=>out[i,:] for i=1:length(columns)))
    end
end

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule, s0::AbstractVector,
                  driving_process::AbstractArray=zeros(0, 0); kwargs...)
    e0 = model.calibration[:exogenous]
    return simulate(model, dr, s0, e0, driving_process; kwargs...)
end

function simulate(model::AbstractNumericModel,  dr::AbstractDecisionRule; kwargs...)
    s0 = model.calibration[:states]
    e0 = model.calibration[:exogenous]
    return simulate(model, dr, s0, e0; kwargs...)
end

using DataFrames

function response(model::AbstractNumericModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol, Impulse::Float64;
                  T::Integer=40)
    index_s = findfirst(model.symbols[:exogenous], shock_name)
    m_simul = zeros(length(model.exogenous.mu), T)
    m_simul[index_s,:] = response(T, Impulse)

    sims = simulate(model, dr, s0, m_simul, stochastic = false; n_exp=1, T=T)
    sim = sims[1,:,:]
    columns = cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
    return DataFrame(Dict(columns[i]=>sim[i,:] for i=1:length(columns)))
end

function response(model::AbstractNumericModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol;
                  T::Integer=40)
    index_s = findfirst(model.symbols[:exogenous], shock_name)
    Impulse = zeros(0)
    if isempty(Impulse)
      Impulse = sqrt(diag(model.exogenous.Sigma)[index_s])
    end
    return response(model, dr, s0, shock_name, Impulse; T=T)
end


function response(model::AbstractNumericModel,  dr::AbstractDecisionRule,
                  shock_name::Symbol;
                  kwargs...)
    s0 = model.calibration[:states]
    return response(model,  dr, s0, shock_name; kwargs...)
end
