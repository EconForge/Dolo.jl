function evaluate_definitions(model, simul::AxisArray{Tf,3}, params=model.calibration[:parameters]) where Tf

    p_ = SVector(params...)

    # @assert axisnames(simul) == (:N,:V,:T)
    T = length( simul[Axis{:T}].val )
    @assert simul[Axis{:T}].val == 1:T

    vars = cat(model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls]; dims=1)

    sim = simul[Axis{:V}(vars)]

    T = length(sim[Axis{:T}].val)

    past = permutedims(sim[Axis{:T}([1;1:T-1])], [2,1,3])
    present = permutedims(sim[Axis{:T}(1:T)], [2,1,3])
    future = permutedims(sim[Axis{:T}([2:T;T])], [2,1,3])

    n_v,N,T = size(past)

    x_past = reshape(reinterpret(Point{n_v}, vec(past.data)), (T*N,))
    x_present = reshape(reinterpret(Point{n_v}, vec(present.data)), (T*N,))
    x_future = reshape(reinterpret(Point{n_v}, vec(future.data)), (T*N,))

    y_ = evaluate_definitions(model, x_past, x_present, x_future, p_)

    auxiliaries = [Dolang.arg_name(e) for e in keys(model.definitions)]
    n_y = length(auxiliaries)

    data = permutedims( reshape(reinterpret(Float64, vec(y_)), (n_y,N,T)), [2,1,3])

    array = AxisArray(copy(data), Axis{:N}(1:N), Axis{:V}(auxiliaries), Axis{:T}(1:T))

end

function evaluate_definitions(model, _simul::AxisArray{__T,2}, params=model.calibration[:parameters]) where __T

    p_ = SVector(params...)

    all_v = _simul[Axis{:V}].val
    all_t = _simul[Axis{:T}].val

    # TODO: why are we even doing that ?
    T = length(all_t)
    if all(all_t .== 0:(T-1))
        all_t = 1:T
    end
    if all(all_t .!= 1:T)
        msg = "Can only evaluate definitions if the T Axis goes from 1:T or 0:(T-1)"
        error(msg)
    end

    # expect dimensions to be (:V,:T) otherwise we transpose
    if axisnames(_simul) == (:T,:V)
        data = _simul.data'
    else
        data = _simul.data
    end

    simul = AxisArray(data, Axis{:V}(all_v), Axis{:T}(all_t))

    vars = cat(model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls]; dims=1)

    sim = simul[Axis{:V}(vars)]

    T = length(sim[Axis{:T}].val)

    past = sim[Axis{:T}([1; 1:T-1])]
    present = sim[Axis{:T}(1:T)]
    future = sim[Axis{:T}([2:T; T])]

    n_v = size(past, 1)

    x_past = reshape(reinterpret(Point{n_v}, vec(past.data)), (T,))
    x_present = reshape(reinterpret(Point{n_v}, vec(present.data)), (T,))
    x_future = reshape(reinterpret(Point{n_v}, vec(future.data)), (T,))

    y_ = evaluate_definitions(model, x_past, x_present, x_future, p_)

    auxiliaries = [Dolang.arg_name(e) for e in keys(model.definitions)]
    n_y = length(auxiliaries)

    data = reshape(reinterpret(Float64, vec(y_)), (n_y, T))

    array = AxisArray(copy(data), Axis{:V}(auxiliaries), simul[Axis{:T}])

end

########################################################################
# For Float

function simulate(model::AbstractModel, dr::AbstractDecisionRule,
                  driving_process::AbstractArray{Float64,3}; s0::AbstractVector=model.calibration[:states])

    # driving_process: (ne, N, T)
    driving_process = convert(Array{Float64,3}, driving_process) # in case arg is an axisarray

    # extract data from model
    calib = model.calibration
    params = SVector(calib[:parameters]...)

    N = size(driving_process, 2)
    T = size(driving_process, 3)
    epsilons = permutedims(driving_process, [2, 1, 3]) # (N, ne, T)
    # calculate initial controls using decision rule

    x0 = dr(epsilons[1, :, 1], s0)

    # get number of states and controls
    ns = length(s0)
    nx = length(x0)
    nsx = nx+ns

    s_simul = Array{Float64}(undef, N, ns, T)
    x_simul = Array{Float64}(undef, N, nx, T)
    for i in 1:N
      s_simul[i, :, 1] = s0
      x_simul[i, :, 1] = x0
    end


    for t in 1:T
        s = copy(to_LOP(s_simul[:, :, t]))
        m = copy(to_LOP(epsilons[:, :, t]))
        x = dr(m, s)
        x_simul[:, :, t] = from_LOP(x)
        if t < T
          M = copy(to_LOP(epsilons[:, :, t+1]))
          ss = transition(model, m, s, x, M, params)
          s_simul[:,:,t+1] = from_LOP(ss)
        end
    end
    sim = cat(epsilons, s_simul, x_simul; dims=2)::Array{Float64,3}

    Ac = cat(model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls]; dims=1)
    ll = [Symbol(i) for i in Ac]

    sim_aa = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))

    if length(model.definitions)>0
        sim_def= evaluate_definitions(model, sim_aa, model.calibration[:parameters])
        return merge(sim_aa,sim_def)
    else
        return sim_aa
    end

end


###########################################################################
# For Int (MC)
function simulate(model::AbstractModel, dr::AbstractDecisionRule,
                  driving_process::AbstractMatrix{Int}, dp_process::DiscreteMarkovProcess;
                  s0::AbstractVector=model.calibration[:states])

    # driving_process: (ne, N, T)

    # extract data from model
    # calib = model.calibration
    params = model.calibration[:parameters]

    N = size(driving_process, 2)
    T = size(driving_process, 1)
    epsilons = zeros(Int, N, 1, T)
    for _n in 1:N
        epsilons[_n, 1, :] = driving_process[:, _n]
    end

    # calculate initial controls using decision rule
    x0 = dr(epsilons[1, :, 1], s0)
    # get number of states and controls
    ns = length(s0)
    nx = length(x0)
    nsx = nx+ns
    nm = length(model.calibration[:exogenous])

    s_simul = Array{Float64}(undef, N, ns, T)
    x_simul = Array{Float64}(undef, N, nx, T)
    for i in 1:N
        s_simul[i, :, 1] = s0
        x_simul[i, :, 1] = x0
    end

    for t in 1:T
        s = s_simul[:, :, t]
        m = epsilons[:, :, t]

        m_ind=cat(m; dims=1)[:, 1]
        m_val= dp_process.values[m_ind, :]
        x = dr(m_ind, s)
        x_simul[:, :, t] = x
        if t < T
            M = epsilons[:, :, t+1]
            M_cat = cat(M, dims=1)[:, 1]
            M_val = dp_process.values[M_cat, :]
            s_simul[:,:,t+1] = transition(model, m_val, s, x, M_val, params)
        end
    end

    epsilons_values = zeros(N,nm,T)
    for n=1:N
      eps = dp_process.values[epsilons[n,:,:],:]
      epsilons_values[n,:,:] = permutedims(eps, [1, 3, 2])
    end

    sim = cat(epsilons, epsilons_values, s_simul, x_simul; dims=2)::Array{Float64,3}

    model_sym = :mc_process
    Ac = cat(model_sym, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls], dims=1)
    ll = [Symbol(i) for i in Ac]
    sim_aa = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))

    if length(model.definitions)>0
        sim_def=evaluate_definitions(model, sim_aa, model.calibration[:parameters])
        return merge(sim_aa,sim_def)
    else
        return sim_aa
    end

end


"""
This function simulates a model given a decision rule.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dr`: Solved decision rule.
* optional:
  * `driving process` simulated series of a model's exogenous process.
  If  `driving process` is not provided, then:
  * `dprocess`: exogenous processes; default: model.exogenous.
  * `s0::ListOfPoints`: List of initial state variable values; default: model.calibration[:states]
  * `N`: number of simulations; default: 1.
  * `T`: number of periods of simulations; default: 40.
  * `i0`: initial state of a discretized exogenous process; default: default_index(dprocess).
  or
  * `m0::ListOfPoints`: List of initial values of a continuous exogenous process; default: model.calibration[:exogenous].
# Returns
* `simulate`: simulated time series.
"""
function simulate(model::AbstractModel, dr::AbstractDecisionRule, dprocess::DiscreteMarkovProcess;
                  i0::Int=default_index(dprocess), s0::AbstractVector=model.calibration[:states], N::Int=1, T::Int=40)
    driving_process = simulate(dprocess; N=N, T=T, i0=i0)
    return simulate(model, dr, driving_process, dprocess; s0=s0)
end


function simulate(model::AbstractModel, dr::AbstractDecisionRule, dprocess::ContinuousProcess;
                  m0::AbstractVector=model.calibration[:exogenous], s0::AbstractVector=model.calibration[:states],
                  N::Int=1, T::Int=40)
    driving_process = simulate(dprocess, m0; N=N, T=T)
    return simulate(model, dr, driving_process; s0 = s0)
end


function simulate(model::AbstractModel, dr::AbstractDecisionRule; kwargs...)
    dprocess = model.exogenous
    return simulate(model, dr, dprocess; kwargs...)
end


"""
This is the one we document.
"""
simulate
################################################################################
## Impulse response functions


function response(model::AbstractModel,  dr::AbstractDecisionRule;
                  s0::AbstractVector, e1::AbstractVector, T::Int=40)
    m_sim = response(model.exogenous, e1; T=T)
    m_simul = reshape(m_sim, size(m_sim,1), 1, size(m_sim,2))
    sim = simulate(model, dr, m_simul; s0=s0)
    sim[1, :, :] # This is now an AxisArray which seems just fine !
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
    s0::AbstractVector, e1::AbstractVector; T::Int=40)

    response(model,  dr; s0=s0, e1=e1, T=T)
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  e1::AbstractVector; T::Int=40)
    s0 = model.calibration[:states]
    response(model, dr; s0=s0, e1=e1, T=T)
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol; T::Int=40)
    index_s = something(findfirst(isequal(shock_name), model.symbols[:exogenous]), 0)
    # e1 = zeros(length(model.exogenous.mu))
    e1 = zeros(length(model.calibration[:exogenous]))
    Impulse = sqrt(diag(model.exogenous.Σ)[index_s])
    e1[index_s] = Impulse
    response(model, dr; s0=s0, e1=e1, T=T)
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  shock_name::Symbol; kwargs...)
    s0 = model.calibration[:states]
    response(model, dr, s0, shock_name;  kwargs...)
end


"""
Function "response" computes the impulse response functions with several major options:
- the user can provide a vector with the first values of the model's exogenous processes, e1.
- the user can provide a name of the shock of interest and the size of the shock_name.
- the user can provide only a name of the shock of interest. The size of the shock is assumed to be a one standard deviation given in the yaml file.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `dr`: Solved decision rule.
* `e1::ListOfPoints`: List of initial model's exogenous processes values.
* If e1 is not provided, then:
  * `shock_name`: the name of the shock of interest.
* optional:
  * `Impulse`: the size of the shock; default: one standard deviation.
  * `s0::ListOfPoints`: List of initial state variable values; default: model.calibration[:states]
# Returns
* `response`: Impulse response function.
"""
function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol, Impulse::Float64; T::Int=40)
    index_s = something(findfirst(isequal(shock_name), model.symbols[:exogenous]), 0)
    e1 = zeros(length(model.calibration[:exogenous]))
    e1[index_s] = Impulse
    s0 = model.calibration[:states]
    response(model, dr; s0=s0, e1=e1, T=T)
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  shock_name::Symbol, Impulse::Float64; T::Int=40)
    s0 = model.calibration[:states]
    response(model, dr, s0, shock_name, Impulse; T=T)
end


"""
Function "tabulate" produces a 2-dimensional AxisArray{Float64,2,...}  with 2 axes : V (containing all the variables of the model) and the name of the state variable chosen in input. 
"""
function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  bounds::Array{Float64,1}, s0::AbstractVector,
                  m0::Union{Int,AbstractVector}; n_steps::Int=100)

    index = findfirst(isequal(state), model.symbols[:states])
    Svalues = range(bounds[1], stop=bounds[2], length=n_steps)
    svec = vcat([e' for e in fill(s0, n_steps)]...)
    svec[:, index] = Svalues

    if isa(dr, AbstractDecisionRule{UnstructuredGrid,UCGrid})
        model_sym = :mc_process
    else
        xvec = dr(m0, svec)
        mm = vcat([e' for e in fill(m0, n_steps)]...)
        l1 = [mm, svec, xvec]
        tb = hcat([e' for e in l1']...)
        model_sym = model.symbols[:exogenous]
    end

    l2 = cat(model_sym , model.symbols[:states], model.symbols[:controls]; dims=1)
    tab_AA = AxisArray(reshape(tb,1,size(tb)...), Axis{:T}(1:1), Axis{:N}(1:n_steps), Axis{:V}(l2))

    ## add definitions
    tab_AAA = permutedims(tab_AA, [2,3,1])

    if length(model.definitions)>0
        tab_defs = evaluate_definitions(model, tab_AAA)
        tab_ = merge(tab_AAA, tab_defs)[Axis{:T}(1)]
    else
        tab_ = tab_AAA[Axis{:T}(1)]
    end
    #change axis names
    res = AxisArray(tab_.data, Axis{state}(tb[:, index+1]), tab_[Axis{:V}])
    res' # so that we can index it directly
end

function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector, m0::Union{Int,AbstractVector};  n_steps=100)
    index = findfirst(isequal(state), model.symbols[:states])
    bounds = [model.domain.d2.min[index], model.domain.d2.max[index]]
    tabulate(model, dr, state, bounds, s0, m0;  n_steps=100)
end


function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector; n_steps=100)

    if isa(dr, AbstractDecisionRule{UnstructuredGrid,UCGrid})
        m0 = 1
    else
        m0 = model.calibration[:exogenous]
    end
    tabulate(model, dr, state, s0, m0;  n_steps=100)
end


function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol; n_steps=100)
    s0 = model.calibration[:states]
    tabulate(model, dr, state, s0;  n_steps=100)
end

"""
That's how we tabulate functions.
"""
