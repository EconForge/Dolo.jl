###########################################################################
# For Float

function simulate(model::AbstractModel, dr::AbstractDecisionRule, s0::AbstractVector,
                  driving_process::AbstractArray{Float64,3})

    # driving_process: (ne, N, T)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]

    N = size(driving_process, 2)
    T = size(driving_process, 3)
    epsilons = permutedims(driving_process, [2, 1, 3]) # (N, ne, T)

    # calculate initial controls using decision rule
    x0 = dr(epsilons[1, :, 1], s0)
    # get number of states and controls
    ns = length(s0)
    nx = length(x0)
    nsx = nx+ns

    s_simul = Array{Float64}(N, ns, T)
    x_simul = Array{Float64}(N, nx, T)
    for i in 1:N
      s_simul[i, :, 1] = s0
      x_simul[i, :, 1] = x0
    end

    for t in 1:T
        s = view(s_simul, :, :, t)
        m = view(epsilons, :, :, t)
        x = dr(m, s)
        x_simul[:, :, t] = x
        if t < T
          M = view(epsilons, :, :, t+1)
          ss = view(s_simul, :, :, t+1)
          transition!(model, ss, m, s, x, M, params)
        end
    end
    sim = cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}

    Ac = cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
    ll = [Symbol(i) for i in Ac]

    sim_aa = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))
    # return sim_aa
    # println(sim_aa)
    sim_def= evaluate_definitions(model, sim_aa, model.calibration[:parameters])
    return merge(sim_aa,sim_def)

end


###########################################################################
# For Int (MC)
function simulate(model::AbstractModel, dr::AbstractDecisionRule, s0::AbstractVector,
                  driving_process::AbstractMatrix{Int}, dp_process::DiscreteMarkovProcess)

    # driving_process: (ne, N, T)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]

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

    s_simul = Array{Float64}(N, ns, T)
    x_simul = Array{Float64}(N, nx, T)
    for i in 1:N
        s_simul[i, :, 1] = s0
        x_simul[i, :, 1] = x0
    end

    for t in 1:T
        s = view(s_simul, :, :, t)
        m = view(epsilons, :, :, t)

        m_ind=cat(1, m)[:, 1]
        m_val= dp_process.values[m_ind, :]
        x = dr(m_ind, s)
        x_simul[:, :, t] = x
        if t < T
            M = view(epsilons, :, :, t+1)
            M_cat = cat(1, M)[:, 1]
            M_val = dp_process.values[M_cat, :]
            ss = view(s_simul, :, :, t+1)
            transition!(model, ss, m_val, s, x, M_val, params)
        end
    end

    epsilons_values = zeros(N,nm,T)
    for n=1:N
      eps = dp_process.values[epsilons[n,:,:],:]
      epsilons_values = permutedims(eps, [1, 3, 2])
    end

    sim = cat(2, epsilons, epsilons_values, s_simul, x_simul)::Array{Float64,3}

    model_sym = :mc_process
    Ac = cat(1, model_sym, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
    ll = [Symbol(i) for i in Ac]
    sim_aa = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))
    sim_def=evaluate_definitions(model, sim_aa, model.calibration[:parameters])
    return merge(sim_aa,sim_def)

end

##############################################################################
# Sub-cases for |Floats|
function simulate(model::AbstractModel, dr::AbstractDecisionRule,
                  driving_process::AbstractArray{Float64,3})

    s0 = model.calibration[:states]
    return simulate(model, dr, s0, driving_process)
end

## methods which simulate the process
##

# Unless stochastic and return_indexes are always == true, it will work
function simulate(model::AbstractModel, dr::AbstractDecisionRule, s0::AbstractVector,
                  m0::Union{Int,AbstractVector}; N=1, T=40)
      driving_process = simulate(model.exogenous, N, T, m0)
    return simulate(model, dr, s0, driving_process)
end


function simulate(model::AbstractModel, dr::AbstractDecisionRule, s0::AbstractVector;
                  N=1, T=40)
        m0 = model.calibration[:exogenous]
        driving_process = simulate(model.exogenous, N, T, m0)
    return simulate(model, dr, s0, driving_process)
end

function simulate(model::AbstractModel,  dr::AbstractDecisionRule; N=1, T=40)
    s0 = model.calibration[:states]
    return simulate(model, dr, s0; N=N, T=T)
end
##############################################################################
# Sub-cases for |Int|

function simulate(model::AbstractModel, dr::AbstractDecisionRule,
                  driving_process::AbstractArray{Int64,2}, dp_process::Dolo.DiscreteMarkovProcess)
    s0 = model.calibration[:states]
    return simulate(model, dr, s0, driving_process, dp_process)
end

function simulate(model::AbstractModel, dr::AbstractDecisionRule, dp_process::Dolo.DiscreteMarkovProcess;
                  N::Int=1, T::Int = 40 , m0::Int = 1)
    driving_process = simulate(dp_process, N, T, m0)
    return simulate(model, dr, driving_process, dp_process)
end

function simulate(model::AbstractModel, dr::AbstractDecisionRule;
                  N::Int=1, T::Int = 40 , m0::Int = 1)
    dp_process= model.exogenous
    return simulate(model, dr, dp_process; N=N, T=T, m0 = m0)
end

function simulate(model::AbstractModel, dr::AbstractDecisionRule,
                  driving_process::AbstractArray{Int64,2};
                  N::Int=1, T::Int = 40 , m0::Int = 1)
    dp_process= model.exogenous
    return simulate(model, dr, driving_process, dp_process; N=N, T=T, m0 = m0)
end


################################################################################
## Impulse response functions

"""
 Function "response" computes the IRFs with several major options:
 - the user can provide a vector with the first values of the model's exogenous processes, e1.
 - the user can provide a name of the shock of interest and the size of the shock_name.
 - the user can provide only a name of the shock of interest. The size of the shock is assumed to be a one standard deviation given in the yaml file.

 # Arguments
 * `model::NumericModel`: Model object that describes the current model environment.
 * `dr`: Solved decision rule.
 * `e1`: the first values of the model's exogenous processes.
   or
   * `shock_name`: the name of the shock of interest.
   * `Impulse`: the size of the shock.
 * `s0`: the values of the state variables, optional.
 # Returns
 * `response`: Impulse response function.
 """


function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, e1::AbstractVector; T::Int=40)
    m_sim = response(model.exogenous, e1; T=T)
    m_simul = reshape(m_sim, 1, size(m_sim)...)
    sim = simulate(model, dr, s0, m_simul)
    sim[1, :, :] # This is now an AxisArray which seems just fine !
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  e1::AbstractVector; T::Int=40)
    s0 = model.calibration[:states]
    response(model, dr, s0, e1; T=T)
end


function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol; T::Int=40)
    index_s = findfirst(model.symbols[:exogenous], shock_name)
    # e1 = zeros(length(model.exogenous.mu))
    e1 = zeros(length(model.calibration[:exogenous]))
    Impulse = sqrt(diag(model.exogenous.Sigma)[index_s])
    e1[index_s] = Impulse
    response(model, dr, s0, e1; T=T)
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  shock_name::Symbol; kwargs...)
    s0 = model.calibration[:states]
    response(model, dr, s0, shock_name;  kwargs...)
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol, Impulse::Float64; T::Int=40)
    index_s = findfirst(model.symbols[:exogenous], shock_name)
    e1 = zeros(length(model.calibration[:exogenous]))
    e1[index_s] = Impulse
    s0 = model.calibration[:states]
    response(model, dr, s0, e1; T=T)
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  shock_name::Symbol, Impulse::Float64; T::Int=40)
    s0 = model.calibration[:states]
    response(model, dr, s0, shock_name, Impulse; T=T)
end



function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  bounds::Array{Float64,1}, s0::AbstractVector,
                  m0::Union{Int,AbstractVector}; n_steps::Int=100)

    index = findfirst(model.symbols[:states], state)
    Svalues = linspace(bounds[1], bounds[2], n_steps)
    svec = vcat([e' for e in fill(s0, n_steps)]...)
    svec[:, index] = Svalues
    m = m0  # why creating m?
    xvec = dr(m0, svec)
    mm = vcat([e' for e in fill(m, n_steps)]...)
    l1 = [mm, svec, xvec]
    tb = hcat([e' for e in l1']...)

    if isa(dr, DecisionRule{UnstructuredGrid,CartesianGrid})
        model_sym = :mc_process
    else
        model_sym = model.symbols[:exogenous]
    end

    l2 = cat(1, model_sym , model.symbols[:states], model.symbols[:controls])
    tab_AA = AxisArray(tb, Axis{state}(tb[:, index+1]), Axis{:V}(l2))
    tab_AA'
end



function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector, m0::Union{Int,AbstractVector};  n_steps=100)
    index = findfirst(model.symbols[:states], state)
    bounds = [dr.grid_endo.min[index], dr.grid_endo.max[index]]
    tabulate(model, dr, state, bounds, s0, m0;  n_steps=100)
end


function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector; n_steps=100)

    if isa(dr, DecisionRule{UnstructuredGrid,CartesianGrid})
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

#
# using PyPlot
# function plot(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
#                   bounds::Array{Float64,1}, s0::AbstractVector, m0::AbstractVector,
#                   plot_controls::Vector{Symbol};  n_steps=100)
#
#
#     df = tabulate(model, dr, state, bounds, s0, m0;  n_steps=100)
#
#     for j in plot_controls
#       fig = PyPlot.figure(j, figsize=(3, 3))
#       fig
#       PyPlot.plot(df[state], df[j], label=j)
#       PyPlot.legend()
#       # PyPlot.xlabel('state = {} | mstate = {}'.format(state, i0))
#     end
#
# end
