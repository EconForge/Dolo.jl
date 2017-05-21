using DataFrames
using PyPlot
using AxisArrays

function simulate(model::AbstractModel, dr::AbstractDecisionRule,
                  s0::AbstractVector, driving_process::Union{Array{Int64,2},Array{Float64,3}})

    # driving_process: (ne,N,T)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]

    if typeof(driving_process)==Array{Int64,2}
          N = size(driving_process,2)
          T = size(driving_process,1)
          epsilons = zeros(Int,1,N,T)
          for i in 1:N
            epsilons[1,i,:]=driving_process[:,i]
          end
    else
      epsilons = driving_process
      N = size(driving_process,2)
      T = size(driving_process,3)
    end


    epsilons = permutedims(epsilons, [2,1,3]) # (N,ne,T)

    # calculate initial controls using decision rule
    x0 = dr(epsilons[1,:,1],s0)
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
        if typeof(driving_process)==Array{Float64,3}
            m_val = m
            x = dr(m,s)
        else
            m_ind=cat(1,m)[:,1]
            m_val= model.exogenous.values[m_ind,:]
            x = dr(m_ind,s)
        end
        x_simul[:, :, t] = x
        if t < T
          M = view(epsilons, :, :, t+1)
          if typeof(driving_process)==Array{Float64,3}
              M_val = M
          else
              M_cat=cat(1,M)[:,1]
              M_val= model.exogenous.values[M_cat,:]
          end
          ss = view(s_simul, :, :, t+1)
          Dolo.transition!(model, (ss), (m_val), (s), (x), (M_val), params)
          # s_simul[:, :, t+1] = ss
        end
    end

    sim = cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}

    if typeof(driving_process)==Array{Float64,3}
        model_sym=model.symbols[:exogenous]
    else
        model_sym=:mc_process
    end

    Ac= cat(1, model_sym, model.symbols[:states], model.symbols[:controls])
    ll=[Symbol(i) for i in Ac]
    AA = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))

    return AA
end



function simulate(model::AbstractModel, dr::AbstractDecisionRule,
                  driving_process::Union{Array{Int64,2},Array{Float64,3}})
    s0 = model.calibration[:states]
    return simulate(model, dr, s0, driving_process)
end



##
## methods which simulate the process
##

# Unless stochastic and return_indexes are always == true, it will work
function simulate(model::AbstractModel, dr::AbstractDecisionRule, s0::AbstractVector,
                  m0::Union{Int,AbstractVector}; N=1, T=40, option=true)
    driving_process = simulate(model.exogenous, N, T, m0)
    return simulate(model, dr, s0, driving_process)
end

function simulate(model::AbstractModel, dr::AbstractDecisionRule, s0::AbstractVector;
                  N=1, T=40, option=true)
    if typeof(dr)==Dolo.DecisionRule{Dolo.UnstructuredGrid,Dolo.CartesianGrid}
        m0 = 1
    else
        m0 = model.calibration[:exogenous]
    end
    driving_process = simulate(model.exogenous, N, T, m0)
    return simulate(model, dr, s0, driving_process)
end

function simulate(model::AbstractModel,  dr::AbstractDecisionRule; kwargs...)
    s0 = model.calibration[:states]
    # m0 = model.calibration[:exogenous]
    return simulate(model, dr, s0; kwargs...)
end

##
## Impulse response functions
##

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, e1::AbstractVector; T::Integer=40)
    m_sim=response(model.exogenous, e1; T=T)
    m_simul=reshape(m_sim, 1, size(m_sim)...)
    sim = simulate(model, dr, s0, m_simul)
    return sim[1,:,:] # This is now an AxisArray which seems just fine !
end

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol; T::Integer=40)
    index_s = findfirst(model.symbols[:exogenous], shock_name)
    e1 = zeros(length(model.exogenous.mu))
    Impulse = sqrt(diag(model.exogenous.Sigma)[index_s])
    e1[index_s] = Impulse
    return response(model, dr, s0, e1; T=T)
end

import DataFrames

function response(model::AbstractModel,  dr::AbstractDecisionRule,
                  shock_name::Symbol; kwargs...)
    s0 = model.calibration[:states]
    return response(model, dr, s0, shock_name;  kwargs...)
end


function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  bounds::Array{Float64,1}, s0::AbstractVector, m0::AbstractVector;  n_steps=100)

    index = findfirst(model.symbols[:states],state)
    Svalues = linspace(bounds[1], bounds[2], n_steps)
    svec = vcat([e' for e in fill(s0, n_steps)]...)
    svec[:,index] = Svalues
    m = m0  # why creating m?
    xvec = dr(m0,svec)
    mm = vcat([e' for e in fill(m, n_steps)]...)
    l1 = [mm, svec, xvec]
    tb = hcat([e' for e in l1']...)

    l2 = cat(1, model.symbols[:exogenous] , model.symbols[:states] , model.symbols[:controls] )
    ll=[string(i) for i in l2]
    d = OrderedDict()
    for i=1:length(ll)
        d[ll[i]] = tb[:,i]
    end

    df = DataFrames.DataFrame(d)

    return df

end



function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector, m0::AbstractVector;  n_steps=100)
    index = findfirst(model.symbols[:states],:z)
    bounds = [dr.grid_endo.min[index], dr.grid_endo.max[index]]
    df = tabulate(model, dr, state,bounds, s0, m0;  n_steps=100)
    return df

end


function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector;  n_steps=100)
    m0 = model.calibration[:exogenous]
    df = tabulate(model, dr, state, s0;  n_steps=100)
    return df
end

function tabulate(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol;  n_steps=100)
    s0 = model.calibration[:states]
    df = tabulate(model, dr, state, s0;  n_steps=100)
    return df
end



function plot(model::AbstractModel, dr::AbstractDecisionRule, state::Symbol,
                  bounds::Array{Float64,1}, s0::AbstractVector, m0::AbstractVector,
                  plot_controls::Vector{Symbol};  n_steps=100)


    df = tabulate(model, dr, state,bounds, s0, m0;  n_steps=100)

    for j in plot_controls
      fig = PyPlot.figure(j,figsize=(3,3))
      fig
      PyPlot.plot(df[state], df[j], label=j)
      PyPlot.legend()
      # PyPlot.xlabel('state = {} | mstate = {}'.format(state, i0))
    end

end
