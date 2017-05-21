using DataFrames
using PyPlot
using AxisArrays




function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule,
                       s0::AbstractVector, driving_process)
                      #  s0::AbstractVector, driving_process::Array{Float64,3})

    # driving_process: (ne,N,T)

    # extract data from model
    calib = model.calibration
    params = calib[:parameters]

    epsilons = permutedims(driving_process, [2,1,3]) # (N,ne,T)
    N = size(epsilons,1)
    T = size(epsilons,3)

    # calculate initial controls using decision rule
    x0 = dr(epsilons[1,:,1],s0)
    # get number of states and controls
    ns = length(s0)
    nx = length(x0)
    nsx = nx+ns

    s_simul = Array(Float64, N, ns, T)
    x_simul = Array(Float64, N, nx, T)
    for i in 1:N
      s_simul[i, :, 1] = s0
      x_simul[i, :, 1] = x0
    end

    for t in 1:T
        s = view(s_simul, :, :, t)
        m = view(epsilons, :, :, t)
        x = dr(m,s)
        x_simul[:, :, t] = x
        if t < T
          M = view(epsilons, :, :, t+1)
          ss = view(s_simul, :, :, t+1)
          Dolo.transition!(model, (ss), (m), (s), (x), (M), params)
          # s_simul[:, :, t+1] = ss
        end
    end

    sim = cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}
    Ac= cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
    ll=[Symbol(i) for i in Ac]
    AA = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))

    return AA

end

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule,
                  s0::AbstractVector, driving_process::Union{Array{Float64,2},Array{Float64,3}})

    # driving_process: (T,ne)
    epsilons = reshape(driving_process, 1, size(driving_process)...) # (N,T,ne)
    epsilons = permutedims(epsilons,[3,1,2])
    sim = simulate(model, dr, s0, epsilons)
    out = sim[1,:,:]
    columns = cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
    return DataFrame(
            merge(
                Dict(:t=>0:(size(out,2)-1)),
                Dict(columns[i]=>out[i,:] for i=1:length(columns))
            )
        )
end

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule,
                  driving_process::Union{Array{Float64,2},Array{Float64,3}})
    s0 = model.calibration[:states]
    return simulate(model, dr, s0, driving_process)
end

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule,
                  driving_process::Array{Int64,3})
    s0 = model.calibration[:states]
    return simulate(model, dr, s0, driving_process)
end

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule,
                  s0::AbstractVector, driving_process::Array{Int64,3})

    # driving_process: (T,ne)
    epsilons = reshape(driving_process, 1, size(driving_process)...) # (N,T,ne)
    epsilons = permutedims(epsilons,[3,1,2])
    sim = simulate(model, dr, s0, epsilons)
    out = sim[1,:,:]
    columns = cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
    return DataFrame(
            merge(
                Dict(:t=>0:(size(out,2)-1)),
                Dict(columns[i]=>out[i,:] for i=1:length(columns))
            )
        )
end

##
## methods which simulate the process
##

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule, s0::AbstractVector, m0::AbstractVector; N=1, T=40, stochastic=true)
                  driving_process = simulate(model.exogenous, N, T, m0; stochastic=stochastic)
    return simulate(model, dr, s0, driving_process)
end

function simulate(model::AbstractNumericModel, dr::AbstractDecisionRule, s0::AbstractVector;
                  N=1, T=40, stochastic=true)
    m0 = model.calibration[:exogenous]
    driving_process = simulate(model.exogenous, N, T, m0; stochastic=stochastic)
    return simulate(model, dr, s0, driving_process)
end

function simulate(model::AbstractNumericModel,  dr::AbstractDecisionRule; kwargs...)
    s0 = model.calibration[:states]
    m0 = model.calibration[:exogenous]
    return simulate(model, dr, s0, m0; kwargs...)
end

##
## Impulse response functions
##

function response(model::AbstractNumericModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, e1::AbstractVector; T::Integer=40)
    m_simul=response(model.exogenous, e1; T=T)
    if typeof(m_simul)==Array{Float64,3}
      m_simul=(vec(m_simul)')'
    else
      m_simul = m_simul'
    end
    sim = simulate(model, dr, s0, m_simul)
    return sim
end

function response(model::AbstractNumericModel,  dr::AbstractDecisionRule,
                  s0::AbstractVector, shock_name::Symbol; T::Integer=40)
    index_s = findfirst(model.symbols[:exogenous], shock_name)
    e1 = zeros(length(model.exogenous.mu))
    Impulse = sqrt(diag(model.exogenous.Sigma)[index_s])
    e1[index_s] = Impulse
    return response(model, dr, s0, e1; T=T)
end

import DataFrames

function response(model::AbstractNumericModel,  dr::AbstractDecisionRule,
                  shock_name::Symbol; kwargs...)
    s0 = model.calibration[:states]
    return response(model, dr, s0, shock_name;  kwargs...)
end


function tabulate(model::AbstractNumericModel, dr::AbstractDecisionRule, state::Symbol,
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



function tabulate(model::AbstractNumericModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector, m0::AbstractVector;  n_steps=100)
    index = findfirst(model.symbols[:states],:z)
    bounds = [dr.grid_endo.min[index], dr.grid_endo.max[index]]
    df = tabulate(model, dr, state,bounds, s0, m0;  n_steps=100)
    return df

end


function tabulate(model::AbstractNumericModel, dr::AbstractDecisionRule, state::Symbol,
                  s0::AbstractVector;  n_steps=100)
    m0 = model.calibration[:exogenous]
    df = tabulate(model, dr, state, s0;  n_steps=100)
    return df
end

function tabulate(model::AbstractNumericModel, dr::AbstractDecisionRule, state::Symbol;  n_steps=100)
    s0 = model.calibration[:states]
    df = tabulate(model, dr, state, s0;  n_steps=100)
    return df
end



function plot(model::AbstractNumericModel, dr::AbstractDecisionRule, state::Symbol,
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
