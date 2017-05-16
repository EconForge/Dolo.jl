
path = Pkg.dir("Dolo")

# Pkg.build("QuantEcon")
import Dolo

using AxisArrays
using Unitful
import Unitful: s, ms, µs



filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)


# model.exogenous
# m0 = model.calibration[:exogenous]
# #############################################################3
# ### Starting a function, simulate indexes
# n_states = size(model.exogenous.values,1)
#
# model.exogenous.values
#
#
#
# horizon = 50
# n_exp=10
# simul = zeros(Int, horizon, n_exp)
# simul[1,:] = 1
# simul
# # rnd = random.rand(horizon* n_exp).reshape((horizon,n_exp))
# rnd = rand(horizon, n_exp)
#
# cumuls = cumsum(model.exogenous.transitions, 2)
#
# function choice(x, n, cumul)
#     i = 1
#     running = true
#     # while (i<n) && running
#     while (i<n) && running
#         if x < cumul[i]
#             running = false
#         else
#             i += 1
#         end
#     end
#     return i
# end
#
#
# for t in 1:(horizon-1)
#   for j in 1:n_exp
#     s = simul[t,j]
#     p = cumuls[s,:]
#     v = choice(rnd[t,j], n_states, p)
#     simul[t+1,j] = v
#   end
#   println(t)
# end
#
# simul
###############   function ##############
function choice(x, n, cumul)
    i = 1
    running = true
    # while (i<n) && running
    while (i<n) && running
        if x < cumul[i]
            running = false
        else
            i += 1
        end
    end
    return i
end

function simulate(nodes::Array{Float64,2}, transitions::Array{Float64,2},N::Int, T::Int; i0::Int=1)
      n_states = size(nodes,1)

      simul = zeros(Int, T, N)
      simul[1,:] = 1
      rnd = rand(T, N)
      cumuls = cumsum(transitions, 2)

      for t in 1:(T-1)
        for j in 1:N
          s = simul[t,j]
          p = cumuls[s,:]
          v = choice(rnd[t,j], n_states, p)
          simul[t+1,j] = v
        end
      end

      return simul

end


###############################################################################
### Simulate shocks
# Simulate indexes
N=10
T= 50
Index_mc = simulate(model.exogenous.values, model.exogenous.transitions, 10, 50;i0=1)
ind2 = zeros(Int,1,N,T)
for i in 1:N
  ind2[1,i,:]=Index_mc[:,i]
end
ind2[1,1,:]==Index_mc[:,1]
Index_mc::AbstractArray
ind2
s0=model.calibration[:states]
Dolo.simulate(model, dr, s0, ind2)




calib = model.calibration
params = calib[:parameters]

epsilons = permutedims(ind2, [2,1,3]) # (N,ne,T)
N = size(epsilons,1)
T = size(epsilons,3)
S0=model.calibration[:states]
# calculate initial controls using decision rule
x0 = dr(epsilons[1,:,1],S0)
s0=S0
ns = length(s0)
nx = length(x0)
nsx = nx+ns

s_simul = Array(Float64, N, ns, T)
x_simul = Array(Float64, N, nx, T)
for i in 1:N
  s_simul[i, :, 1] = s0
  x_simul[i, :, 1] = x0
end


# t=1
# s = view(s_simul, :, :, t)
# m = view(epsilons, :, :, t)
# m_cat=cat(1,m)[:,1]
# x = dr(m_cat,s)
# # s_cat = cat(1,s)
#
# # m_cat
# typeof(dr)
# println()
# x = vcat( [dr(m_cat[j], s_cat[j,:])' for j=1:size(s_cat,1)]... )
#
# m_val= model.exogenous.values[m_cat,:]
#
# x_simul[:, :, t] = x
# s
# x
# ss
# t=1
#
# M = view(epsilons, :, :, t+1)
# M_cat= cat(1,M)[:,1]
# M_val= model.exogenous.values[M_cat,:]
# ss = view(s_simul, :, :, t+1)
# Dolo.transition!(model, (ss), (m_val), (s), (x), (M_val), params)


for t in 1:T
    s = view(s_simul, :, :, t)
    m = view(epsilons, :, :, t)
    m_cat=cat(1,m)[:,1]
    m_val= model.exogenous.values[m_cat,:]
    x = dr(m_cat,s)
    x_simul[:, :, t] = x
    if t < T
      M = view(epsilons, :, :, t+1)
      M_cat= cat(1,M)[:,1]
      M_val= model.exogenous.values[M_cat,:]
      ss = view(s_simul, :, :, t+1)
      Dolo.transition!(model, (ss), (m_val), (s), (x), (M_val), params)
      # s_simul[:, :, t+1] = ss
    end
end
epsilons
s_simul
x_simul
sim = cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}

model_sym=model.symbols[:exogenous]
model_sym=:mc_process


Ac= cat(1, model_sym, model.symbols[:states], model.symbols[:controls])
ll=[Symbol(i) for i in Ac]
AA = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))
N
T=50







Dolo.simulate(model, dr, S0, ind2)


Shocks= model.exogenous.values[Index_mc]
Shocks[:,1]
S2 = zeros(1,N,T)
for i in 1:N
  S2[1,i,:]=Shocks[:,i]
end
S2[1,1,:]==Shocks[:,1]



function simulate(process::Dolo.DiscreteMarkovProcess, N::Int, T::Int; i0=1)
      Index_mc = simulate(process.values, process.transitions, N, T; i0=1)
      Shocks = model.exogenous.values[Index_mc]
      S2 = zeros(1,N,T)
      for i in 1:N
        S2[1,i,:]=Shocks[:,i]
      end
      return S2
end


Shocks_mc= simulate(model.exogenous, 10,50)

s0=model.calibration[:states]
S0::AbstractVector


Dolo.simulate(model, dr, S0, Shocks_mc)
model::Dolo.AbstractNumericModel
dr:: Dolo.AbstractDecisionRule

####################################################################################

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model2 = Dolo.yaml_import(filename)
@time dr2 = Dolo.time_iteration(model2, verbose=true, maxit=10000)

s0=model2.calibration[:states]
m0 = model2.calibration[:exogenous]
driving_process = Dolo.simulate(model2.exogenous, 10, 40, m0)

# permutedims(driving_process, [2,1,3])
Dolo.simulate(model2, dr2, driving_process)



# driving_process: (ne,N,T)

# extract data from model
calib = model2.calibration
params = calib[:parameters]

epsilons = permutedims(driving_process, [2,1,3]) # (N,ne,T)
N = size(epsilons,1)
T = size(epsilons,3)

# calculate initial controls using decision rule
x0 = dr2(epsilons[1,:,1],s0)
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

t=1
s = view(s_simul, :, :, t)
m = view(epsilons, :, :, t)

if
    m_cat=cat(1,m)[:,1]
    m_val= model2.exogenous.values[m_cat,:]
end

x = dr2(m,s)
x_simul[:, :, t] = x
if t < T
  M = view(epsilons, :, :, t+1)
  M_cat= cat(1,M)[:,1]
  M_val= model2.exogenous.values[M_cat,:]
  ss = view(s_simul, :, :, t+1)
  Dolo.transition!(model2, (ss), (m_val), (s), (x), (M_val), params)
  # s_simul[:, :, t+1] = ss
end



if t < T
  M = view(epsilons, :, :, t+1)
  ss = view(s_simul, :, :, t+1)
  Dolo.transition!(model2, (ss), (m), (s), (x), (M), params)
end




for t in 1:T
    s = view(s_simul, :, :, t)
    m = view(epsilons, :, :, t)
    x = dr2(m,s)
    x_simul[:, :, t] = x
    if t < T
      M = view(epsilons, :, :, t+1)
      ss = view(s_simul, :, :, t+1)
      Dolo.transition!(model2, (ss), (m), (s), (x), (M), params)
      # s_simul[:, :, t+1] = ss
    end
end

sim = cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}
Ac= cat(1, model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls])
ll=[Symbol(i) for i in Ac]
AA = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))
