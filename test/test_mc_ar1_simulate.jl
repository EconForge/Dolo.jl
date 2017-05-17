
path = Pkg.dir("Dolo")

# Pkg.build("QuantEcon")
import Dolo

using AxisArrays
using Unitful
import Unitful: s, ms, µs



filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)

###############################################################################
### Simulate a model with a Markov Chain

N=10
T= 50
Index_mc = Dolo.simulate(model.exogenous, N, T, 1)
# ind2 = zeros(Int,1,N,T)
# for i in 1:N
#   ind2[1,i,:]=Index_mc[:,i]
# end

s0=model.calibration[:states]

sim = Dolo.simulate(model, dr, s0, Index_mc)

Tvalues = linspace(1, sim[Axis{:T}][end], sim[Axis{:T}][end])

a_sort = sort(collect(sim[Axis{:V}(:k)]),1)
a_median = a_sort[Int(size(sim,1)*0.5),:]


import PyPlot
plt = PyPlot;
fig = plt.figure("Markov Chaine")
for i in 1:sim[Axis{:N}][end]
    plt.plot(Tvalues, sim[Axis{:N}(i), Axis{:V}(:k)], color="red", alpha = .05)
end
plt.title("Simulations");
plt.xlabel("Horizon");
plt.ylabel("Capital");


plt.plot(Tvalues, a_median, marker="o")
#
#
# calib = model.calibration
# params = calib[:parameters]
#
# epsilons = permutedims(ind2, [2,1,3]) # (N,ne,T)
# N = size(epsilons,1)
# T = size(epsilons,3)
# S0=model.calibration[:states]
# # calculate initial controls using decision rule
# x0 = dr(epsilons[1,:,1],S0)
# s0=S0
# ns = length(s0)
# nx = length(x0)
# nsx = nx+ns
#
# s_simul = Array(Float64, N, ns, T)
# x_simul = Array(Float64, N, nx, T)
# for i in 1:N
#   s_simul[i, :, 1] = s0
#   x_simul[i, :, 1] = x0
# end
#
#
# t=1
# s = view(s_simul, :, :, t)
# m = view(epsilons, :, :, t)
# m_cat=cat(1,m)[:,1]
# x = dr(m_cat,s)
# s_cat = cat(1,s)
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
#
#
# for t in 1:T
#     s = view(s_simul, :, :, t)
#     m = view(epsilons, :, :, t)
#     m_cat=cat(1,m)[:,1]
#     m_val= model.exogenous.values[m_cat,:]
#     x = dr(m_cat,s)
#     x_simul[:, :, t] = x
#     if t < T
#       M = view(epsilons, :, :, t+1)
#       M_cat= cat(1,M)[:,1]
#       M_val= model.exogenous.values[M_cat,:]
#       ss = view(s_simul, :, :, t+1)
#       Dolo.transition!(model, (ss), (m_val), (s), (x), (M_val), params)
#       # s_simul[:, :, t+1] = ss
#     end
# end
# epsilons
# s_simul
# x_simul
# sim = cat(2, epsilons, s_simul, x_simul)::Array{Float64,3}
#
# model_sym=model.symbols[:exogenous]
# model_sym=:mc_process
#
#
# Ac= cat(1, model_sym, model.symbols[:states], model.symbols[:controls])
# ll=[Symbol(i) for i in Ac]
# AA = AxisArray(sim, Axis{:N}(1:N), Axis{:V}(ll), Axis{:T}(1:T))
# N
# T=50
#
# Shocks= model.exogenous.values[Index_mc]
# Shocks[:,1]
# S2 = zeros(1,N,T)
# for i in 1:N
#   S2[1,i,:]=Shocks[:,i]
# end
# S2[1,1,:]==Shocks[:,1]

# function simulate(process::Dolo.DiscreteMarkovProcess, N::Int, T::Int; i0=1)
#       Index_mc = simulate(process.values, process.transitions, N, T; i0=1)
#       Shocks = model.exogenous.values[Index_mc]
#       S2 = zeros(1,N,T)
#       for i in 1:N
#         S2[1,i,:]=Shocks[:,i]
#       end
#       return S2
# end


Shocks_mc= Dolo.simulate(model.exogenous, 10,50, 1; return_indexes=false)

s0=model.calibration[:states]

# Check it won't work with values
Dolo.simulate(model, dr, s0, Shocks_mc)

# Check it works without providing States
sim_mc2 = Dolo.simulate(model, dr, Index_mc)
sim_mc2==sim

# Check it works without providing a driving process
sim_mc3 = Dolo.simulate(model, dr, s0; N=10, T=50)


# Check it works without providing a driving process
sim_mc4 = Dolo.simulate(model, dr, s0, 1; N=10,T=50)

import PyPlot
plt = PyPlot;
fig = plt.figure("Markov Chaine check")
for i in 1:sim_mc4[Axis{:N}][end]
    plt.plot(Tvalues, sim_mc4[Axis{:N}(i), Axis{:V}(:k)], color="red", alpha = .05)
end
plt.title("Simulations");
plt.xlabel("Horizon");
plt.ylabel("Capital");

# Check Simulation work for AR1 process

####################################################################################

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model2 = Dolo.yaml_import(filename)
@time dr2 = Dolo.time_iteration(model2, verbose=true, maxit=10000)

s0=model2.calibration[:states]
m0 = model2.calibration[:exogenous]
driving_process = Dolo.simulate(model2.exogenous, N, T, m0)

# permutedims(driving_process, [2,1,3])
sim_ar = Dolo.simulate(model2, dr2, driving_process)

Tvalues = linspace(1, sim_ar[Axis{:T}][end], sim_ar[Axis{:T}][end])

a_sort = sort(collect(sim_ar[Axis{:V}(:k)]),1)
a_median = a_sort[Int(size(sim_ar,1)*0.5),:]


import PyPlot
plt = PyPlot;
fig = plt.figure("Simulate a model with AR1")
for i in 1:sim_ar[Axis{:N}][end]
    plt.plot(Tvalues, sim_ar[Axis{:N}(i), Axis{:V}(:k)], color="red", alpha = .05)
end
plt.title("Simulations");
plt.xlabel("Horizon");
plt.ylabel("Capital");


plt.plot(Tvalues, a_median, marker="o")

sim_ar2 = Dolo.simulate(model2, dr2, driving_process)
sim_ar2==sim_ar

# Check it works without providing a driving process
sim_ar3 = Dolo.simulate(model2, dr2, s0)


# Check it works without providing a driving process
sim_ar4 = Dolo.simulate(model2, dr2, s0, m0; N=100)
Tvalues = linspace(1, sim_ar4[Axis{:T}][end], sim_ar4[Axis{:T}][end])


import PyPlot
plt = PyPlot;
fig = plt.figure("AR1 check")
for i in 1:sim_ar4[Axis{:N}][end]
    plt.plot(Tvalues, sim_ar4[Axis{:N}(i), Axis{:V}(:k)], color="red", alpha = .05)
end
plt.title("Simulations");
plt.xlabel("Horizon");
plt.ylabel("Capital");
