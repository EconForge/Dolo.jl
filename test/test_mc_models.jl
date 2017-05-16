
path = Pkg.dir("Dolo")

# Pkg.build("QuantEcon")
import Dolo

using AxisArrays
using Unitful
import Unitful: s, ms, µs



filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)


model.exogenous
m0 = model.calibration[:exogenous]
#############################################################3
### Starting a function, simulate indexes
n_states = size(model.exogenous.values,1)

model.exogenous.values



horizon = 50
n_exp=10
simul = zeros(Int, horizon, n_exp)
simul[1,:] = 1
simul
# rnd = random.rand(horizon* n_exp).reshape((horizon,n_exp))
rnd = rand(horizon, n_exp)

cumuls = cumsum(model.exogenous.transitions, 2)

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


for t in 1:(horizon-1)
  for j in 1:n_exp
    s = simul[t,j]
    p = cumuls[s,:]
    v = choice(rnd[t,j], n_states, p)
    simul[t+1,j] = v
  end
  println(t)
end

simul
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
Index_mc = simulate(model.exogenous.values, model.exogenous.transitions, 10, 50;i0=1)
ind2 = zeros(1,N,T)
for i in 1:N
  ind2[1,i,:]=Index_mc[:,i]
end
ind2[1,1,:]==Index_mc[:,1]
Index_mc::AbstractArray

Dolo.simulate(model, dr, Index_mc)


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


Dolo.simulate(model, dr, Shocks_mc)
model::Dolo.AbstractNumericModel
dr:: Dolo.AbstractDecisionRule


filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model2 = Dolo.yaml_import(filename)
@time dr2 = Dolo.time_iteration(model2, verbose=true, maxit=10000)

m0 = model2.calibration[:exogenous]
driving_process = Dolo.simulate(model2.exogenous, 10, 40, m0)

permutedims(driving_process, [2,1,3])
Dolo.simulate(model2, dr2, driving_process)
