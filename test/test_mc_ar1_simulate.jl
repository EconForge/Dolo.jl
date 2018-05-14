

# Pkg.build("QuantEcon")
import Dolo

path = Dolo.pkg_path
using AxisArrays

filename = joinpath(path,"examples","models","rbc_mc.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)
@time dr = Dolo.improved_time_iteration(model, verbose=true, maxit=10000)

N=10
T= 50
Index_mc = Dolo.simulate(model.exogenous, N, T, 1)

s0=model.calibration[:states]
s0=repmat(s0,N)
typeof(s0)<:AbstractVector
Index_mc.data

dp_process=model.exogenous

sim = Dolo.simulate(model, dr.dr, Index_mc, model.exogenous)
sim = Dolo.simulate(model, dr.dr; i0=2)
sim = Dolo.simulate(model, dr.dr)

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

Shocks_mc= Dolo.simulate(model.exogenous, 10,50, 1; return_indexes=false)

s0=model.calibration[:states]

# Check it won't work with values
# Dolo.simulate(model, dr, s0, Shocks_mc)
# Check


# Check it works without providing States
sim_mc2 = Dolo.simulate(model, dr, Index_mc, model.exogenous)
sim_mc2==sim

# Check it works without providing a driving process
sim_mc3 = Dolo.simulate(model, dr, model.exogenous; N=10, T=50, m0=2)

# Check it works without providing a driving process
sim_mc4 = Dolo.simulate(model, dr, s0, model.exogenous)
Tvalues = linspace(1, sim_mc4[Axis{:T}][end], sim_mc4[Axis{:T}][end])

import PyPlot
plt = PyPlot;
fig = plt.figure("Markov Chaine check")
for i in 1:sim_mc4[Axis{:N}][end]
    plt.plot(Tvalues, sim_mc4[Axis{:N}(i), Axis{:V}(:k)], color="red", alpha = .05)
end
plt.title("Simulations");
plt.xlabel("Horizon");
plt.ylabel("Capital");


################################################################################
# Tabulation
s0=model.calibration[:states]
m0=1

Dolo.tabulate(model, dr, :k, s0, m0)
bounds = [4.0, 10.0]

Dolo.tabulate(model, dr, :k, bounds, model.calibration[:states], m0)

# Check Simulation work for AR1 process

####################################################################################

filename = joinpath(path,"examples","models","rbc.yaml")
model2 = Dolo.yaml_import(filename)
@time dr2 = Dolo.time_iteration(model2, verbose=true, maxit=10000)

driving_process = Dolo.simulate(model2.exogenous, N, T)
dprocess = model2.exogenous
sim_ar = Dolo.simulate(model2, dr2.dr, driving_process)
sim_ar = Dolo.simulate(model2, dr2.dr, dprocess; m0=[0.1], N=20)
sim_ar = Dolo.simulate(model2, dr2.dr)


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
s0=model2.calibration[:states]
sim_ar3 = Dolo.simulate(model2, dr2, s0)

# Check it works without providing a driving process
m0 = model2.calibration[:exogenous]
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

################################################################################
# Tabulation
s0=model2.calibration[:states]
m0 = model2.calibration[:exogenous]

Tab= Dolo.tabulate(model2, dr2, :k, s0, m0)
Dolo.tabulate(model2, dr2, :k)


Tab[:n]

# Check Simulation work for AR1 process when discretized with MC

####################################################################################


filename = joinpath(path,"examples","models","rbc.yaml")
model3 = Dolo.yaml_import(filename)

n_states=5
mc_ar = Dolo.discretize_mc(model3.exogenous;N=n_states)

@time dr_armc = Dolo.time_iteration(model3, mc_ar, verbose=true, maxit=10000)

# If you solve a model with a continuous exogenous process using MC disretization, then to simulate you need first proved a driving process originated form this new MC. model.exogenous
# At this point there is no option in the code. If you do not specify any driving porcess(series of shocks), then Dolo takes model.exogenous to do it itself, but then it can't combine dr(i,s) with values of shocks. Shall it be added?
Index_mc = Dolo.simulate(mc_ar, N, T, 1)

s0=model3.calibration[:states]
sim_armc = Dolo.simulate(model3, dr_armc, s0, Index_mc, mc_ar)
sim_armc_2 = Dolo.simulate(model3, dr_armc, s0, mc_ar)
sim_armc_3 = Dolo.simulate(model3, dr_armc, mc_ar)
