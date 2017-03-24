
path = Pkg.dir("Dolo")

Pkg.build("QuantEcon")
import Dolo

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
# model2 = Dolo.yaml_import(filename2)
N = 1
T=40
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)


filename2 = joinpath(path,"examples","models","rbc_dtcc_iid_2ar1.yaml")
model2 = Dolo.yaml_import(filename2)
@time dr2 = Dolo.time_iteration(model2, verbose=true, maxit=10000)
s0 = model2.calibration[:states]
e0 = model2.calibration[:exogenous]

# index_s = findfirst(model.symbols[:exogenous], :e_z)

irf=Dolo.response(model2, dr2, s0, e0, :e_d, 0.3; T=40)
irf2=Dolo.response(model2, dr2, s0, e0, :e_z, 0.3; T=40)

Dolo.response(model2, dr2, s0, e0, :e_d)
# IRF_2 =Dolo.response(model, dr, e0, :e_z)






sims = Dolo.simulate(model2, dr2, s0, e0, stochastic = true; n_exp=1, T=40 )  






Impulse = 0.3
index_s = findfirst(model2.symbols[:exogenous], :e_d)
isempty(Impulse)
if isempty(Impulse)
  Impulse = sqrt(diag(model2.exogenous.Sigma)[index_s])
end
Impulse

isempty(Impulse)
if !isempty(Impulse)
  Impulse = impulse_response
else
  Impulse = sqrt(diag(model.exogenous.Sigma)[index_s])
end










### Simulation model




# calculate initial controls using decision rule
println(size(epsilons))
x0 = dr(epsilons[1,:,1],s0)

# get number of states and controls
calib = model.calibration
params = calib[:parameters]
driving_process = m_simul
driving_process
!isempty(driving_process)
epsilons=zeros(length(model.exogenous.mu), 1, T)
if !isempty(driving_process)
    for ii in 1:size(driving_process)[2]
      epsilons[:,1,ii] = driving_process[:,ii]
    end
else
      epsilons = Dolo.simulate(model.exogenous, 1, T, e0; stochastic=true)
end
epsilons

ns = length(s0)
nx = length(x0)
nsx = nx+ns
n_exp=1
T=40
# TODO: talk to Pablo and make a decision about this
s_simul = Array(Float64, n_exp, ns, T)
x_simul = Array(Float64, n_exp, nx, T)
for i in 1:n_exp
  s_simul[i, :, 1] = s0
  x_simul[i, :, 1] = x0
end
# NOTE: this will be empty if ny is zero. That's ok. Our call to `cat`  #       below will work either way  y_simul = Array(Float64, n_exp, ny, horizon)
epsilons

t=1
s = copy(view(s_simul, :, :, t))
m = copy(view(epsilons, :, :, t))
x = dr(m,s)
x_simul[:, :, t] = x
t
M = view(epsilons, :, :, t+1)
ss = view(s_simul, :, :, t+1)
println([size(e) for e in [ss,m,s,x,M]])
ss = Dolo.transition!(model, (ss), (m), (s), (x), (M), params)
s_simul[:, :, t+1] = ss
























stochastic = false
irf=true
model
model.exogenous
epsilons = Dolo.simulate(model, model.exogenous, 1, 40, e0, :e_d; stochastic=stochastic, irf=irf)
epsilons = Dolo.simulate(model, model.exogenous, 1, 40, e0, :e_d, [0.03]; stochastic=stochastic, irf=irf)


# stochastic=true
# e0 = model.calibration[:exogenous]
# n_exp=1
# horizon=40
model.exogenous.Sigma
sqrt(model.exogenous.Sigma)

# # simulate exogenous shocks: size (ne.N.T)
# stochastic = false
# epsilons = Dolo.simulate(model.exogenous, n_exp, horizon, e0; stochastic=stochastic)
# epsilons = Dolo.simulate(model.exogenous, n_exp, horizon, e0; stochastic=false, irf=true)
# epsilons = permutedims(epsilons, [2,1,3]) # (N,ne,T)
# irf = n_exp == 1 ? true : false


N = 1
T=40
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)


ind_shock = findfirst(model.symbols[:exogenous], :e_d)
e0 = model.calibration[:exogenous][ind_shock]

 mvn = model.exogenous
mvn2 = model2.exogenous
xx= diag(mvn.Sigma)[ind_shock]
typeof(xx)
# mvn2.Sigma




















irf = Dolo.response(model, dr, e0; horizon = 100)

kirf = irf[:,3]
iirf = irf[:,2]
nirf = irf[:,4]
zirf = irf[:,5]
horizon=100
time = linspace(0,horizon-1,horizon)
using Gadfly

plot(x=time, y=kirf, Geom.point, Geom.line,
     Guide.xlabel("horizon"), Guide.ylabel("Capital"), Guide.title("IRF"))
plot(x=time, y=nirf, Geom.point, Geom.line,Guide.xlabel("horizon"),
     Guide.ylabel("Hours"), Guide.title("IRF"))
plot(x=time, y=iirf, Geom.point, Geom.line, Guide.xlabel("horizon"),
   Guide.ylabel("Investments"), Guide.title("IRF"))
plot(x=time, y=zirf, Geom.point, Geom.line, Guide.xlabel("horizon"),
      Guide.ylabel("AR1"), Guide.title("IRF"))
















# test simulation
res = Dolo.simulation(model, sigma, dr,s0, n_exp, horizon, seed, zeros(0, 0))
res_long = Dolo.simulation(model, sigma, dr,s0, n_exp=0, horizon=1000, seed=42)

res = Dolo.simulation(model, sigma, dr,s0)
res = Dolo.simulation(model, sigma)

kvec = res_long[:,2,:]
ivec = res_long[:,4,:]
nvec = res_long[:,3,:]
horizon=1000
time = linspace(0,horizon-1,horizon)
using Gadfly

plot(x=time, y=kvec, Geom.point, Geom.line,
     Guide.xlabel("horizon"), Guide.ylabel("Capital"), Guide.title("Simulations"))
plot(x=time, y=nvec, Geom.point, Geom.line,Guide.xlabel("horizon"), Guide.ylabel("Hours"), Guide.title("Simulations"))
plot(x=time, y=ivec, Geom.point, Geom.line, Guide.xlabel("horizon"), Guide.ylabel("Investments"), Guide.title("Simulations"))


O = Theme(
    default_color=colorant"orange")

plot(layer(x=time, y=nvec, Geom.point, Geom.line, O),
     layer(x=time, y=ivec, Geom.point, Geom.line),
     Guide.xlabel("horizon"),
     Guide.ylabel("Hours (orange), Investment"),
     Guide.title("Simulations"))


# verbose=true
# verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n" "t" model.symbols[:states][1] model.symbols[:states][2] model.symbols[:controls][1] model.symbols[:controls][2]
# verbose && println(repeat("-", 35))


#verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n"  t round(s[1],2) round(s[2],2) round(x[1],2) round(x[2],2)
