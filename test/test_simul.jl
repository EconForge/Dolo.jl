
path = Pkg.dir("Dolo")

Pkg.build("QuantEcon")
import Dolo

filename = joinpath(path,"examples","models","rbc_dtcc_iid_2ar1.yaml")
# filename2 = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
# model2 = Dolo.yaml_import(filename2)

e0 = model.calibration[:exogenous]
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
# mvn2 = model2.exogenous
# mvn2.Sigma
# diag(mvn.Sigma)[ind_shock]
# mvn.mu[ind_shock]
dist = Distributions.MvNormal([mvn.mu[ind_shock]], diag(mvn.Sigma)[ind_shock])
# dist = Distributions.MvNormal(mvn.mu, diag(mvn.Sigma))
# dist2 = Distributions.MvNormal(mvn2.mu, diag(mvn2.Sigma))
d = length(mvn.mu)
out = zeros(d, N, T)
# x::Symbol = :e_z



for i=1:N
    out[ind_shock,i,1] = e0
end
out

rand(dist, N)
for t =2:T
    out[ind_shock,:,t] = rand(dist, N)
end
out


irf = true
if irf
    out[ind_shock,:,2] = diag(sqrt(mvn.Sigma))[ind_shock]
      for t =3:T
              out[:,:,t] = 0
      end
else
      for t =2:T
          out[ind_shock,:,t] = 0
      end
end

out


f























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
