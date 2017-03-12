

path = Pkg.dir("Dolo")

Pkg.build("QuantEcon")
import Dolo

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
n_exp = 0
sigma = model.calibration.flat[:sig_z]
horizon=40
seed=42
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)
s0=model.calibration[:states]
# You have to specify parameters of the function inside the brakets, strange
res = Dolo.simulation(model, sigma, dr,s0, n_exp, horizon, seed, zeros(0, 0))
res_long = Dolo.simulation(model, sigma, dr,s0, n_exp=0, horizon=100, seed=42)

res = Dolo.simulation(model, sigma, dr,s0)
res = Dolo.simulation(model, sigma)

kvec = res_long[:,2,:]
ivec = res_long[:,4,:]
nvec = res_long[:,3,:]
horizon=100
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
