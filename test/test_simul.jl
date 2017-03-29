
path = Pkg.dir("Dolo")

Pkg.build("QuantEcon")
import Dolo

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
# model2 = Dolo.yaml_import(filename2)
N = 1
T=40
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)

s0 = model.calibration[:states]
e0 = model.calibration[:exogenous]
irf = Dolo.response(model, dr, s0, :e_z, 0.03)
Dolo.response(model, dr, s0, :e_z)
e0=[5]
irf = true
if irf==true
   e0=[4];
end
e0
############## 2 shocks
filename2 = joinpath(path,"examples","models","rbc_dtcc_iid_2ar1.yaml")
model2 = Dolo.yaml_import(filename2)
@time dr2 = Dolo.time_iteration(model2, verbose=true, maxit=10000)
s0 = model2.calibration[:states]
e0 = model2.calibration[:exogenous]

# index_s = findfirst(model.symbols[:exogenous], :e_z)

irf=Dolo.response(model2, dr2, s0, :e_d, 0.3; T=40)
irf2=Dolo.response(model2, dr2, s0, :e_z)

Dolo.response(model2, dr2, s0, :e_d)
# IRF_2 =Dolo.response(model, dr, e0, :e_z)

#########################################################################
# simulation AR1
var = Dolo.VAR1([0.0,0.0],[0.99 0.0; 0.0 0.08],eye(2)*0.01)
n = size(var.M, 1);
n
var.Sigma
zeros(size(var.M, 1))
[0.1, 0]
sqrt(diag(var.Sigma)[1])
Dolo.response(var,1,[0.1,0])


srand(123) # Setting the seed


d = Dolo.MvNormal(var.Sigma);
N =1
T=40
VAR_process=Dolo.simulate(var, N, T)
w= permutedims(VAR_process, [2,1,3])
w[2]
E = VAR_process;
E[:, :, 1]
repeat(mean(E[:, :, 1], 2))
mean(E[:, :, 1])
XN = VAR_process[1];


# Computing moments
Mean_sim = mean(squeeze(mean(VAR_process, 3), 2), 2);
squeeze(mean(VAR_process, 3), 2)
mean(VAR_process, 3)
Mean_sim
#Computing the std of simulted Processes (Covariance matrix)
E1_d = (E[:, :, 1] - repeat(mean(E[:, :, 1], 2), inner=[1, T]));
cov(E1_d[:, 1])    # which is X1_d[:, 1]'*X1_d[:, 1]/T
diag(cov(E1_d[:, :]))   # covariances across simulation
Sigma_sim = mean(diag(cov(E1_d[:, :])) )  # mean of covariances across simulations
# Autocorrelation matrix
X_d=zeros(N, T, 2)
X_d0=zeros(N, T-1, 2);
X_d1=zeros(N, T-1, 2);
R_sim = zeros(2, 2);
for ii in [1, 2]
X_d[:, :, ii] = (XN[:, :, 2] - repeat(mean(XN[:, :, ii], 2), inner=[1, T]));
X_d0[:, :, ii]  = X_d[:, 1:end-1, ii];
X_d1[:, :, ii] = (XN[:, 2:end, ii] - repeat(mean(XN[:, 2:end, ii], 2), inner=[1, T-1]));
R_sim[ii, ii] =  mean(diag( cov(X_d0[:, :, ii], X_d1[:, :, ii] )/sqrt(var(X_d0[:, :, ii]))/sqrt(var(X_d1[:, :, ii]))  ))
  if ii == 2
      R_sim[ii-1, ii]  =   mean(diag( cov(X_d0[:, :, ii-1], X_d1[:, :, ii] )/sqrt(var(X_d0[:, :, ii-1]))/sqrt(var(X_d1[:, :, ii])) ))
      R_sim[ii, ii-1]  =   mean(diag( cov(X_d0[:, :, ii], X_d1[:, :, ii-1] )/sqrt(var(X_d0[:, :, ii]))/sqrt(var(X_d1[:, :, ii-1])) ))
  end
end

return Mean_sim, Sigma_sim, R_sim




















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
