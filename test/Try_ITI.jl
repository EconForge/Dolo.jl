

# Pkg.build("QuantEcon")
import Dolo
path = Dolo.pkg_path
# import Bruteforce_module

###############################################################################
## Markov Chain
###############################################################################
filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
# model = Dolo.Model(joinpath(Dolo.pkg_path, "examples", "models", "rbc_dtcc_mc.yaml"), print_code=true)
model = Dolo.yaml_import(filename)

@time dr_ITI  = Dolo.improved_time_iteration(model; verbose=true, tol = 1e-06, smaxit=100)

# dprocess = Dolo.discretize( model.exogenous )
# init_dr = Dolo.ConstantDecisionRule(model.calibration[:controls])
# @time dr_ITI_2  = Bruteforce_module.improved_time_iteration(model, dprocess,init_dr)

@time dr_TI  = Dolo.time_iteration(model;tol_η=1e-08, maxit=1000)


################################################################################
# function profile_ITI(m)
#     dr_ITI  = Bruteforce_module.improved_time_iteration(m)
#     return dr_ITI
# end
#
#
# using ProfileView
#
# profile_ITI(model)
# Profile.clear()
# @profile profile_ITI(model)
# ProfileView.view()
#



###############################################################################
#Plotting the DR

# ITIT method
df_ITI = Dolo.tabulate(model, dr_ITI.dr, :k)

import PyPlot
plt = PyPlot;
fig = plt.figure("Decision Rule, ITI-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df_ITI[:k], df_ITI[:n])
plt.ylabel("Hours");
plt.xlabel("state = k");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df_ITI[:k], df_ITI[:i])
plt.ylabel("Investment");
plt.xlabel("state = k");
plt.title("Decision Rule");
###############################################
# TI method
df = Dolo.tabulate(model, dr_TI.dr, :k)


fig = plt.figure("Decision Rule, IT-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df[:k], df[:n])
plt.ylabel("Hours");
plt.xlabel("state = k");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df[:k], df[:i])
plt.ylabel("Investment");
plt.xlabel("state = k");
plt.title("Decision Rule");

################################################################################
# RBC with an AR process
################################################################################

filename = joinpath(path,"examples","models","rbc_dtcc_ar1.yaml")
model_ar1 = Dolo.yaml_import(filename)
@time dr_ITI_ar1  = Dolo.improved_time_iteration(model_ar1; verbose=true, tol = 1e-06, smaxit=50)
@time dr_TI_ar1  = Dolo.time_iteration(model_ar1; tol_η=1e-08, maxit=1000)

dprocess = Dolo.discretize( model_ar1.exogenous , [3], [2])
@time dr2_ITI_ar1  = Dolo.improved_time_iteration(model_ar1, dprocess; verbose=true, tol = 1e-06, smaxit=50)
@time dr2_TI_ar1  = Dolo.time_iteration(model_ar1, dprocess; tol_η=1e-08, maxit=1000)


################################################################################
df_ITI = Dolo.tabulate(model_ar1, dr2_ITI_ar1.dr, :k)

import PyPlot
plt = PyPlot;
fig = plt.figure("Decision Rule, ITI-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df_ITI[:k], df_ITI[:n])
plt.ylabel("Hours");
plt.xlabel("state = k");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df_ITI[:k], df_ITI[:i])
plt.ylabel("Investment");
plt.xlabel("state = k");
plt.title("Decision Rule");


df = Dolo.tabulate(model_ar1, dr2_TI_ar1.dr, :k)


fig = plt.figure("Decision Rule, IT-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df[:k], df[:n])
plt.ylabel("Hours");
plt.xlabel("state = k");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df[:k], df[:i])
plt.ylabel("Investment");
plt.xlabel("state = k");
plt.title("Decision Rule");



################################################################################
# RBC with an IID process
################################################################################

filename = joinpath(path,"examples","models","rbc_dtcc_iid.yaml")
model_iid = Dolo.yaml_import(filename)
@time dr_ITI_iid  = Bruteforce_module.improved_time_iteration(model_iid; verbose=true, tol = 1e-06, smaxit=50)
@time dr_TI_iid  = Dolo.time_iteration(model_iid; tol_η=1e-08, maxit=1000)




df_ITI = Dolo.tabulate(model_iid, dr_ITI_iid.dr, :k)

import PyPlot
plt = PyPlot;
fig = plt.figure("Decision Rule, ITI-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df_ITI[:k], df_ITI[:n])
plt.ylabel("Hours");
plt.xlabel("state = k");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df_ITI[:k], df_ITI[:i])
plt.ylabel("Investment");
plt.xlabel("state = k");
plt.title("Decision Rule");


df = Dolo.tabulate(model_iid, dr_TI_iid.dr, :k)


fig = plt.figure("Decision Rule, IT-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df[:k], df[:n])
plt.ylabel("Hours");
plt.xlabel("state = k");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df[:k], df[:i])
plt.ylabel("Investment");
plt.xlabel("state = k");
plt.title("Decision Rule");

################################################################################
# RBC with an AR1 with 2 shocks
################################################################################

# filename = joinpath(path,"examples","models","rbc_dtcc_ar1_2shocks.yaml")
# model_ar1_2 = Dolo.yaml_import(filename)
# @time dr_ITI_ar1_2 = Bruteforce_module.improved_time_iteration(model_ar1_2; verbose=true, tol = 1e-06, smaxit=50)
# @time dr_TI_ar_2  = Dolo.time_iteration(model_ar1_2; tol_η=1e-08, maxit=1000)
#



################################################################################
# RBC with a sudden stop. Complementarities
################################################################################

filename = joinpath(path,"examples","models","sudden_stop.yaml")
model_Sudden = Dolo.yaml_import(filename)
@time dr_ITI_s  = Dolo.improved_time_iteration(model_Sudden; complementarities = true, verbose=true, tol = 1e-06, smaxit=50)
@time dr_TI_s  = Dolo.time_iteration(model_Sudden; tol_η=1e-08, maxit=1000)




df_ITI = Dolo.tabulate(model_Sudden, dr_ITI_s.dr, :l)

import PyPlot
plt = PyPlot;
fig = plt.figure("Decision Rule, ITI-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df_ITI[:l], df_ITI[:b])
plt.ylabel("b(t-1)");
plt.xlabel("state = b");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df_ITI[:l], df_ITI[:lam])
plt.ylabel("b/c");
plt.xlabel("state = b");
plt.title("Decision Rule");

#
df = Dolo.tabulate(model_Sudden, dr_TI_s.dr, :l)


fig = plt.figure("Decision Rule, IT-method",figsize=(8.5,5))

plt.subplot(1,2,1)
plt.plot(df[:l], df[:b])
plt.ylabel("Hours");
plt.xlabel("state = k");
plt.title("Decision Rule");

plt.subplots_adjust(wspace=1)
# plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(df[:l], df[:lam])
plt.ylabel("Investment");
plt.xlabel("state = k");
plt.title("Decision Rule");
