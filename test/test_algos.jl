path = Pkg.dir("Dolo")

import Dolo


fn = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
model_mc = Dolo.yaml_import(fn)

drc = Dolo.ConstantDecisionRule(model_mc.calibration[:controls])
# @time dr0, drv0 = Dolo.solve_policy(model_mc, drc) #, verbose=true, maxit=10000 )
@time dr = Dolo.time_iteration(model_mc, verbose=true, maxit=10000)
#
# @time drv = Dolo.evaluate_policy(model_mc, dr, verbose=true, maxit=10000)
#
# @time drd = Dolo.time_iteration_direct(model_mc, dr, verbose=true, maxit=500)
#
# # compare with prerecorded values
# kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
# nvec = [dr(1,[k])[1] for k in kvec]
# ivec = [dr(1,[k])[2] for k in kvec]
# # compare  time_iteration_direct
# nvec_d = [drd(1,[k])[1] for k in kvec]
# ivec_d = [drd(1,[k])[2] for k in kvec]
# @assert maxabs(nvec_d-nvec)<1e-4
#
# # compare  vfi
# nvec_0 = [dr0(1,[k])[1] for k in kvec]
# ivec_0 = [dr0(1,[k])[2] for k in kvec]
# @assert maxabs(nvec_0-nvec)<1e-4


# let's redo when model is stable !
# ivec_test = [0.295977,  0.257538,  0.21566,  0.173564,  0.132103,  0.0915598,  0.0520067,  0.0134661,  7.01983e-6, 3.40994e-17]
# nvec_test = [ 0.391997,  0.348033,  0.318369,  0.296276,  0.278821,  0.264487,  0.25239 ,  0.241974,  0.236604,  0.233779 ]
# @assert maximum(abs(ivec-ivec_test))<1e-5
# @assert maximum(abs(nvec-nvec_test))<1e-5


path = Pkg.dir("Dolo")
fn = joinpath(path,"examples","models","rbc_dtcc_iid.yaml")
model = Dolo.yaml_import(fn)

@time dr = Dolo.perturbate(model)
drc = Dolo.ConstantDecisionRule(model.calibration[:controls])
@time dr = Dolo.time_iteration(model, maxit=100, verbose=true)
#
# @time dr0, drv0 = Dolo.solve_policy(model, drc) #;, verbose=true, maxit=1000 )
# @time drd = Dolo.time_iteration_direct(model, maxit=1000, verbose=true)
# # @time dr = Dolo.time_iteration_direct(model, dr, maxit=500, verbose=true)
# @time drv = Dolo.evaluate_policy(model, dr, verbose=true)
#
# kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
# nvec = [dr(1,[k])[1] for k in kvec]
# ivec = [dr(1,[k])[2] for k in kvec]
# nvec_d = [drd(1,[k])[1] for k in kvec]
# ivec_d = [drd(1,[k])[2] for k in kvec]
# nvec_0 = [dr0(1,[k])[1] for k in kvec]
# ivec_0 = [dr0(1,[k])[2] for k in kvec]
#
# @assert maxabs(nvec_d-nvec)<1e-5
# @assert maxabs(nvec_0-nvec)<1e-5 # not satisfied right now (see tol. of optimizer)


# AR1 model: this one should be exactly equivalent to rbc_dtcc_ar1
fn = Pkg.dir("Dolo","examples","models","rbc_dtcc_ar1.yaml")
model = Dolo.yaml_import(fn)
dp = Dolo.discretize(model.exogenous)
@time dr = Dolo.perturbate(model)
cdr =Dolo.CachedDecisionRule(dr,dp)
@time dr = Dolo.time_iteration(model, model.exogenous, cdr)

@time dr = Dolo.time_iteration(model)
