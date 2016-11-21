path = Pkg.dir("Dolo")


import Dolo

filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
# filename = joinpath(path,"examples","models","sudden_stop.yaml")
model_mc = Dolo.yaml_import(filename)

drc = Dolo.ConstantDecisionRule(model_mc.calibration[:controls])
@time dr0, drv0 = Dolo.solve_policy(model_mc, drc, verbose=true, maxit=10000 )


@time dr = Dolo.time_iteration(model_mc, verbose=true, maxit=10000)
@time drv = Dolo.evaluate_policy(model_mc, dr, verbose=true, maxit=10000)

@time drd = Dolo.time_iteration_direct(model_mc, dr, verbose=true, maxit=500)

# compare with prerecorded values
kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
nvec = [Dolo.evaluate(dr,1,[k])[1] for k in kvec]
ivec = [Dolo.evaluate(dr,1,[k])[2] for k in kvec]
# compare  time_iteration_direct
nvec_d = [Dolo.evaluate(drd,1,[k])[1] for k in kvec]
ivec_d = [Dolo.evaluate(drd,1,[k])[2] for k in kvec]
@assert maximum(abs(nvec_d-nvec))<1e-

# compare  vfi
nvec_0 = [Dolo.evaluate(dr0,1,[k])[1] for k in kvec]
ivec_0 = [Dolo.evaluate(dr0,1,[k])[2] for k in kvec]
@assert maximum(abs(nvec_0-nvec))<1e-4


# let's redo when model is stable !
# ivec_test = [0.295977,  0.257538,  0.21566,  0.173564,  0.132103,  0.0915598,  0.0520067,  0.0134661,  7.01983e-6, 3.40994e-17]
# nvec_test = [ 0.391997,  0.348033,  0.318369,  0.296276,  0.278821,  0.264487,  0.25239 ,  0.241974,  0.236604,  0.233779 ]
# @assert maximum(abs(ivec-ivec_test))<1e-5
# @assert maximum(abs(nvec-nvec_test))<1e-5


# this one needs a lower value of beta or a better initial guess
filename = joinpath(path,"examples","models","rbc_dtcc_iid.yaml")
model = Dolo.yaml_import(filename)

drc = Dolo.ConstantDecisionRule(model.calibration[:controls])
@time dr0, drv0 = Dolo.solve_policy(model, drc, verbose=true, maxit=1000 )


@time dr = Dolo.time_iteration(model, maxit=1000, verbose=true)
@time drd = Dolo.time_iteration_direct(model, maxit=1000, verbose=true)
# @time dr = Dolo.time_iteration_direct(model, dr, maxit=500, verbose=true)
@time drv = Dolo.evaluate_policy(model, dr, verbose=true)

kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
nvec = [Dolo.evaluate(dr,1,[k])[1] for k in kvec]
ivec = [Dolo.evaluate(dr,1,[k])[2] for k in kvec]
nvec_d = [Dolo.evaluate(drd,1,[k])[1] for k in kvec]
ivec_d = [Dolo.evaluate(drd,1,[k])[2] for k in kvec]
nvec_0 = [Dolo.evaluate(dr0,1,[k])[1] for k in kvec]
ivec_0 = [Dolo.evaluate(dr0,1,[k])[2] for k in kvec]

@assert maximum(abs(nvec_d-nvec))<1e-5
@assert maximum(abs(nvec_0-nvec))<1e-5 # not satisfied right now (see tol. of optimizer)




# does not work yet
filename = joinpath(path,"examples","models","rbc_dtcc_ar1.yaml")
model = Dolo.yaml_import(filename)

Dolo.discretize(model.exogenous)

@time dr = Dolo.time_iteration(model)
@time drv = Dolo.evaluate_policy(model, dr, verbose=true)
