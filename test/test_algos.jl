import Dolo
path = Pkg.dir("Dolo")



filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
model_mc = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model_mc, verbose=false)
# @time dr = Dolo.time_iteration_direct(model_mc, verbose=true)
@time drv = Dolo.evaluate_policy(model_mc, dr, verbose=true)
@time drd = Dolo.time_iteration_direct(model_mc, dr, verbose=true, maxit=500)

# compare with prerecorded values
kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
nvec = [Dolo.evaluate(dr,1,[k])[1] for k in kvec]
ivec = [Dolo.evaluate(dr,1,[k])[2] for k in kvec]
ivec_test = [0.295977,  0.257538,  0.21566,  0.173564,  0.132103,  0.0915598,  0.0520067,  0.0134661,  7.01983e-6, 3.40994e-17]
nvec_test = [ 0.391997,  0.348033,  0.318369,  0.296276,  0.278821,  0.264487,  0.25239 ,  0.241974,  0.236604,  0.233779 ]
@assert maximum(abs(ivec-ivec_test))<1e-5
@assert maximum(abs(nvec-nvec_test))<1e-5

# compare time_iteration and time_iteration_direct
nvec_d = [Dolo.evaluate(drd,1,[k])[1] for k in kvec]
ivec_d = [Dolo.evaluate(drd,1,[k])[2] for k in kvec]
@assert maximum(abs(nvec_d-nvec))<1e-5



# this one needs a lower value of beta or a better initial guess
filename = joinpath(path,"examples","models","rbc_dtcc_iid.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration_direct(model, verbose=true)
@time dr = Dolo.time_iteration(model, verbose=true)
@time drv = Dolo.evaluate_policy(model, dr, verbose=true)

@time dr = Dolo.time_iteration_direct(model, dr, verbose=true)

kvec = linspace(dr.grid.min[1],dr.grid.max[1],10)
nvec = [Dolo.evaluate(dr,1,[k])[1] for k in kvec]
ivec = [Dolo.evaluate(dr,1,[k])[2] for k in kvec]



kvec = linspace(8,12)
nvec = [Dolo.evaluate(dr,1,[u])[1] for u in kvec]
nvec


# does not work yet
filename = joinpath(path,"examples","models","rbc_dtcc_ar1.yaml")
model = Dolo.yaml_import(filename)

Dolo.discretize(model.exogenous)

@time dr = Dolo.time_iteration(model)
@time drv = Dolo.evaluate_policy(model, dr, verbose=true)
