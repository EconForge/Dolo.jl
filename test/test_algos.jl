import Dolo
path = Pkg.dir("Dolo")



filename = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
model_mc = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model_mc, verbose=false)
@time drv = Dolo.evaluate_policy(model_mc, dr, verbose=true)


# this one needs a lower value of beta or a better initial guess
filename = joinpath(path,"examples","models","rbc_dtcc_iid.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model, verbose=true)
@time drv = Dolo.evaluate_policy(model, dr, verbose=true)
@time dr = Dolo.time_iteration_direct(model, true)




# does not work yet
filename = joinpath(path,"examples","models","rbc_dtcc_ar1.yaml")
model = Dolo.yaml_import(filename)

Dolo.discretize(model.exogenous)

@time dr = Dolo.time_iteration(model)
@time drv = Dolo.evaluate_policy(model, dr, verbose=true)
