path = Pkg.dir("Dolo")

import Dolo

fn = joinpath(path,"examples","models","rbc_dtcc_mc.yaml")
model_mc = Dolo.yaml_import(fn)

drc = Dolo.ConstantDecisionRule(model_mc.calibration[:controls])
@time dr0, drv0 = Dolo.solve_policy(model_mc, drc; n_eval=1000) #, verbose=true, maxit=10000 )




# two loops:
# - VFI:
    # stops when err_v M tol_v

# optimize then inner_loop
# - inner: update until:
        - it_inner > maxit_inner
        - err_h < tol_h
