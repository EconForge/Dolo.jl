# import Dolo
#
# url = "examples/models/rbc_dtcc_mc.yaml"
# typeof(url)
#
s

# model = Dolo.SModel(url)
# data = model.data
#
# @time Dolo.get_name(model)
# @time Dolo.get_symbols(model)
# @time Dolo.get_definitions(model)
# @time Dolo.get_equations(model)
# @time Dolo.get_infos(model)
# @time Dolo.get_options(model)
# @time Dolo.get_definitions(model)
#
# # apparently the time is dominated by the solution of the calibration
# # the result should be cached
#
# @time calib = Dolo.get_calibration(model)
# @time domain = Dolo.get_domain(model)
# @time exo = Dolo.get_exogenous(model)
# @time grid = Dolo.get_grid(model)
#
# Dolo.set_calibration(model, :k, 0.1)
# Dolo.get_domain(model)
#

csubs(:(a + b + b(1)), :b, :(c+d(1)) == :(a+c+d(1)+c(1)+d(2))


import Dolo

url = "examples/models/rbc_dtcc_mc.yaml"
typeof(url)

nmodel = Dolo.Model(url; print_code=false)

# nmodel.name
nmodel.calibration[:states, :controls]
# nmodel.grid
# nmodel.exogenous
nmodel.equations[:arbitrage]
m,s,x,p = nmodel.calibration[:exogenous,:states,:controls,:parameters]
Dolo.transition(nmodel, m,s,x,m,p)

Dolo.arbitrage(nmodel, m,s,x,m,s,x,p)


dr = Dolo.time_iteration(nmodel)
sol = Dolo.solve_policy(nmodel, dr)


nmodel.grid
