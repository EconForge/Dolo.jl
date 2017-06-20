
import Dolo
import PyPlot
path = Dolo.pkg_path

filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)

n_steps=100
s0 = model.calibration[:states]
m0 = model.calibration[:exogenous]
index = findfirst(model.symbols[:states],:z)
bounds = [dr.grid_endo.min[index], dr.grid_endo.max[index]]


df = Dolo.tabulate(model, dr, :z, bounds, s0, m0)
df = Dolo.tabulate(model, dr, :z, s0, m0)

cn =[:i, :n]
cn[1]
Svalues = linspace(bounds[1], bounds[2], n_steps)
PyPlot.plot(Svalues, df[cn[2]])

for j in cn
  fig = PyPlot.figure(j,figsize=(3,3))
  fig
  PyPlot.plot(Svalues, df[j], label=j)
  PyPlot.legend()
  # PyPlot.xlabel('state = {} | mstate = {}'.format(state, i0))
end

plot_dr = Dolo.tabulate(model, dr, :k, bounds, s0, m0, cn)
