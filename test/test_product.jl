using Gadfly
import Dolo
using AxisArrays

path = Dolo.pkg_path


model = Dolo.yaml_import(
        joinpath(path,
         "examples","models",
         "rbc_catastrophe.yaml")
)

Dolo.discretize(Dolo.MarkovChain, model.exogenous)

dp = Dolo.discretize(model.exogenous)

sol = Dolo.time_iteration(model, dp)

sim = Dolo.simulate(model, sol.dr, dp, T=100)
sim
plot(y=sim[Axis{:N}(1), Axis{:V}(:k)],Geom.line)

plot(y=sim[Axis{:N}(1), Axis{:V}(:mc_process)],Geom.line)

dp












Gadfly.plot(y=)
