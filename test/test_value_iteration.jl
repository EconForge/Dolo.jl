import Dolo

model = Dolo.yaml_import("rbc_dtcc_iid.yaml")

# define
include("tmp_module.jl")
import temp

# get grid for endogenous
gg = model.options.grid
gg.a
gg.orders             # this is imported from the yaml file
grid = CartesianGrid(gg.a, gg.b, gg.orders) # temporary compatibility


# now look at exogenous
# everything below should work for any type of process (normal, markov chain, ar1, ...)

ptype = "MC"
if ptype == "IID"
    process = MvNormal(0.01)
elseif ptype == "MC"
    values = [ -0.01  0.01 ]'
    transitions = [
        0.9 0.1;
        0.1 0.9
    ]
    process = DiscreteMarkovProcess(transitions, values)
end

dprocess = discretize(process)

model.calibration.flat[:beta]

model.calibration[:parameters]
#beta=model.calibration.flat[:beta]
β = 0.5

drv_answer = value_iteration(model, dprocess, process, grid, β, true, 300)

println(drv_answer)
using Gadfly
kvec = linspace(8,12,100)

drv_answer.grid
yvec = [drv_answer(1,[k])[1] for k in kvec]
plot(x=kvec, y=yvec, Geom.line)
