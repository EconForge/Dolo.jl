import Dolo
include("tmp_module.jl")
import temp
include("steady_state.jl")
dolo_dir = Dolo.pkg_path

model = Dolo.yaml_import("rbc_dtcc_iid.yaml")


process = nothing
ptype = "MC"
if ptype == "IID"
    process = MvNormal(0.01)
elseif ptype == "MC"

    values = [ 0.0 ]'
    transitions = [
        1.0
    ]'

    process = DiscreteMarkovProcess(transitions, values)
end

dp = discretize(process)


steady_state(model, model.calibration)

@time dr = time_iteration(model, process, verbose=true)

ivals = [dr(1, [i])[1] for i=kvec ]
nvals = [dr(1, [i])[2] for i=kvec ]


using Gadfly

plot(x=kvec, y=ivals)

plot(x=kvec, y=nvals)

kkvec = [dr(1,[k])[1] for k in kvec]
plot(x=kvec, y=yvec, Geom.line)
