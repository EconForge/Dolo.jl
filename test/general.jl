using Dolo
const NL = Dolo

root_dir = pkgdir(Dolo)
# model = yaml_import("$(root_dir)/examples/ymodels/rbc_mc.yaml")
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")

model.calibration

model2 = NL.recalibrate(model; β=0.8)

dmodel = NL.discretize(model)

sol = NL.time_iteration(dmodel; verbose=false)

tab = tabulate(model, sol.dr, :k;)

sim = NL.simulate(model, sol.dr)



model = include("$(root_dir)/examples/ymodels/rbc_iid.jl")

s0 = NL.draw(model.states)

x = SVector( 0.2, 0.3)

gen = NL.τ(model, s0, x)

NL.calibrated(NL.QP, model, :states)
# NL.calibrated(NL.QP, model)

NL.QP(model.states, [0.2, 0.4])

    

NL.QP(model.states, SVector(0.0, 5.5))


sol = Dolo.time_iteration(model;  verbose=false);
Dolo.time_iteration(model;  verbose=false, improve=true, improve_wait=0, improve_K=500);


# here the problem is the DFun is not initialized with the right
# variables

tab = tabulate(model, sol.dr, :k)

s0 = NL.draw(model.states)

sim = NL.simulate(model, sol.dr, s0)




model2 = include("$(root_dir)/examples/ymodels/consumption_savings.jl")

@time Dolo.time_iteration(model; verbose=false);
@time Dolo.time_iteration(model; verbose=false, improve=true, improve_K=1000);

# @time Dolo.vfi(model; improve=true);





model = include("$(root_dir)/examples/ymodels/rbc_ar1.jl")

Dolo.time_iteration(model; verbose=false);
Dolo.time_iteration(model; verbose=false, improve=true, improve_K=1000);

Dolo.time_iteration(model; verbose=false, improve=true, improve_K=1000, interp_method=:cubic);

# @time Dolo.vfi(model; improve=true);
