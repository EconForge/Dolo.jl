using Dolo

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_iid.jl")

isbits(model)


Dolo.time_iteration(model; verbose=false);
Dolo.time_iteration(model; verbose=false, improve=true, improve_K=1000);
@time Dolo.time_iteration(model; verbose=false, improve=true, improve_K=1000, interpolation=:cubic);



dmodel = Dolo.discretize(model)

wksp = Dolo.time_iteration_workspace(dmodel; interp_mode=:linear);
res = Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000, engine=:default);


wksp = Dolo.time_iteration_workspace(dmodel; interp_mode=:linear);
res = Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000, engine=:default, trace=true);



@time Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000, engine=:default);

@time Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000);


@time Dolo.time_iteration(dmodel; verbose=true, improve=true, improve_wait=0, improve_K=500);

@time Dolo.time_iteration(dmodel; verbose=true, improve=true, improve_wait=0, improve_K=500, interpolation=:cubic);




wksp = Dolo.time_iteration_workspace(dmodel; interp_mode=:cubic);
Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000);
@time Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000);


@time Dolo.time_iteration(dmodel; verbose=true, improve=true, improve_wait=0, improve_K=500);
@time Dolo.time_iteration(dmodel; verbose=true, improve=true, improve_wait=0, improve_K=500);

