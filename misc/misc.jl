using Dolo

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")



@time Dolo.time_iteration(model; verbose=true, engine=:cpu);
@time Dolo.time_iteration(model; verbose=true, improve=true, improve_K=1000, improve_wait=0);
@time Dolo.time_iteration(model; verbose=true, improve=true, improve_K=2, interpolation=:cubic);




using Dolo

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_iid.jl")


@time Dolo.time_iteration(model; verbose=false);
@time Dolo.time_iteration(model; verbose=false, engine=:cpu);
# @time Dolo.time_iteration(model; verbose=false, engine=:gpu);


@time Dolo.time_iteration(model; verbose=true, improve=true, improve_K=1000, improve_wait=0);
@time Dolo.time_iteration(model; verbose=true, improve=true, improve_K=2, interpolation=:cubic);

####
####
####

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_ar1.jl")


dmodel = Dolo.discretize(model)
wksp = Dolo.time_iteration_workspace(dmodel; interp_mode=:linear);
res = Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000, engine=:default);
res

wksp = Dolo.time_iteration_workspace(dmodel; interp_mode=:linear);
res = Dolo.time_iteration(dmodel, wksp; verbose=true, improve=false, improve_K=1000, engine=:default, trace=true);
res


@time Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000, engine=:default);

@time Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000);


@time Dolo.time_iteration(dmodel; verbose=true, improve=true, improve_wait=0, improve_K=500);




wksp = Dolo.time_iteration_workspace(dmodel; interp_mode=:cubic);
Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000);
@time Dolo.time_iteration(dmodel, wksp; verbose=false, improve=false, improve_K=1000);


@time Dolo.time_iteration(dmodel; verbose=true, improve=true, improve_wait=0, improve_K=500);
@time Dolo.time_iteration(dmodel; verbose=true, improve=true, improve_wait=0, improve_K=500);

