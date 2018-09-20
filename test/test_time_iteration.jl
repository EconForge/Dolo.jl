using Dolo

dolo_dir = Dolo.pkg_path

model = Dolo.yaml_import("examples/models/rbc_iid.yaml")

process = Dolo.MvNormal(0.001)

dp = Dolo.discretize(process)




@time sol = time_iteration(model, dp; verbose=true, maxit=5)

@time dd = Dolo.improved_time_iteration(model, dp, sol.dr; verbose=true)




ivals = [dr(1, [i])[1] for i=kvec ]
nvals = [dr(1, [i])[2] for i=kvec ]


using Gadfly

plot(x=kvec, y=ivals)

plot(x=kvec, y=nvals)

kkvec = [dr(1,[k])[1] for k in kvec]
plot(x=kvec, y=yvec, Geom.line)
