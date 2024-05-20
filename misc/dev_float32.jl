using Dolo

root_dir = pkgdir(Dolo)
model = include("$(root_dir)/examples/ymodels/rbc_mc.jl")

dm = Dolo.discretize(model, Dict(:endo=>[10000000]) )

Dolo.convert