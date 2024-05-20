using Dolo

root_dir = pkgdir(Dolo)
model32 = include("$(root_dir)/misc/rbc_float32.jl")


dm32 = Dolo.discretize(model, Dict(:endo=>[10000000]) )

Dolo.convert