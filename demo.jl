model = include("examples/ymodels/rbc_mc.jl")

dm = Dolo.discretize(model)

@time Dolo.time_iteration(dm, verbose=falsefalse);

@time Dolo.time_iteration(dm; verbose=true, improve=true);

@time wk = Dolo.time_iteration_workspace(dm);

