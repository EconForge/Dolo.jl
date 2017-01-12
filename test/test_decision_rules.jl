include("tmp_module.jl")

include("test_processes.jl")
# this loads:
# - mc: markov chain
# - mvn: multivariate normal law


grid = CartesianGrid([0.0,1.0], [0.1, 0.4], [5,4])
gg = nodes(grid)
n_x = 3
values = rand(size(gg,1), n_x)
vv = [values for i=1:size(mc.values,1)]
dr = DecisionRule(mc, grid, vv)
set_values(dr,vv) # filter coefficients

dr(1, [0.1 0.5])
s0 = [0.2, 0.5]'
res =  dr(1, s0)


grid = CartesianGrid([0.0,1.0], [0.1, 0.4], [5,4])
gg = nodes(grid)
n_x = 3
values = rand(size(gg,1), n_x)
vv = [values]
dr = DecisionRule(mvn, grid, vv)
set_values(dr, values) # filter coefficients
s0 = [0.2, 0.5]'
res =  dr(1, s0)
dr(1, [0.1 0.5])
