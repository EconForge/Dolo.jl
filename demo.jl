a = 4

using Dolo


model = include("examples/ymodels/consumption_savings.jl")

Dolo.time_iteration(model)


