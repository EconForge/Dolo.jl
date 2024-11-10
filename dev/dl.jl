using Dolo

model = include("../examples/ymodels/consumption_savings_iid.jl")


# dmodel = Dolo.discretize(model)

# wk = Dolo.time_iteration_workspace(dmodel)

# r = Dolo.F(dmodel, wk.x0, wk.φ)

import Dolo: F

function F(model::Dolo.AModel,s,φ)



end

Dolo.rand(model.states)