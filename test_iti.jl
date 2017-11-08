import Dolo

import Dolo: n_nodes, n_inodes, nodes, CachedDecisionRule
import Dolo: invert_jac

using ProfileView

model = Dolo.yaml_import("examples/models/rbc_dtcc_mc.yaml")
# model = Dolo.yaml_import("examples/models/sudden_stop.yaml")

@time dr = Dolo.improved_time_iteration(model, verbose=true, method=:gmres, complementarities=true)

function time_it(t=100)
    l = []
    for i=1:t
        sol = Dolo.improved_time_iteration(model, verbose=false, method=:gmres, complementarities=true)
        push!(l, sol)
    end
    return true
end
Profile.clear()
@profile time_it(10)
ProfileView.view()

@time dr = Dolo.improved_time_iteration(model, verbose=false, method=:gmres)

@time dr = Dolo.improved_time_iteration(model, verbose=false, method=:iti)
