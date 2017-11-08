import Dolo

import Dolo: n_nodes, n_inodes, nodes, CachedDecisionRule
import Dolo: invert_jac

model = Dolo.yaml_import("examples/models/rbc_dtcc_iid.yaml")
# model = Dolo.yaml_import("examples/models/sudden_stop.yaml")

@time dr = Dolo.improved_time_iteration(model, verbose=true, method=:iti, complementarities=true)

@time dr = Dolo.improved_time_iteration(model, verbose=false, method=:gmres)

@time dr = Dolo.improved_time_iteration(model, verbose=false, method=:iti)
