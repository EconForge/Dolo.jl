import Dolo

import Dolo: n_nodes, n_inodes, nodes, CachedDecisionRule
import Dolo: invert_jac


model = Dolo.yaml_import("examples/models/rbc_dtcc_mc.yaml")
dp = Dolo.discretize(model.exogenous)

m_ss = model.calibration[:exogenous]
x_ss = model.calibration[:controls]
s_ss = model.calibration[:states]

@time sol = Dolo.time_iteration(model,verbose=true,complementarities=true, verbose=true, dampen=0.1)

Dolo.improved_time_iteration(model, sol.dr,verbose=true)
@time Dolo.improved_time_iteration(model,verbose=true,complementarities=false, verbose=false)

module InitDR
    import Dolo: AbstractDiscretizedProcess, Point, node, AbstractDecisionRule, Grid, EmptyGrid
    using StaticArrays


end

import InitDR

function iidr(m,s)
    k = s[1]
    kss = 9.3549
    n = 0.33-0.00*(k-kss)
    i = 0.234+0.000*(k-kss)
    SVector(n,i)
end


cfdr = InitDR.CFunDR(iidr, dp)

Dolo.improved_time_iteration(model, verbose=true, complementarities=false)

Dolo.improved_time_iteration(model, cfdr, verbose=true, complementarities=false)


Dolo.improved_time_iteration(model, cfdr, verbose=true)





methods(InitDR.CFunDR)

typeof(iidr)
typeof(dp)<:Dolo.AbstractDiscretizedProcess

cfdr(1, SVector(s_ss...))











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
