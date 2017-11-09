import Dolo

import Dolo: n_nodes, n_inodes, nodes, CachedDecisionRule
import Dolo: invert_jac
using StaticArrays





module InitDR
    import Dolo: AbstractDiscretizedProcess, Point, node, AbstractDecisionRule, Grid, EmptyGrid
    using StaticArrays

    struct CFunDR{S,T,nx} <: AbstractDecisionRule{S,T,nx}
        fun::Function
        dprocess::AbstractDiscretizedProcess
    end
    CFunDR(fun::Function, dprocess::AbstractDiscretizedProcess) = CFunDR{typeof(dprocess.grid), EmptyGrid, 1}(fun, dprocess)
    (cfdr::CFunDR)(i::Int, x::Point{d}) where d = cfdr.fun(node(Point, cfdr.dprocess,i),x)
    (cfdr::CFunDR)(i::Int, x::Vector{Point{d}}) where d = [cfdr(i,e) for e in x]
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










model = Dolo.yaml_import("examples/models/rbc_dtcc_mc.yaml")
dp = Dolo.discretize(model.exogenous)

m_ss = model.calibration[:exogenous]
x_ss = model.calibration[:controls]
s_ss = model.calibration[:states]


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
