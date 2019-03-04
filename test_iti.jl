import Dolo

import Dolo: n_nodes, n_inodes, nodes, CachedDecisionRule
import Dolo: invert_jac

model = Dolo.yaml_import("examples/models/rbc_mc.yaml")
dp = Dolo.discretize(model.exogenous)

using StaticArrays
import OrderedCollections
import Dolang
fff = Dolang.FlatFunctionFactory(OrderedCollections.OrderedDict{Symbol,Union{Expr, Number, Symbol}}(:_k__0_=>:((1 - _delta_) * _k_m1_ + _i_m1_)), OrderedCollections.OrderedDict(
:m=>[:_z_m1_, :_z2_m1_],:s=>[:_k_m1_],:x=>[:_n_m1_, :_i_m1_],:M=>[:_z__0_, :_z2__0_],
:p=>[:_beta_, :_sigma_, :_eta_, :_chi_, :_delta_, :_alpha_, :_rho_, :_zbar_, :_sig_z_]),
 Symbol[:_k__0_], OrderedCollections.OrderedDict{Symbol,Union{Expr, Number, Symbol}}(), :transition)

fun_code = Dolang.gen_generated_gufun(fff)

s0 = SVector(model.calibration[:states]...)
x0 = SVector(model.calibration[:controls]...)
m0 = SVector(model.calibration[:exogenous]...)
p = SVector(model.calibration[:parameters]...)

fun(m0, s0, x0, m0, p)



@time sol = Dolo.time_iteration(model, verbose=false, complementarities=false)
@time sol = Dolo.time_iteration(model, sol.dr, verbose=false, complementarities=false)

@time sol = Dolo.time_iteration(model, verbose=false, complementarities=true)

@time sol = Dolo.improved_time_iteration(model, verbose=true, complementarities=false, method=:gmres)
@time sol = Dolo.improved_time_iteration(model, verbose=true, complementarities=true, method=:gmres)


@time sol = Dolo.improved_time_iteration(model, verbose=true, complementarities=false, method=:iti)

@time sol = Dolo.improved_time_iteration(model, verbose=true, complementarities=false, method=:gmres)

@time solv = Dolo.evaluate_policy(model, sol.dr, verbose=false)

@time solv = Dolo.value_iteration(model, verbose=true)


using ProfileView


xl = Profile.init()

yl = (1000000, 0.1)

Profile.init(n=yl[1], delay=yl[2])
Profile.clear()
@profile solv = Dolo.value_iteration(model, verbose=true)
ProfileView.view()


Profile.clear()
@profile sol = Dolo.time_iteration(model, verbose=false, complementarities=true)
ProfileView.view()


Profile.clear()
@profile sol = Dolo.improved_time_iteration(model, verbose=false, complementarities=true, method=:gmres)
ProfileView.view()

Profile.clear()
@profile sol = Dolo.time_iteration_direct(model, verbose=false)
ProfileView.view()



@time sol = Dolo.improved_time_iteration(model, verbose=true, complementarities=true, method=:gmres)


sim = Dolo.simulate(model, sol.dr)
dp = Dolo.discretize(model.exogenous)
grid = model.grid

@time Dolo.evaluate_policy(model, dp, grid, sol.dr, verbose=false)

SVector(model.calibration[:parameters]...)


@time sol = Dolo.time_iteration(model,verbose=true,complementarities=true, verbose=true, dampen=0.1)

Dolo.improved_time_iteration(model, sol.dr,verbose=true)
@time Dolo.improved_time_iteration(model,verbose=true,complementarities=false, verbose=false)


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
