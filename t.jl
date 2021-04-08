import Dolo: Model

import Dolo
import Dolo: time_iteration
import Dolo: EmptyDomain


model = Model("examples/models/rbc.yaml")


m,s,x,p = model.calibration[:exogenous, :states, :controls, :parameters]



Dolo.transition(model, m,s,x,m,p)
Dolo.perturb(model)



import Dolo: Domain

Dolo.time_iteration(model)


function get_infos(model::AModel)
    if "infos" in keys(model.data)
        return model.data["infos"]
    else
        return Dict()
    end
end



# ####################################################
# ###### TODO   conversion to SymExpr incorrect ######
# ####################################################

# model = Model("examples/models/rbc.yaml")

import Dolang
# ff_a = get_factory(model, "transition")

ff, ff_l, ff_u = Dolo.get_factory(model, "arbitrage")

source_l = Dolang.gen_generated_gufun(ff_l)

l = Base.eval(Dolang, source_l)


# m,s,x,p = model.calibration[:exogenous, :states, :controls, :parameters]

# source_g = Dolang.gen_generated_gufun(ff_a)
# source_f = Dolang.gen_generated_gufun(ff_s)

# g = Base.eval(Dolang, source_g)
# f = Base.eval(Dolang, source_f)

# g(m,s,x,m,p)
# f(m,s,x,m,s,x,p)

