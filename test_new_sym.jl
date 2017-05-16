
import Dolo


url = "examples/models/rbc_dtcc_mc.yaml"
typeof(url)

model = Dolo.SModel(url)
data = model.data

@time Dolo.get_name(model)
@time Dolo.get_symbols(model)
@time Dolo.get_definitions(model)
@time Dolo.get_equations(model)
@time Dolo.get_infos(model)
@time Dolo.get_options(model)
@time Dolo.get_definitions(model)

# apparently the time is dominated by the solution of the calibration
# the result should be cached

@time calib = Dolo.get_calibration(model)
@time domain = Dolo.get_domain(model)
@time exo = Dolo.get_exogenous(model)
@time grid = Dolo.get_grid(model)


Dolo.set_calibration(model, :k, 0.1)
Dolo.get_domain(model)

nmodel = Dolo.Model(url)
nmodel.name
nmodel.calibration[:states, :controls]
nmodel.grid
nmodel.exogenous

nmodel.equations[:controls_ub]
keys(nmodel.equations)
nmodel.equations[:arbitrage]
nmodel.equations[:controls_lb]

import Dolang
import DataStructures

# assume dolo style model is default, special case others
function _get_args(sm::Dolo.Model, spec)
    # get args
    args = DataStructures.OrderedDict{Symbol,Vector{Tuple{Symbol,Int}}}()
    for (grp, shift, nm) in spec[:eqs]
        grp == "parameters" && continue

        syms = sm.symbols[Symbol(grp)]
        args[Symbol(nm)] = Tuple{Symbol,Int}[(v, shift) for v in syms]
    end
    args
end

type StupidType end

function Dolang.FunctionFactory(sm::Dolo.Model, func_nm::Symbol)

    spec = Dolo.RECIPES[:dtcc][:specs][func_nm]
    eqs = sm.equations[func_nm]

    # get targets
    target = get(spec, :target, [nothing])[1]
    has_targets = !(target === nothing)
    targets = has_targets ? sm.symbols[Symbol(target)] : Symbol[]

    # get other stuff
    args = _get_args(sm, spec)
    params = sm.symbols[:parameters]
    dispatch = StupidType
    # dispatch = Dolo._numeric_mod_type(sm)

    definitions = Dolo.get_definitions(model)

    Dolang.FunctionFactory(dispatch, eqs, args, params, targets=targets,
                    defs=definitions, funname=func_nm)

end
