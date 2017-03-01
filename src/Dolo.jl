# __precompile__()

module Dolo

using MacroTools
using DataStructures: OrderedDict
using YAML: load_file, load
using Requests: get
using NLsolve
using Optim
using QuantEcon
using Distributions: MvNormal, Distribution
using Dolang
using Dolang: _to_expr
import ForwardDiff
import YAML
using Compat; import Compat: String, view

# export AbstractModel, AbstractSymbolicModel, AbstractNumericModel, ASM, ANM,
#        SymbolicModel, NumericModel,
#        FlatCalibration, GroupedCalibration, ModelCalibration, TaylorExpansion,
#        RECIPES,

       # model functions
export arbitrage, transition, auxiliary, value, expectation,
       direct_response, controls_lb, controls_ub, arbitrage_2, felicity,

       # mutating version of model functions
       arbitrage!, transition!, auxiliary!, value!, expectation!,
       direct_response, controls_lb!, controls_ub!, arbitrage_2!, felicity!

       # dolo functions
export yaml_import, eval_with, evaluate, evaluate!, model_type, name, filename, id

export time_iteration, value_iteration, steady_state_residuals

# set up core types
abstract AbstractModel{ID}
abstract AbstractSymbolicModel{ID} <: AbstractModel{ID}
abstract AbstractNumericModel{ID} <: AbstractModel{ID}

typealias ASM AbstractSymbolicModel
typealias ANM AbstractNumericModel


id{ID}(::AbstractModel{ID}) = ID

# recursively make all keys at any layer of nesting a symbol
# included here instead of util.jl so we can call it on RECIPES below
_symbol_dict(x) = x
@compat _symbol_dict(d::Associative) =
    Dict{Symbol,Any}([(Symbol(k), _symbol_dict(v)) for (k, v) in d])

const src_path = dirname(@__FILE__)
const pkg_path = dirname(src_path)
const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))

# define these _functions_ so users can add their own _methods_ to extend them
for f in [:arbitrage, :transition, :auxiliary, :value, :expectation,
          :direct_response, :controls_lb, :controls_ub, :arbitrage_2,
          :arbitrage!, :transition!, :auxiliary!, :value!, :expectation!,
          :direct_response, :controls_lb!, :controls_ub!, :arbitrage_2!]
    eval(Expr(:function, f))
end


include("numeric/newton.jl")
include("numeric/grids.jl")
include("numeric/processes.jl")
include("numeric/decision_rules.jl")
include("numeric/taylor_series.jl")
include("numeric/simulations.jl")
include("numeric/derivatives.jl")

include("util.jl")
include("symbolic.jl")
include("calibration.jl")
include("minilang.jl")
include("numeric.jl")
include("model_import.jl")

include("algos/steady_state.jl")
include("algos/time_iteration.jl")
include("algos/time_iteration_direct.jl")
include("algos/value_iteration.jl")
include("algos/perturbation.jl")


end # module
