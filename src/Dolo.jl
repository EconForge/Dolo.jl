# __precompile__()

module Dolo

using MacroTools
using DataStructures: OrderedDict
using YAML: load_file, load
using Requests: get
using NLsolve
using QuantEcon
using Distributions: MvNormal, Distribution
using Dolang
using Dolang: _to_expr
import ForwardDiff
import YAML
using Compat; import Compat: String, view

export AbstractModel, AbstractSymbolicModel, AbstractNumericModel, ASM, ANM,
       SymbolicModel, NumericModel,
       FlatCalibration, GroupedCalibration, ModelCalibration, TaylorExpansion,
       RECIPES,

       # model functions
       arbitrage, transition, auxiliary, value, expectation,
       direct_response, controls_lb, controls_ub, arbitrage_2, felicity,

       # mutating version of model functions
       arbitrage!, transition!, auxiliary!, value!, expectation!,
       direct_response, controls_lb!, controls_ub!, arbitrage_2!, felicity!,

       # dolo functions
       yaml_import, eval_with, evaluate, evaluate!, model_type, name, filename,
       linear_solve, simulate, id

# set up core types
abstract AbstractModel{ID,kind}
abstract AbstractSymbolicModel{ID,kind} <: AbstractModel{ID,kind}
abstract AbstractNumericModel{ID,kind} <: AbstractModel{ID,kind}

typealias AbstractDTCSCC{ID} AbstractNumericModel{ID,:dtcscc}
typealias AbstractDTMSCC{ID} AbstractNumericModel{ID,:dtmscc}

typealias ASM AbstractSymbolicModel
typealias ANM AbstractNumericModel

abstract AbstractDecisionRule

id{ID}(::AbstractModel{ID}) = ID
model_type{_,kind}(::AbstractModel{_,kind}) = kind

# recursively make all keys at any layer of nesting a symbol
# included here instead of util.jl so we can call it on RECIPES below
_symbol_dict(x) = x
@compat _symbol_dict(d::Associative) =
    Dict{Symbol,Any}(Symbol(k) => _symbol_dict(v) for (k, v) in d)

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


include("util.jl")
include("symbolic.jl")
include("calibration.jl")
include("minilang.jl")
include("numeric.jl")
include("model_import.jl")

include("numeric/taylor_series.jl")
include("numeric/simulations.jl")
include("numeric/derivatives.jl")

include("algos/dtcscc.jl")

end # module
