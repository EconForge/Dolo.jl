module Dolo

using MacroTools
using DataStructures: OrderedDict
using YAML: load_file, load
using Requests: get
using NLsolve
using QuantEcon
using Distributions: MvNormal
import ForwardDiff

export AbstractModel, AbstractSymbolicModel, AbstractNumericModel, ASM, ANM,
       SymbolicModel, DTCSCCModel, DTMSCCModel,
       FlatCalibration, GroupedCalibration, ModelCalibration, TaylorExpansion,
       RECIPES,

       # model functions
       arbitrage, transition, auxiliary, value, expectation, direct_response,
       controls_lb, controls_ub, arbitrage_2,

       # mutating version of model functions
       arbitrage!, transition!, auxiliary!, value!, expectation!,
       direct_response, controls_lb!, controls_ub!, arbitrage_2!,

       # dolo functions
       yaml_import, eval_with, evaluate, evaluate!, model_type, name, filename,
       linear_solve, simulate, id

# set up core types
abstract AbstractModel{ID,kind}
abstract AbstractSymbolicModel{ID,kind} <: AbstractModel{ID,kind}
abstract AbstractNumericModel{ID,kind} <: AbstractModel{ID,kind}

id{ID}(::AbstractModel{ID}) = ID
model_type{_,kind}(::AbstractModel{_,kind}) = kind

typealias ASM AbstractSymbolicModel
typealias ANM AbstractNumericModel

abstract AbstractDecisionRule

_symbol_dict(x) = x
_symbol_dict(d::Associative) =
    Dict{Symbol,Any}([(symbol(k), _symbol_dict(v)) for (k, v) in d])

const src_path = dirname(@__FILE__)
const pkg_path = dirname(src_path)
const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))

include("util.jl")
include("parser.jl")
include("symbolic.jl")
include("calibration.jl")
include("numeric.jl")
include("model_import.jl")

include("numeric/taylor_series.jl")
include("numeric/simulations.jl")

include("algos/dtcscc.jl")

end # module
