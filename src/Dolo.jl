module Dolo

using MacroTools
using DataStructures: OrderedDict
using YAML: load_file, load
using Requests: get
using NLsolve
import ForwardDiff

export AbstractModel, AbstractSymbolicModel, AbstractNumericModel, ASM, ANM,
       AbstractDoloFunctor, SymbolicModel, DTCSCCModel, DTMSCCModel,
       FlatCalibration, GroupedCalibration, ModelCalibration, TaylorExpansion,
       RECIPES,

       # functions
       yaml_import, eval_with, evaluate, evaluate!, model_type, name, filename,
       DecisionRule, MixedDecisionRule, set_values,
       # defined in dtmscc
       mlinspace, constant_guess, evaluate_policy, solve_policy, residuals, time_iteration

# set up core types
abstract AbstractModel
abstract AbstractSymbolicModel <: AbstractModel
abstract AbstractNumericModel <: AbstractModel

typealias ASM AbstractSymbolicModel
typealias ANM AbstractSymbolicModel

abstract AbstractDoloFunctor

_symbol_dict(x) = x
_symbol_dict(d::Associative) =
    Dict{Symbol,Any}([(symbol(k), _symbol_dict(v)) for (k, v) in d])

const src_path = dirname(@__FILE__)
const pkg_path = dirname(src_path)
const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))

include("util.jl")
include("parser.jl")
include("model_types.jl")
include("model_import.jl")

include("numeric/misc.jl")
include("numeric/taylor_series.jl")
include("numeric/serial.jl")
include("numeric/policies.jl")

include("algos/dtcscc.jl")
include("algos/dtmscc.jl")

end # module
