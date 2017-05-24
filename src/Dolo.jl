# __precompile__()

module Dolo

using MacroTools
using DataStructures: OrderedDict
using YAML: load_file, load
using Requests: get
using NLsolve
using Optim
using QuantEcon
using Dolang
using Dolang: _to_expr
# using Distributions: MvNormal, Distribution
import Distributions
import ForwardDiff
import YAML
using Compat; import Compat: String, view

       # model functions
export arbitrage, transition, auxiliary, value, expectation,
       direct_response, controls_lb, controls_ub, arbitrage_2, felicity,

       # mutating version of model functions
       arbitrage!, transition!, auxiliary!, value!, expectation!,
       direct_response, controls_lb!, controls_ub!, arbitrage_2!, felicity!

       # dolo functions
export yaml_import, eval_with, evaluate, evaluate!, model_type, name, filename, id

export time_iteration, value_iteration, steady_state_residuals, simulation

export ModelCalibration, FlatCalibration, GroupedCalibration

# set up core types
abstract AbstractSymbolicModel{ID}
# abstract AbstractSymbolicModel{ID} <: AbstractModel{ID}
abstract AbstractModel{ID} <: AbstractSymbolicModel{ID}

typealias ASModel AbstractSymbolicModel
typealias AModel AbstractModel

id{ID}(::AbstractModel{ID}) = ID

#
# # duplicate hierarchy to experiment with
# @compat abstract type ASModel{ID} end                # symbolic model
# @compat abstract type AModel{ID} <: ASModel{ID} end      # numeric model
#



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
include("numeric/simulations.jl")
include("numeric/derivatives.jl")

include("util.jl")
include("symbolic.jl")
include("calibration.jl")
include("minilang.jl")
include("numeric.jl")
include("model.jl")
include("printing.jl")

include("algos/steady_state.jl")
include("algos/time_iteration.jl")
include("algos/time_iteration_direct.jl")
include("algos/value_iteration.jl")
include("algos/perturbation.jl")
include("algos/simulation.jl")


end # module
