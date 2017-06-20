# __precompile__(false)

module Dolo

# model import utils
using DataStructures: OrderedDict
import YAML; using YAML: load_file, load
using Requests: get

# solvers
using NLsolve
using Optim

# QuantEcon
using QuantEcon; import QuantEcon: simulate
const QE = QuantEcon

# Dolang
using Dolang
using Dolang: _to_expr, inf_to_Inf, solution_order, solve_triangular_system

# Numerical Tools
using MacroTools  # used for eval_with
import Distributions

# Simulation/presentation
using AxisArrays

# Compat across julia versions
using Compat; import Compat: String, view


# exports
       # model functions
export arbitrage, transition, auxiliary, value, expectation,
       direct_response, controls_lb, controls_ub, arbitrage_2, felicity,

       # mutating version of model functions
       arbitrage!, transition!, auxiliary!, value!, expectation!,
       direct_response, controls_lb!, controls_ub!, arbitrage_2!, felicity!

       # dolo functions
export yaml_import, eval_with, evaluate, evaluate!, model_type, name, filename, id

export time_iteration, improved_time_iteration, value_iteration, residuals, simulation

export ModelCalibration, FlatCalibration, GroupedCalibration

# set up core types
@compat abstract type AbstractSymbolicModel{ID} end
@compat abstract type AbstractModel{ID} <: AbstractSymbolicModel{ID} end

@compat const ASModel = AbstractSymbolicModel
@compat const AModel = AbstractModel

id{ID}(::AbstractModel{ID}) = ID

# recursively make all keys at any layer of nesting a symbol
# included here instead of util.jl so we can call it on RECIPES below
_symbol_dict(x) = x
_symbol_dict(d::Associative) =
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

include("numeric/splines/splines.jl")
import .splines

include("numeric/newton.jl")
include("numeric/grids.jl")
include("numeric/processes.jl")
include("numeric/decision_rules.jl")
include("numeric/derivatives.jl")

include("util.jl")
include("symbolic.jl")
include("calibration.jl")
include("minilang.jl")
include("model.jl")
include("printing.jl")

include("algos/steady_state.jl")
include("algos/time_iteration.jl")
include("algos/improved_time_iteration.jl")
include("algos/time_iteration_direct.jl")
include("algos/value_iteration.jl")
include("algos/perturbation.jl")
include("algos/simulation.jl")


end # module
