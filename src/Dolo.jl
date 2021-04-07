__precompile__(false)

module Dolo

using Printf

import Dolang: SymExpr, list_syms

# model import utils
using DataStructures: OrderedDict
import YAML; using YAML: load_file, load

import IterTools
import HTTP
using StaticArrays

# solvers
using NLsolve
using Optim

# QuantEcon
using QuantEcon; import QuantEcon: simulate
const QE = QuantEcon

# Dolang
using Dolang
using Dolang: _to_expr, inf_to_Inf, solution_order, solve_triangular_system, _get_oorders
import Dolang: Language, add_language_elements!, FromGreek

# Numerical Tools
using MacroTools  # used for eval_with
import Distributions
import BasisMatrices
const BM = BasisMatrices

# Simulation/presentation
using AxisArrays
using StringDistances

# Compat across julia versions
using Compat; import Compat: String, view, @__dot__

using StaticArrays
using IterativeSolvers

# Functions from base we extend
# import Base.A_mul_B!
import Base.size
import Base.eltype
import Base.*

using LinearAlgebra


# exports
       # model functions
export arbitrage, transition, auxiliary, value, expectation,
       direct_response, controls_lb, controls_ub, arbitrage_2, felicity,

       # mutating version of model functions
       arbitrage!, transition!, auxiliary!, value!, expectation!,
       direct_response, controls_lb!, controls_ub!, arbitrage_2!, felicity!,
       evaluate_definitions

       # dolo functions
export lint, yaml_import, eval_with, evaluate, evaluate!, model_type, name, filename, id, features, set_calibration!

export time_iteration, improved_time_iteration, value_iteration, residuals,
        response, simulate, perfect_foresight, time_iteration_direct, find_deterministic_equilibrium, perturb, tabulate

export ModelCalibration, FlatCalibration, GroupedCalibration
export AbstractModel, AbstractDecisionRule, Model

# set up core typesr
# abstract type AbstractSymbolicModel{ID} end
# abstract type AbstractModel{ID} <: AbstractSymbolicModel{ID} end

abstract type AbstractModel{ID} end

const AModel = AbstractModel

id(::AbstractModel{ID}) where {ID} = ID

Expression = Union{Expr, Symbol, Float64, Int64}
 
# conventions for list of points
Point{d} = SVector{d,Float64}
Value{n} = SVector{n,Float64}
ListOfPoints{d} = Vector{Point{d}}
ListOfValues{n} = Vector{Value{n}}


# recursively make all keys at any layer of nesting a symbol
# included here instead of util.jl so we can call it on RECIPES below
_symbol_dict(x) = x
_symbol_dict(d::AbstractDict) =
    Dict{Symbol,Any}([(Symbol(k), _symbol_dict(v)) for (k, v) in d])

const src_path = dirname(@__FILE__)
const pkg_path = dirname(src_path)
const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))

# define these _functions_ so users can add their own _methods_ to extend them
for f in [:arbitrage, :transition, :auxiliary, :value, :expectation,
          :direct_response, :controls_lb, :controls_ub, :arbitrage_2,
          :arbitrage!, :transition!, :auxiliary!, :value!, :expectation!,
          :direct_response, :controls_lb!, :controls_ub!, :arbitrage_2!]
    Core.eval(Dolo, Expr(:function, f))
end

include("numeric/splines/splines.jl")
import .splines
import .splines: eval_UC_spline, eval_UC_spline!, prefilter!

include("util.jl")

include("numeric/complementarities.jl")
include("numeric/newton.jl")
include("numeric/grids.jl")
include("numeric/processes.jl")

minilang = Language(Dict())
add_language_elements!(minilang, Dict(
    "!Normal"=>Normal,
    "!MarkovChain"=>MarkovChain,
    "!Product"=>Product,
    "!PoissonProcess"=>PoissonProcess,
    "!DeathProcess"=>DeathProcess,
    "!AgingProcess"=>AgingProcess,
    "!VAR1"=>VAR1,
))

include("linter.jl")
include("calibration.jl")
include("minilang.jl")
include("new_model.jl")
# include("model.jl")
# include("symbolic.jl")
# include("printing.jl")

include("numeric/decision_rules/core.jl")
include("numeric/decision_rules/csplines.jl")
include("numeric/decision_rules/constructor.jl")
include("numeric/decision_rules/compat.jl")
include("numeric/decision_rules/smolyak.jl")
include("numeric/decision_rules/complete.jl")

include("algos/steady_state.jl")
include("algos/time_iteration.jl")
include("algos/improved_time_iteration.jl")
include("algos/time_iteration_direct.jl")
include("algos/value_iteration.jl")
include("algos/perturbation.jl")
include("algos/simulation.jl")
include("algos/perfect_foresight.jl")

include("algos/ergodic.jl")

end # module
