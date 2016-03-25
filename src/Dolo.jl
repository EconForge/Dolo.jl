module Dolo

using MacroTools
using DataStructures: OrderedDict
using YAML: load_file

export AbstractModel, AbstractSymbolicModel, AbstractNumericModel, ASM, ANM,
       AbstractDoloFunctor, SymbolicModel, DTCSCCModel, DTMSCCModel

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


end # module
