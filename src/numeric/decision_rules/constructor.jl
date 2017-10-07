# cubic is default

abstract type DRSelector end
struct Linear <: DRSelector end
struct Cubic <: DRSelector end

struct CompletePolnomial{order} <: DRSelector end
struct Smolyak <: DRSelector end
struct Chebyshev <: DRSelector end
struct BSpline{order} <: DRSelector end
struct PWLinear <: DRSelector end

# incomplete type inference:
DecisionRule(exo_grid, endo_grid, i::Int64) = DecisionRule(exo_grid, endo_grid, Val{i})
DecisionRule(exo_grid, endo_grid, i::Int64, interp_type) = DecisionRule(exo_grid, endo_grid, Val{i}, interp_type)


# cubic interpolation
function DecisionRule(exo_grid::EmptyGrid, endo_grid::CartesianGrid, ::Type{Val{nx}}, ::Type{Cubic}) where nx
    CubicDR(exo_grid, endo_grid, Val{nx})
end

function DecisionRule(exo_grid::UnstructuredGrid, endo_grid::CartesianGrid, ::Type{Val{nx}}, ::Type{Cubic}) where nx
    CubicDR(exo_grid, endo_grid, Val{nx})
end

function DecisionRule(exo_grid::CartesianGrid, endo_grid::CartesianGrid, ::Type{Val{nx}}, ::Type{Cubic}) where nx
    CubicDR(exo_grid, endo_grid, Val{nx})
end

# Smolyak interpolation
function DecisionRule(
        exo_grid::S, endo_grid::SmolyakGrid, ::Type{Val{nx}}, ::Type{Smolyak}
    ) where S <: Union{EmptyGrid,<:UnstructuredGrid} where nx
    SmolyakDR(exo_grid, endo_grid, Val{nx})
end

# Chebyshev interpolation
function DecisionRule(
        exo_grid::S, endo_grid::CartesianGrid, ::Type{Val{nx}}, ::Type{Chebyshev}
    ) where S <: Union{EmptyGrid,<:UnstructuredGrid} where nx
    BasisMatricesDR(exo_grid, endo_grid, Val{nx}, BM.ChebParams)
end

# BSpline interpolation
function DecisionRule(
        exo_grid::S, endo_grid::CartesianGrid, ::Type{Val{nx}}, ::Type{BSpline{order}}
    ) where S <: Union{EmptyGrid,<:UnstructuredGrid} where nx where order
    BasisMatricesDR(exo_grid, endo_grid, Val{nx}, BM.SplineParams, (order,))
end

# Piecewise Linear interpolation
function DecisionRule(
        exo_grid::S, endo_grid::CartesianGrid, ::Type{Val{nx}}, ::Type{PWLinear}
    ) where S <: Union{EmptyGrid,<:UnstructuredGrid} where nx
    BasisMatricesDR(exo_grid, endo_grid, Val{nx}, BM.LinParams)
end

# Complete polynomials
function DecisionRule(
        exo_grid::S, endo_grid::Grid, ::Type{Val{nx}}, ::Type{CompletePolnomial{order}}
    ) where S <: Union{EmptyGrid,<:UnstructuredGrid} where nx where order
    CompletePolyDR(exo_grid, endo_grid, Val{nx}, order)
end

# default to order 3
function DecisionRule(
        exo_grid::S, endo_grid::Grid, ::Type{Val{nx}}, ::Type{<:CompletePolnomial}
    ) where S <: Union{EmptyGrid,<:UnstructuredGrid} where nx
    CompletePolyDR(exo_grid, endo_grid, Val{nx})
end

# generic constructors

function DecisionRule(exo_grid::Grid, endo_grid::Grid, vals::Vector{ListOfPoints{n_x}}) where n_x
    dr = DecisionRule(exo_grid, endo_grid, Val{n_x}, Cubic)
    set_values!(dr, vals)
    dr
end

function DecisionRule(exo_grid::Grid, endo_grid::Grid, vals::Vector{ListOfPoints{n_x}}, tt) where n_x
    dr = DecisionRule(exo_grid, endo_grid, Val{n_x}, tt)
    set_values!(dr, vals)
    dr
end

# #####
# ##### Cached Decision Rules (do we really need them ?)
# #####

type CachedDecisionRule{T,S}
    dr::T
    process::S
end

function CachedDecisionRule(
        process::AbstractDiscretizedProcess, grid::Grid, n_x::Int, ::Type{T}=Cubic
    ) where T <: DRSelector
    CachedDecisionRule(DecisionRule(process.grid, grid, n_x, T), process)
end

function CachedDecisionRule(
        process::AbstractDiscretizedProcess, grid::Grid, values, ::Type{T}=Cubic
    ) where T <: DRSelector
    CachedDecisionRule(DecisionRule(process.grid, grid, values, T), process)
end

set_values!(cdr::CachedDecisionRule, v) = set_values!(cdr.dr, v)
nodes(cdr::CachedDecisionRule) = nodes(cdr.dr)

# defaults

(cdr::CachedDecisionRule)(v::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,v)
(cdr::CachedDecisionRule)(m::Union{AbstractVector,AbstractMatrix}, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr, m, s)
(cdr::CachedDecisionRule)(i::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,node(cdr.process, i), s)
(cdr::CachedDecisionRule)(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,inode(cdr.process, i, j), s)

(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:CartesianGrid}, <:DiscreteMarkovProcess})(i::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,i, s)
(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:CartesianGrid}, <:DiscreteMarkovProcess})(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,j, s)
