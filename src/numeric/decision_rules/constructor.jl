# cubic is default

struct Linear end
struct Cubic end

struct CompletePolnomial{order} end
struct Smolyak end

# incomplete type inference:
DecisionRule(exo_grid, endo_grid, i::Int64) = DecisionRule(exo_grid, endo_grid, Val{i})
DecisionRule(exo_grid, endo_grid, i::Int64, interp_type) = DecisionRule(exo_grid, endo_grid, Val{i}, interp_type)


# cubic interpolation
function DecisionRule(exo_grid::EmptyGrid, endo_grid::UCGrid, ::Type{Val{nx}}, ::Type{Cubic}) where nx
    CubicDR(exo_grid, endo_grid, Val{nx})
end

function DecisionRule(exo_grid::UnstructuredGrid, endo_grid::UCGrid, ::Type{Val{nx}}, ::Type{Cubic}) where nx
    CubicDR(exo_grid, endo_grid, Val{nx})
end

function DecisionRule(exo_grid::UCGrid, endo_grid::UCGrid, ::Type{Val{nx}}, ::Type{Cubic}) where nx
    CubicDR(exo_grid, endo_grid, Val{nx})
end

# Smolyak interpolation
function DecisionRule(
        exo_grid::S, endo_grid::SmolyakGrid, ::Type{Val{nx}}, ::Type{Smolyak}
    ) where S <: Union{EmptyGrid,<:UnstructuredGrid} where nx
    SmolyakDR(exo_grid, endo_grid, Val{nx})
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
    grid = ProductGrid(exo_grid, endo_grid)
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

mutable struct CachedDecisionRule{T,S}
    dr::T
    process::S
end

CachedDecisionRule(process::AbstractDiscretizedProcess, grid::Grid, n_x::Int) =
    CachedDecisionRule(DecisionRule(process.grid, grid, n_x), process)

CachedDecisionRule(process::AbstractDiscretizedProcess, grid::Grid, values) =
    CachedDecisionRule(DecisionRule(process.grid, grid, values), process)

set_values!(cdr::CachedDecisionRule, v) = set_values!(cdr.dr, v)

# defaults

(cdr::CachedDecisionRule)(v::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr,v)
(cdr::CachedDecisionRule)(m::Union{<:Point,<:ListOfPoints}, s::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr, m, s)

# Interpolating rules
(cdr::CachedDecisionRule{<:AbstractDecisionRule{T,<:Grid}, <:AbstractDiscretizedProcess})(i::Int, s::Union{<:Point,<:ListOfPoints}) where T<:Union{UCGrid, RandomGrid} = evaluate(cdr.dr,node(Point, cdr.process, i), s)

(cdr::CachedDecisionRule{<:AbstractDecisionRule{T,<:Grid}, <:AbstractDiscretizedProcess})(i::Int, j::Int, s::Union{<:Point,<:ListOfPoints}) where T<:Union{UCGrid, RandomGrid} = evaluate(cdr.dr,inode(Point, cdr.process, i, j), s)


(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:EmptyGrid,<:Grid}, <:AbstractDiscretizedProcess})(i::Int, j::Int, s::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr, s)

(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:EmptyGrid,<:Grid}, <:AbstractDiscretizedProcess})(i::Int, s::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr, s)

# maybe keep only the i,j,s calls.

(cdr::CachedDecisionRule{<:AbstractDecisionRule{T,<:Grid}, <:DiscreteMarkovProcess})(i::Int, s::Union{<:Point,<:ListOfPoints}) where T<:Union{UCGrid, RandomGrid} = evaluate(cdr.dr,node(Point, cdr.process, i), s)
(cdr::CachedDecisionRule{<:AbstractDecisionRule{T,<:Grid}, <:DiscreteMarkovProcess})(i::Int, j::Int, s::Union{<:Point,<:ListOfPoints}) where T<:Union{UCGrid, RandomGrid} = evaluate(cdr.dr,inode(Point, cdr.process, i, j), s)

(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:Grid}, <:DiscreteMarkovProcess})(i::Int, s::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr, i, s)
(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:Grid}, <:DiscreteMarkovProcess})(i::Int, j::Int, s::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr, j, s)

(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:Grid}, <:DiscreteMarkovProcess})(::Val{(0,3)}, i::Int, s::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr, Val((0,2)), i, s)
(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:Grid}, <:DiscreteMarkovProcess})(::Val{(0,3)}, i::Int, j::Int, s::Union{<:Point,<:ListOfPoints}) = evaluate(cdr.dr, Val((0,2)), j, s)


#
# # We don't do compat for CachedDR cuz it's a nigthmare.
#
# (cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:Grid}, <:DiscreteMarkovProcess})(i::Int, s::Union{Vector{Float64},Matrix{Float64}}) = evaluate(cdr.dr, i, s)
# (cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:Grid}, <:DiscreteMarkovProcess})(i::Int, j::Int, s::Union{Vector{Float64},Matrix{Float64}}) = evaluate(cdr.dr, j, s)
#
#
# (cdr::CachedDecisionRule)(v::Union{Vector{Float64},Matrix{Float64}}) = evaluate(cdr.dr,v)
# (cdr::CachedDecisionRule)(m::Union{Vector{Float64},Matrix{Float64}}, s::Union{Vector{Float64},Matrix{Float64}}) = evaluate(cdr.dr, m, s)
# (cdr::CachedDecisionRule)(i::Int, s::Union{Vector{Float64},Matrix{Float64}}) = evaluate(cdr.dr,node(cdr.process, i), s)
# (cdr::CachedDecisionRule)(i::Int, j::Int, s::Union{Vector{Float64},Matrix{Float64}}) = evaluate(cdr.dr,inode(cdr.process, i, j), s)
