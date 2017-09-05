# cubic is default

abstract type Linear end
abstract type Cubic end

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

CachedDecisionRule(process::AbstractDiscretizedProcess, grid::Grid, n_x::Int) =
    CachedDecisionRule(DecisionRule(process.grid, grid, n_x), process)

CachedDecisionRule(process::AbstractDiscretizedProcess, grid::Grid, values) =
    CachedDecisionRule(DecisionRule(process.grid, grid, values), process)

set_values!(cdr::CachedDecisionRule, v) = set_values!(cdr.dr, v)

# defaults

(cdr::CachedDecisionRule)(v::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,v)
(cdr::CachedDecisionRule)(m::Union{AbstractVector,AbstractMatrix}, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr, m, s)
(cdr::CachedDecisionRule)(i::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,node(cdr.process, i), s)
(cdr::CachedDecisionRule)(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,inode(cdr.process, i, j), s)

(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:CartesianGrid}, <:DiscreteMarkovProcess})(i::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,i, s)
(cdr::CachedDecisionRule{<:AbstractDecisionRule{<:UnstructuredGrid,<:CartesianGrid}, <:DiscreteMarkovProcess})(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = evaluate(cdr.dr,j, s)
