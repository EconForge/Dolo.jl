@compat abstract type AbstractDecisionRule{S<:Grid,T<:Grid,nx} end

function Base.show(io::IO, dr::AbstractDecisionRule)
    println(io, typeof(dr))
end

# ---------------------- #
# Constant decision rule #
# ---------------------- #

@compat type ConstantDecisionRule{nx} <: AbstractDecisionRule{EmptyGrid,EmptyGrid,nx}
    constants::Vector{Float64}
end

function ConstantDecisionRule(constants::Vector{Float64})
    nx = length(constants)
    ConstantDecisionRule{nx}(constants)
end

(dr::ConstantDecisionRule)(x::AbstractVector) = dr.constants
(dr::ConstantDecisionRule)(x::AbstractMatrix) = repmat(dr.constants', size(x, 1), 1)
(dr::ConstantDecisionRule)(x::AbstractVector, y::AbstractVector) = dr.constants
(dr::ConstantDecisionRule)(x::AbstractVector, y::AbstractMatrix) = repmat( dr.constants', size(y, 1), 1)
(dr::ConstantDecisionRule)(x::AbstractMatrix, y::AbstractVector) = repmat( dr.constants', size(x, 1), 1)
(dr::ConstantDecisionRule)(x::AbstractMatrix, y::AbstractMatrix) = repmat( dr.constants', size(x, 1), 1)
(dr::ConstantDecisionRule)(i::Int, x::Union{AbstractVector,AbstractMatrix}) = dr(x)
(dr::ConstantDecisionRule)(i::Int, j::Int, x::Union{AbstractVector,AbstractMatrix}) = dr(x)

# ------------------------------ #
# 2-dimensional Taylor Expansion #
# ------------------------------ #

@compat type BiTaylorExpansion{nx} <: AbstractDecisionRule{EmptyGrid,EmptyGrid,nx}
    m0::Vector{Float64}
    s0::Vector{Float64}
    x0::Vector{Float64}
    x_m::Matrix{Float64}
    x_s::Matrix{Float64}
end

(dr::BiTaylorExpansion)(m::AbstractVector, s::AbstractVector) = dr.x0 + dr.x_m*(m-dr.m0) + dr.x_s*(s-dr.s0)
(dr::BiTaylorExpansion)(m::AbstractMatrix, s::AbstractVector) = vcat([(dr(m[i, :], s))' for i=1:size(m, 1) ]...)
(dr::BiTaylorExpansion)(m::AbstractVector, s::AbstractMatrix) = vcat([(dr(m, s[i, :]))' for i=1:size(s, 1) ]...)
(dr::BiTaylorExpansion)(m::AbstractMatrix, s::AbstractMatrix) = vcat([(dr(m[i, :], s[i, :]))' for i=1:size(m, 1) ]...)

# ------------ #
# DecisionRule #
# ------------ #

type DecisionRule{S<:Grid,T<:Grid,nx,Titp} <: AbstractDecisionRule{S,T,nx}
    grid_exo::S
    grid_endo::T
    itp::Titp

    function (::Type{DecisionRule{S,T,nx,Titp}}){S,T,Titp,nx}(grid_exo::S, grid_endo::T, itp::Titp)
        new{S,T,nx,Titp}(grid_exo, grid_endo, itp)
    end
end

function DecisionRule{S,T,nx,Titp}(grid_exo::S, grid_endo::T, itp::Titp, ::Union{Val{nx},Type{Val{nx}}})
    DecisionRule{typeof(grid_exo),typeof(grid_endo),nx,typeof(itp)}(grid_exo, grid_endo, itp)
end

function set_values!(dr::T, values::Array{Float64,2}) where T <: DecisionRule{<:EmptyGrid}
    set_values!(dr, [values])
end

## Common routines for grid_exo <: EmptyGrid
function DecisionRule(grid_exo::EmptyGrid, grid_endo, values::Vector{Matrix{Float64}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, Val{n_x})
    set_values!(dr, values)
    return dr
end

(dr::DecisionRule{<:EmptyGrid})(z::AbstractMatrix) = evaluate(dr, z)
(dr::DecisionRule{<:EmptyGrid})(z::AbstractVector) = vec(dr(z'))
(dr::DecisionRule{<:EmptyGrid})(i::Int, x::Union{AbstractVector,AbstractMatrix}) = dr(x)
(dr::DecisionRule{<:EmptyGrid})(x::Union{AbstractVector,AbstractMatrix}, y::Union{AbstractVector,AbstractMatrix}) = dr(y)

## Common routines for grid_exo <: UnstructuredGrid
function DecisionRule(grid_exo::UnstructuredGrid, grid_endo, values::Array{Array{Float64,2}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, Val{n_x})
    set_values!(dr, values)
    return dr
end

(dr::DecisionRule{<:UnstructuredGrid})(i::Int, y::AbstractMatrix) = evaluate(dr, i, y)
# (dr::DecisionRule{<:UnstructuredGrid})(i::AbstractMatrix{Int}, y::AbstractMatrix) =
# vcat( [evaluate(dr, i[j, 1], y[j, :]) for j=1:size(y, 1)]... )

(dr::DecisionRule{<:UnstructuredGrid})(i::Int, y::AbstractVector) = dr(i, reshape(y, 1, length(y)))[:]
(dr::DecisionRule{<:UnstructuredGrid})(i::AbstractVector{Int}, y::AbstractMatrix) =
  vcat( [dr(i[j], y[j, :])' for j=1:size(y, 1)]... )
(dr::DecisionRule{<:UnstructuredGrid})(i::AbstractVector{Int}, y::AbstractVector) =
    vcat( [dr(i[j], y[j, :])' for j=1:size(y, 1)]... )


#####
##### Cached Decision Rules (do we really need them ?)
#####
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

(cdr::CachedDecisionRule)(v::Union{AbstractVector,AbstractMatrix}) = cdr.dr(v)
(cdr::CachedDecisionRule)(m::Union{AbstractVector,AbstractMatrix}, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(m, s)
(cdr::CachedDecisionRule)(i::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(node(cdr.process, i), s)
(cdr::CachedDecisionRule)(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(inode(cdr.process, i, j), s)

@compat (cdr::CachedDecisionRule{<:DecisionRule{<:UnstructuredGrid}, DiscreteMarkovProcess})(i::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(i, s)
@compat (cdr::CachedDecisionRule{<:DecisionRule{<:UnstructuredGrid}, DiscreteMarkovProcess})(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(j, s)

# --------------- #
# Helper function #
# --------------- #

function filter_mcoeffs(a::Array{Float64,1}, b::Array{Float64,1}, n::Array{Int,1}, mvalues::Array{Float64})
    n_x = size(mvalues)[end]
    vals = reshape(mvalues, n..., n_x)
    coeffs = zeros(n_x, (n+2)...)
    ii = [Colon() for i=1:(ndims(vals)-1)]
    for i_x in 1:n_x
        tmp = splines.filter_coeffs(a, b, n, vals[ii..., i_x])
        coeffs[i_x, ii...] = tmp
    end
    return coeffs
end
