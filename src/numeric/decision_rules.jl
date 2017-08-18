@compat abstract type AbstractDecisionRule{S<:Grid,T<:Grid} end

function Base.show(io::IO, dr::AbstractDecisionRule)
    println(io, typeof(dr))
end

# abstract AbstractCachedDecisionRule{S,T} <: AbstractDecisionRule{S,T}
# we don't implement that using subtypes anymore

# ---------------------- #
# Constant decision rule #
# ---------------------- #

@compat type ConstantDecisionRule <: AbstractDecisionRule{EmptyGrid,EmptyGrid}
    constants::Vector{Float64}
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

@compat type BiTaylorExpansion <: AbstractDecisionRule{EmptyGrid,EmptyGrid}
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

# -------------------------- #
# Cubic spline DecisionRules #
# -------------------------- #

type DecisionRule{S<:Grid,T<:Grid,Titp} <: AbstractDecisionRule{S,T}
    grid_exo::S
    grid_endo::T
    n_x::Int # number of values (could be a method if there was a nice method name)
    itp::Titp
end

function set_values!(dr::DecisionRule{EmptyGrid}, values::Array{Float64,2})
    set_values!(dr, [values])
end

@compat const CubicSplineDR{S<:Grid,T<:Grid} = DecisionRule{S,T,<:Vector{<:Array{Float64}}}
@compat const SmolyakDR{S<:Grid} = DecisionRule{S,SmolyakGrid,Vector{Matrix{Float64}}}
@compat const CompletePolyDR{S<:Grid} = DecisionRule{S,RandomGrid,Vector{Matrix{Float64}}}

#####
##### 1-argument decision rule
#####

## Splines!

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, n_x::Int)
    orders = grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...)]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::CubicSplineDR{EmptyGrid, CartesianGrid}, values::Array{Array{Float64,2},1})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    orders = dr.grid_endo.n
    dr.itp = [filter_mcoeffs(a, b, orders, v) for v in values]
end

function evaluate(dr::CubicSplineDR{EmptyGrid, CartesianGrid}, z::AbstractMatrix)
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

## Smolyak!

function DecisionRule(grid_exo::EmptyGrid, grid_endo::SmolyakGrid, n_x::Int)
    coeffs = [Array{Float64}(n_nodes(grid_endo), n_x)]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::SmolyakDR{EmptyGrid}, values::Vector{Matrix{Float64}})
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], dr.grid_endo.B_nodes, values[i])
    end
end

function evaluate(dr::SmolyakDR{EmptyGrid}, z::AbstractMatrix)
    B = BM.evalbase(dr.grid_endo.smol_params, z)
    B*dr.itp[1]
end

## Complete polynomials

function DecisionRule(grid_exo::EmptyGrid, grid_endo::RandomGrid, n_x::Int)
    n_s = size(nodes(grid_endo), 2)
    coeffs = [Array{Float64}(BM.n_complete(n_s, 3), n_x)]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::CompletePolyDR{EmptyGrid}, values::Vector{Matrix{Float64}})
    B_grid = BM.complete_polynomial(nodes(dr.grid_endo), 3)
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], B_grid, values[i])
    end
end

function evaluate(dr::CompletePolyDR{EmptyGrid}, z::AbstractMatrix)
    B = BM.complete_polynomial(z, 3)
    B*dr.itp[1]
end

## Common routines for grid_exo <: EmptyGrid
function DecisionRule(grid_exo::EmptyGrid, grid_endo, values::Vector{Matrix{Float64}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

(dr::DecisionRule{EmptyGrid})(z::AbstractMatrix) = evaluate(dr, z)
(dr::DecisionRule{EmptyGrid})(z::AbstractVector) = vec(dr(z'))
(dr::DecisionRule{EmptyGrid})(i::Int, x::Union{AbstractVector,AbstractMatrix}) = dr(x)
(dr::DecisionRule{EmptyGrid})(x::Union{AbstractVector,AbstractMatrix}, y::Union{AbstractVector,AbstractMatrix}) = dr(y)


####
#### 2 continous arguments d.r.
####

function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = [grid_exo.n; grid_endo.n]
    coeffs = [zeros(n_x, (orders+2)...)]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end
#
function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, values::Array{Array{Float64,2},1})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(dr::Dolo.AbstractDecisionRule{Dolo.CartesianGrid, Dolo.CartesianGrid}, values::Array{Array{Float64,2},1})
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    orders = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    n_x = size(values[1], 2)
    vv = zeros(prod(dr.grid_exo.n), prod(dr.grid_endo.n), n_x)
    for n=1:size(vv, 1)
        # ah, if only Julia had chosen C-order arrays !
        vv[n, :, :] = values[n]
    end
    vv = reshape(vv, prod(dr.grid_endo.n)*prod(dr.grid_exo.n), n_x)
    dr.itp = [Dolo.filter_mcoeffs(a, b, orders, vv)]
end

function evaluate(dr::AbstractDecisionRule{CartesianGrid,CartesianGrid}, z::AbstractMatrix)
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    n = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    cc = dr.itp[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end


(dr::DecisionRule{CartesianGrid, CartesianGrid})(z::AbstractMatrix) = evaluate(dr, z)
(dr::DecisionRule{CartesianGrid, CartesianGrid})(z::AbstractVector) = dr(z')[:]
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::AbstractVector, y::AbstractVector) = dr(cat(1, x, y))
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::AbstractMatrix, y::AbstractMatrix) = dr([x y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::AbstractVector, y::AbstractMatrix) = dr([repmat(x', size(y, 1), 1) y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(i::Int, y::Union{AbstractVector,AbstractMatrix}) = dr(node(dr.grid_exo, i), y)


####
#### 2 continous arguments d.r.
####

## Splines!

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...) for i in 1:n_nodes(grid_exo)]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::CubicSplineDR{UnstructuredGrid,CartesianGrid}, values::Array{Matrix{Float64},1})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    orders = dr.grid_endo.n
    dr.itp = [filter_mcoeffs(a, b, orders, vals) for vals in values]
end


function evaluate(dr::CubicSplineDR{UnstructuredGrid,CartesianGrid}, i::Int, z::AbstractMatrix)
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[i]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

## Smolyak!

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::SmolyakGrid, n_x::Int)
    coeffs = [Array{Float64}(n_nodes(grid_endo), n_x) for i in 1:n_nodes(grid_exo)]
    DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::SmolyakDR{UnstructuredGrid}, values::Vector{Matrix{Float64}})
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], dr.grid_endo.B_nodes, values[i])
    end
end

function evaluate(dr::SmolyakDR{UnstructuredGrid}, i::Int, z::AbstractMatrix)
    BM.evalbase(dr.grid_endo.smol_params, z) * dr.itp[i]
end

## Complete polynomials

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::RandomGrid, n_x::Int)
    n_s = size(nodes(grid_endo), 2)
    coeffs = [
        Array{Float64}(BM.n_complete(n_s, 3), n_x)
        for i in 1:n_nodes(grid_exo)
    ]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::CompletePolyDR{UnstructuredGrid}, values::Vector{Matrix{Float64}})
    B_grid = BM.complete_polynomial(nodes(dr.grid_endo), 3)
    for i in 1:length(values)
        A_ldiv_B!(dr.itp[i], B_grid, values[i])
    end
end

function evaluate(dr::CompletePolyDR{UnstructuredGrid}, i::Int, z::AbstractMatrix)
    B = BM.complete_polynomial(z, 3)
    B*dr.itp[i]
end

## Common routines for grid_exo <: UnstructuredGrid
function DecisionRule(grid_exo::UnstructuredGrid, grid_endo, values::Array{Array{Float64,2}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

(dr::DecisionRule{UnstructuredGrid})(i::Int, y::AbstractMatrix) = evaluate(dr, i, y)
# (dr::DecisionRule{UnstructuredGrid})(i::AbstractMatrix{Int}, y::AbstractMatrix) =
# vcat( [evaluate(dr, i[j, 1], y[j, :]) for j=1:size(y, 1)]... )

(dr::DecisionRule{UnstructuredGrid})(i::Int, y::AbstractVector) = dr(i, reshape(y, 1, length(y)))[:]
(dr::DecisionRule{UnstructuredGrid})(i::AbstractVector{Int}, y::AbstractMatrix) =
  vcat( [dr(i[j], y[j, :])' for j=1:size(y, 1)]... )
(dr::DecisionRule{UnstructuredGrid})(i::AbstractVector{Int}, y::AbstractVector) =
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

@compat (cdr::CachedDecisionRule{<:DecisionRule{UnstructuredGrid}, DiscreteMarkovProcess})(i::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(i, s)
@compat (cdr::CachedDecisionRule{<:DecisionRule{UnstructuredGrid}, DiscreteMarkovProcess})(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(j, s)

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
