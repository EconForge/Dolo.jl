abstract AbstractDecisionRule{S,T}

# abstract AbstractCachedDecisionRule{S,T} <: AbstractDecisionRule{S,T}
# we don't implement that using subtypes anymore

type DecisionRule{S,T} <: AbstractDecisionRule{S,T}
    grid_exo::S
    grid_endo::T
    n_x::Int # number of values (could be a method if there was a nice method name)
    coefficients::Array{Array{Float64}}
end


type ConstantDecisionRule
    constants::Vector{Float64}
end

(dr::ConstantDecisionRule)(x::Vector{Float64}) = dr.constants
(dr::ConstantDecisionRule)(x::Matrix{Float64}) = repmat(dr.constants',size(x,1),1)
(dr::ConstantDecisionRule)(x::Vector{Float64}, y::Vector{Float64}) = dr.constants
(dr::ConstantDecisionRule)(x::Vector{Float64}, y::Matrix{Float64}) = repmat( dr.constants', size(y,1), 1)
(dr::ConstantDecisionRule)(x::Matrix{Float64}, y::Vector{Float64}) = repmat( dr.constants', size(x,1), 1)
(dr::ConstantDecisionRule)(x::Matrix{Float64}, y::Matrix{Float64}) = repmat( dr.constants', size(x,1), 1)
(dr::ConstantDecisionRule)(i::Int64, x::Union{Vector{Float64},Matrix{Float64}}) = dr(x)
(dr::ConstantDecisionRule)(i::Int64, j::Int64, x::Union{Vector{Float64},Matrix{Float64}}) = dr(x)


type BiTaylorExpansion <: AbstractDecisionRule
    m0::Vector{Float64}
    s0::Vector{Float64}
    x0::Vector{Float64}
    x_m::Matrix{Float64}
    x_s::Matrix{Float64}
end

(dr::BiTaylorExpansion)(m::Vector{Float64}, s::Vector{Float64}) = dr.x0 + dr.x_m*(m-dr.m0) + dr.x_s*(s-dr.s0)
(dr::BiTaylorExpansion)(m::Matrix{Float64}, s::Vector{Float64}) = vcat( [ (dr(m[i,:],s))' for i=1:size(m,1) ]...)
(dr::BiTaylorExpansion)(m::Vector{Float64}, s::Matrix{Float64}) = vcat( [ (dr(m,s[i,:]))' for i=1:size(s,1) ]...)
(dr::BiTaylorExpansion)(m::Matrix{Float64}, s::Matrix{Float64}) = vcat( [ (dr(m[i,:],s[i,:]))' for i=1:size(m,1) ]...)


function filter_mcoeffs(a::Array{Float64,1},b::Array{Float64,1},n::Array{Int64,1},mvalues::Array{Float64})
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



#####
##### 1-argument decision rule
#####

using splines

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, n_x::Int)
    orders = grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...)]
    return DecisionRule{EmptyGrid, CartesianGrid}(grid_exo, grid_endo, n_x, coeffs)
end

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, values::Array{Array{Float64,2}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(dr::AbstractDecisionRule{EmptyGrid, CartesianGrid}, values::Array{Array{Float64,2},1})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    orders = dr.grid_endo.n
    dr.coefficients = [filter_mcoeffs(a, b, orders, v) for v in values]
end

function set_values!(dr::AbstractDecisionRule{EmptyGrid, CartesianGrid}, values::Array{Float64,2})
    set_values!(dr, [values])
end

function evaluate(dr::AbstractDecisionRule{EmptyGrid, CartesianGrid},z::Matrix{Float64})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.coefficients[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

(dr::DecisionRule{EmptyGrid, CartesianGrid})(z::Matrix{Float64}) = evaluate(dr,z)
(dr::DecisionRule{EmptyGrid, CartesianGrid})(z::Vector{Float64}) = dr(z')[:]
(dr::DecisionRule{EmptyGrid, CartesianGrid})(i::Int64, x::Union{Vector{Float64},Matrix{Float64}}) = dr(x)
(dr::DecisionRule{EmptyGrid, CartesianGrid})(x::Union{Vector{Float64},Matrix{Float64}}, y::Union{Vector{Float64},Matrix{Float64}}) = dr(y)


####
#### 2 continous arguments d.r.
####

function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = grid_exo.n + grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...)]
    return (DecisionRule{CartesianGrid, CartesianGrid})(grid_exo, grid_endo, n_x, coeffs)
end
#
function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, values::Array{Array{Float64,2},1})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(dr::Dolo.AbstractDecisionRule{Dolo.CartesianGrid, Dolo.CartesianGrid}, values::Array{Array{Float64,2},1})
    a = cat(1,dr.grid_exo.min,dr.grid_endo.min)
    b = cat(1,dr.grid_exo.max,dr.grid_endo.max)
    orders = cat(1,dr.grid_exo.n,dr.grid_endo.n)
    n_x = size(values[1],2)
    vv = zeros(prod(dr.grid_exo.n), prod(dr.grid_endo.n), n_x)
    for n=1:size(vv,1)
        # ah, if only Julia had chosen C-order arrays !
        vv[n,:,:] = values[n]
    end
    vv = reshape(vv, prod(dr.grid_endo.n)*prod(dr.grid_exo.n),n_x)
    dr.coefficients = [Dolo.filter_mcoeffs(a,b,orders,vv)]
end

function evaluate(dr::AbstractDecisionRule{CartesianGrid, CartesianGrid}, z::Matrix{Float64})
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    n = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    cc = dr.coefficients[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

(dr::DecisionRule{CartesianGrid, CartesianGrid})(z::Matrix{Float64}) = evaluate(dr,z)
(dr::DecisionRule{CartesianGrid, CartesianGrid})(z::Vector{Float64}) = dr(z')[:]
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::Vector{Float64},y::Vector{Float64}) = dr(cat(1,x,y))
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::Matrix{Float64},y::Matrix{Float64}) = dr([x y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::Vector{Float64},y::Matrix{Float64}) = dr([repmat(x',size(y,1),1) y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(i::Int64,y::Union{Vector{Float64},Matrix{Float64}}) = dr(node(dr.grid_exo,i),y)


####
#### 2 continous arguments d.r.
####

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...) for i in 1:n_nodes(grid_exo)]
    return (DecisionRule{UnstructuredGrid, CartesianGrid})(grid_exo, grid_endo, n_x, coeffs)
end
#
function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, values::Array{Array{Float64,2}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end


function set_values!(dr::AbstractDecisionRule{UnstructuredGrid, CartesianGrid}, values::Array{Matrix{Float64},1})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    orders = dr.grid_endo.n
    dr.coefficients = [filter_mcoeffs(a, b, orders, vals) for vals in values]
end


function evaluate(dr::AbstractDecisionRule{UnstructuredGrid, CartesianGrid}, i::Int64, z::Matrix{Float64})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.coefficients[i]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

(dr::DecisionRule{UnstructuredGrid, CartesianGrid})(i::Int64,y::Matrix{Float64}) = evaluate(dr,i,y)
(dr::DecisionRule{UnstructuredGrid, CartesianGrid})(i::Int64,y::Vector{Float64}) = dr(dr,i,y')[:]




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

set_values!(cdr::CachedDecisionRule, v) = set_values!(cdr.dr,v)

# defaults

(cdr::CachedDecisionRule)(v::Union{Vector{Float64},Matrix{Float64}}) = cdr.dr(v)
(cdr::CachedDecisionRule)(m::Union{Vector{Float64},Matrix{Float64}},s::Union{Vector{Float64},Matrix{Float64}}) = cdr.dr(m,s)
(cdr::CachedDecisionRule)(i::Int64,s::Union{Vector{Float64},Matrix{Float64}}) = cdr.dr(node(cdr.process,i),s)
(cdr::CachedDecisionRule)(i::Int64,j::Int64,s::Union{Vector{Float64},Matrix{Float64}}) = cdr.dr(inode(cdr.process,i,j),s)

(cdr::CachedDecisionRule{DecisionRule{UnstructuredGrid, CartesianGrid}, DiscreteMarkovProcess})(i::Int64,s::Union{Vector{Float64},Matrix{Float64}}) = cdr.dr(i,s)
(cdr::CachedDecisionRule{DecisionRule{UnstructuredGrid, CartesianGrid}, DiscreteMarkovProcess})(i::Int64,j::Int64,s::Union{Vector{Float64},Matrix{Float64}}) = cdr.dr(j,s)
