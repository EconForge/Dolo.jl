
@compat abstract type AbstractDecisionRule{S,T} end


# abstract AbstractCachedDecisionRule{S,T} <: AbstractDecisionRule{S,T}
# we don't implement that using subtypes anymore

type DecisionRule{S,T} <: AbstractDecisionRule{S,T}
    grid_exo::S
    grid_endo::T
    n_x::Int # number of values (could be a method if there was a nice method name)
    coefficients::Array
end

function Base.show(io::IO, dr::AbstractDecisionRule)
    println(typeof(dr))
end


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
(dr::BiTaylorExpansion)(i::Int64, s::AbstractMatrix) = dr(s)


#####
##### 1-argument decision rule
#####

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, n_x::Int)
    orders = grid_endo.n
    coeffs = zeros(SVector{n_x,Float64}, (orders+2)...)
    return DecisionRule{EmptyGrid,CartesianGrid}(grid_exo, grid_endo, n_x, coeffs)
end

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, values::Vector{Point{n_x}}) where n_x
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(ddr::Dolo.DecisionRule{EmptyGrid,CartesianGrid}, V::Array{Value{n_x},d}) where n_x where d
    n = ddr.grid_endo.n
    C = ddr.coefficients
    ind = [2:(n[i]+1) for i=1:length(n)]
    C[ind...] = V
    prefilter!(C)
end

function evaluate(dr::Dolo.DecisionRule{EmptyGrid,CartesianGrid},points::Vector{Point{d}}) where d
    a = SVector{d,Float64}(dr.grid_endo.min)
    b = SVector{d,Float64}(dr.grid_endo.max)
    n = SVector{d,Int64}(dr.grid_endo.n)
    C = dr.coefficients
    return eval_UC_spline(a,b,n,C,points)
end

(dr::DecisionRule{EmptyGrid, CartesianGrid})(z::ListOfPoints{d}) where d = evaluate(dr, z)
(dr::DecisionRule{EmptyGrid, CartesianGrid})(z::Point{d}) where d = evaluate(dr, [z])[1]


# compatibility with multiple exogenous d.r.

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, values::Vector{Vector{Point{n_x}}}) where n_x
    DecisionRule(grid_exo, grid_endo, values[1])
end

function set_values!(ddr::Dolo.DecisionRule{EmptyGrid,CartesianGrid}, V::Vector{Vector{Value{n_x}}}) where n_x
    set_values!(ddr, V[1])
end



# backward compatibility

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, values::Vector{Matrix{Float64}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(ddr::AbstractDecisionRule{EmptyGrid,CartesianGrid}, values::Matrix{Float64})
    n = ddr.grid_endo.n
    n_x = size(values,2)
    V = reshape(values, n..., n_x)
    perm = cat(1,[ndims(V)],1:(length(n)))
    W = reinterpret(SVector{n_x,Float64}, permutedims(V, perm), tuple(n...))
    set_values!(ddr, W)
end

function set_values!(dr::AbstractDecisionRule{EmptyGrid, CartesianGrid}, values::Vector{Matrix{Float64}})
    set_values!(dr, values[1])
end

function evaluate(dr::AbstractDecisionRule{EmptyGrid, CartesianGrid}, z::AbstractMatrix)
    N,d = size(z)
    zz = reinterpret(Point{d},z',(N,))
    out = evaluate(dr, zz)
    n_x = length(out[1])
    reinterpret(Float64, out, (n_x,N))'
end

(dr::DecisionRule{EmptyGrid, CartesianGrid})(z::AbstractMatrix) = evaluate(dr, z)
(dr::DecisionRule{EmptyGrid, CartesianGrid})(z::AbstractVector) = dr(z')[:]
(dr::DecisionRule{EmptyGrid, CartesianGrid})(i::Int, x::Union{AbstractVector,AbstractMatrix}) = dr(x)
(dr::DecisionRule{EmptyGrid, CartesianGrid})(x::Union{AbstractVector,AbstractMatrix}, y::Union{AbstractVector,AbstractMatrix}) = dr(y)

####
#### 2 continous arguments d.r.
####

function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = [grid_exo.n; grid_endo.n]
    coeffs = zeros(Point{n_x}, (orders+2)...)
    return DecisionRule{CartesianGrid, CartesianGrid}(grid_exo, grid_endo, n_x, coeffs)
end
#

function set_values!(dr::Dolo.AbstractDecisionRule{Dolo.CartesianGrid, Dolo.CartesianGrid}, V::Array{Value{n_x},d}) where n_x where d
    C = dr.coefficients
    dims = size(C)
    inds = [2:(i-1) for i in dims]
    C[inds...] = V
    prefilter!(C)
end


function evaluate(dr::AbstractDecisionRule{CartesianGrid,CartesianGrid}, z::ListOfPoints{d}) where d
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    n = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    cc = dr.coefficients
    res = splines.eval_UC_spline(a, b, n, cc, z)
    return res
end


# backward compat

function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, values::Vector{Matrix{Float64}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(dr::Dolo.AbstractDecisionRule{Dolo.CartesianGrid, Dolo.CartesianGrid}, values::Vector{Matrix{Float64}})
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    orders = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    n_x = size(values[1], 2)
    vv = zeros(prod(dr.grid_exo.n), prod(dr.grid_endo.n), n_x)
    for n=1:size(vv, 1)
        # ah, if only Julia had chosen C-order arrays !
        vv[n, :, :] = values[n]
    end
    N = prod(dr.grid_endo.n)*prod(dr.grid_exo.n)
    vv = reshape(vv, prod(dr.grid_endo.n)*prod(dr.grid_exo.n), n_x)
    perm = cat(1,[ndims(vv)],1:(ndims(vv)-1))
    vv = permutedims(vv, perm)
    V = reinterpret(Value{n_x}, vv, (N,))
    set_values!(dr, V)
end

function evaluate(dr::AbstractDecisionRule{CartesianGrid,CartesianGrid}, z::AbstractMatrix)
    N,ns = size(z)
    zz = reinterpret(Point{ns}, z', (N,))
    res = evaluate(dr, zz)
    n_x = length(res[1])
    reinterpret(Float64, res, (n_x, N))'
end

(dr::DecisionRule{CartesianGrid, CartesianGrid})(z::AbstractMatrix) = evaluate(dr, z)
(dr::DecisionRule{CartesianGrid, CartesianGrid})(z::AbstractVector) = dr(vector_to_matrix(z))[:]
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::AbstractVector, y::AbstractVector) = dr(cat(1, x, y))
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::AbstractMatrix, y::AbstractMatrix) = dr([x y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::AbstractVector, y::AbstractMatrix) = dr([repmat(vector_to_matrix(x'), size(y, 1), 1) y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(i::Int, y::Union{AbstractVector,AbstractMatrix}) = dr(node(dr.grid_exo, i), y)


####
#### markov chain exogenous process d.r.
####

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = grid_endo.n
    coeffs = [zeros(Point{n_x}, (orders+2)...) for i in 1:n_nodes(grid_exo)]
    return (DecisionRule{UnstructuredGrid, CartesianGrid})(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::AbstractDecisionRule{UnstructuredGrid,CartesianGrid}, values::Vector{Vector{Point{n_x}}}) where n_x
    orders = dr.grid_endo.n
    inds = [2:(o+1) for o in orders]
    for (i,C) in enumerate(dr.coefficients)
        C[inds...] = reshape(values[i],orders...)
        prefilter!(C)
    end
end

function evaluate(dr::AbstractDecisionRule{UnstructuredGrid,CartesianGrid}, i::Int, z::Vector{Point{d}}) where d
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.coefficients[i]
    res = splines.eval_UC_spline(a, b, n, cc, z)
    return res
end


# backward compatibility

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, values::Vector{Matrix{Float64}})
    n_x = size(values[1],2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(dr::AbstractDecisionRule{UnstructuredGrid,CartesianGrid}, values::Vector{Matrix{Float64}})
    N,n_x = size(values[1])
    vals = [reinterpret(Point{n_x},v',(N,)) for v in values]
    set_values!(dr, vals)
end

function evaluate(dr::AbstractDecisionRule{UnstructuredGrid,CartesianGrid}, i::Int, z::AbstractMatrix)
    N,d = size(z)
    zz = reinterpret(Point{d},z',(N,))
    oo = evaluate(dr,i,zz)
    n_x = length(oo[1])
    o = reinterpret(Float64, oo, (n_x,N))
    return o'
end

(dr::DecisionRule{UnstructuredGrid, CartesianGrid})(i::Int, y::AbstractMatrix) = evaluate(dr, i, y)
(dr::DecisionRule{UnstructuredGrid, CartesianGrid})(i::Int, y::AbstractVector) = dr(i, vector_to_matrix(y))[:]
(dr::DecisionRule{UnstructuredGrid, CartesianGrid})(i::AbstractVector{Int}, y::AbstractMatrix) = vcat( [dr(i[j], y[j, :])' for j=1:size(y, 1)]... )

(dr::DecisionRule{UnstructuredGrid, CartesianGrid})(i::AbstractVector{Int}, y::AbstractVector) = vcat( [dr(i[j], y[j, :])' for j=1:size(y, 1)]... )

#####
##### Cached Decision Rules (do we really need them ?)
#####

type CachedDecisionRule{T,S}
    dr::T
    process::S
end

@compat const AbstractADecisionRule = Union{DecisionRule,CachedDecisionRule}

CachedDecisionRule(process::AbstractDiscretizedProcess, grid::Grid, n_x::Int) =
    CachedDecisionRule(DecisionRule(process.grid, grid, n_x), process)

CachedDecisionRule(process::AbstractDiscretizedProcess, grid::Grid, values) =
    CachedDecisionRule(DecisionRule(process.grid, grid, values), process)

set_values!(cdr::CachedDecisionRule, v) = set_values!(cdr.dr, v)

# defaults

(cdr::CachedDecisionRule)(v::Point{d}) where d = cdr.dr(v)
(cdr::CachedDecisionRule)(v::ListOfPoints{d}) where d = cdr.dr(v)
(cdr::CachedDecisionRule)(i::Int, v::ListOfPoints{d}) where d = cdr.dr(node(cdr.process, i), v)
(cdr::CachedDecisionRule)(i::Int, j::Int, v::ListOfPoints{d}) where d = cdr.dr(inode(cdr.process, i, j), v)


(cdr::CachedDecisionRule)(v::Union{AbstractVector,AbstractMatrix}) = cdr.dr(v)
(cdr::CachedDecisionRule)(m::Union{AbstractVector,AbstractMatrix}, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(m, s)
(cdr::CachedDecisionRule)(i::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(node(cdr.process, i), s)
(cdr::CachedDecisionRule)(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(inode(cdr.process, i, j), s)


(cdr::CachedDecisionRule{DecisionRule{UnstructuredGrid, CartesianGrid}, DiscreteMarkovProcess})(i::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(i, s)
(cdr::CachedDecisionRule{DecisionRule{UnstructuredGrid, CartesianGrid}, DiscreteMarkovProcess})(i::Int, j::Int, s::Union{AbstractVector,AbstractMatrix}) = cdr.dr(j, s)
