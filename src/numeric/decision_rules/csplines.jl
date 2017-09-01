@compat const CubicSplineDR{S<:Grid,T<:Grid,nx} = DecisionRule{S,T,nx,<:Array{<:Point}}

#####
##### 1-argument decision rule
#####

function DecisionRule{nx}(grid_exo::EmptyGrid, grid_endo::CartesianGrid, ::Union{Val{nx},Type{Val{nx}}})
    orders = grid_endo.n
    coeffs = zeros(SVector{n_x,Float64}, (orders+2)...)
    return DecisionRule(grid_exo, grid_endo, coeffs, Val{nx})
end

function set_values!(dr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid},  V::Array{Value{n_x},d}) where n_x where d
    n = ddr.grid_endo.n
    C = ddr.itp
    ind = [2:(n[i]+1) for i=1:length(n)]
    C[ind...] = V
    prefilter!(C)
end

function evaluate(dr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid}, points::Vector{Point{d}}) where d
    a = SVector{d,Float64}(dr.grid_endo.min)
    b = SVector{d,Float64}(dr.grid_endo.max)
    n = SVector{d,Int64}(dr.grid_endo.n)
    C = dr.itp
    return eval_UC_spline(a, b, n, C, points)
end

# compatibility with multiple exogenous d.r.

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, values::Vector{Vector{Point{n_x}}}) where n_x
    DecisionRule(grid_exo, grid_endo, values[1])
end

function set_values!(ddr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid}, V::Vector{Vector{Value{n_x}}}) where n_x
    set_values!(ddr, V[1])
end

# backward compatibility

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, values::Vector{Matrix{Float64}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(ddr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid}, values::Matrix{Float64})
    n = ddr.grid_endo.n
    n_x = size(values,2)
    V = reshape(values, n..., n_x)
    perm = cat(1,[ndims(V)],1:(length(n)))
    W = reinterpret(SVector{n_x,Float64}, permutedims(V, perm), tuple(n...))
    set_values!(ddr, W)
end

function set_values!(dr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid}, values::Vector{Matrix{Float64}})
    set_values!(dr, values[1])
end

function evaluate(dr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid}, z::AbstractMatrix)
    N,d = size(z)
    zz = reinterpret(Point{d},z',(N,))
    out = evaluate(dr, zz)
    n_x = length(out[1])
    reinterpret(Float64, out, (n_x,N))'
end

####
#### 2 CartesianGrid continous arguments d.r.
####

function DecisionRule{nx}(grid_exo::CartesianGrid, grid_endo::CartesianGrid, ::Union{Val{nx},Type{Val{nx}}})
    # hmm kind of silently assuming we have cartesian grid
    orders = [grid_exo.n; grid_endo.n]
    coeffs = zeros(Point{n_x}, (orders+2)...)
    return DecisionRule(grid_exo, grid_endo, coeffs, Val{nx})
end

function set_values!(dr::CubicSplineDR{<:CartesianGrid,<:CartesianGrid}, V::Array{Value{n_x},d}) where n_x where d
    C = dr.itp
    dims = size(C)
    inds = [2:(i-1) for i in dims]
    C[inds...] = V
    prefilter!(C)
end

function evaluate(dr::CubicSplineDR{<:CartesianGrid,<:CartesianGrid}, z::ListOfPoints{d}) where d
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    n = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    cc = dr.itp
    res = splines.eval_UC_spline(a, b, n, cc, z)
    return res
end

## COMPAT
function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, values::Vector{Matrix{Float64}})
    nx = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, Val{nx})
    set_values!(dr, values)
    return dr
end

function set_values!(dr::CubicSplineDR{<:CartesianGrid,<:CartesianGrid,n_x}, values::Vector{Matrix{Float64}})
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    orders = cat(1, dr.grid_exo.n, dr.grid_endo.n)
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

function evaluate(dr::CubicSplineDR{<:CartesianGrid,<:CartesianGrid,n_x}, z::AbstractMatrix)
    N,ns = size(z)
    zz = reinterpret(Point{ns}, z', (N,))
    res = evaluate(dr, zz)
    reinterpret(Float64, res, (n_x, N))'
end


(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(z::AbstractMatrix) = evaluate(dr, z)
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(z::AbstractVector) = vec(dr(vector_to_matrix(z)))
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(x::AbstractVector, y::AbstractVector) = dr(cat(1, x, y))
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(x::AbstractMatrix, y::AbstractMatrix) = dr([x y])
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(x::AbstractVector, y::AbstractMatrix) = dr([vector_to_matrix(x'), size(y, 1), 1) y])
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(i::Int, y::Union{AbstractVector,AbstractMatrix}) = dr(node(dr.grid_exo, i), y)

####
#### UnstructuredGrid × CartesianGrid 2 continous arguments d.r.
####

function DecisionRule{nx}(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, ::Union{Val{nx},Type{Val{nx}}})
    # hmm kind of silently assuming we have cartesian grid
    orders = grid_endo.n
    coeffs = [zeros(Point{n_x}, (orders+2)...) for i in 1:n_nodes(grid_exo)]
    return DecisionRule(grid_exo, grid_endo, coeffs, Val{nx})
end

function set_values!(dr::CubicSplineDR{<:UnstructuredGrid,<:CartesianGrid}, values::Vector{Vector{Point{n_x}}}) where n_x
    orders = dr.grid_endo.n
    inds = [2:(o+1) for o in orders]
    for (i,C) in enumerate(dr.itp)
        C[inds...] = reshape(values[i],orders...)
        prefilter!(C)
    end
end

function evaluate(dr::CubicSplineDR{<:UnstructuredGrid,<:CartesianGrid}, i::Int, z::Vector{Point{d}}) where d
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[i]
    splines.eval_UC_spline(a, b, n, cc, z)
end

# COMPAT

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, values::Vector{Matrix{Float64}})
    n_x = size(values[1],2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end

function set_values!(dr::CubicSplineDR{<:UnstructuredGrid,<:CartesianGrid,n_x}, values::Vector{Matrix{Float64}})
    N = size(values[1], 1)
    vals = [reinterpret(Point{n_x},v',(N,)) for v in values]
    set_values!(dr, vals)
end

function evaluate(dr::CubicSplineDR{<:UnstructuredGrid,<:CartesianGrid,n_x}, i::Int, z::AbstractMatrix)
    N,d = size(z)
    zz = reinterpret(Point{d},z',(N,))
    oo = evaluate(dr,i,zz)
    o = reinterpret(Float64, oo, (n_x, N))
    return o'
end
