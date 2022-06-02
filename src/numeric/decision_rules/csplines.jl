struct CubicDR{S<:Grid, T<:UCGrid, nx, d} <: AbstractDecisionRule{S, T, nx}
    grid::ProductGrid{S, T}
    grid_exo::S
    grid_endo::T
    itp::Vector{Array{Value{nx}, d}}
end


(dr::CubicDR)(args...) = evaluate(dr, args...)

#####
##### 1-argument decision rule
#####

function CubicDR(exo_grid::EmptyGrid{d1}, endo_grid::UCGrid{d2}, i::Union{Val{nx}, Type{Val{nx}}}) where d2 where d1 where nx
    dims = endo_grid.n .+ 2
    c = [zeros(Value{nx},dims...)]
    grid = ProductGrid(exo_grid, endo_grid)
    CubicDR{EmptyGrid{d1}, UCGrid{d2}, nx, d2}(grid, exo_grid, endo_grid, c)
end

function set_values!(dr::CubicDR{EmptyGrid{d1},UCGrid{d},n_x,d},  V::MSM{Value{n_x}}) where n_x where d where d1
    set_values!(dr, V.views)
end


function set_values!(dr::CubicDR{EmptyGrid{d1},UCGrid{d},n_x,d},  V::AbstractVector{<:AbstractVector{Value{n_x}}}) where n_x where d where d1
    n = dr.grid_endo.n
    C = dr.itp[1]
    ind = [2:(n[i]+1) for i=1:length(n)]
    # this segfaults, don't know why
    # C[ind...] = reshape(V[1], n)
    C[ind...] = V[1]
    prefilter!(C)
end

function evaluate(dr::CubicDR{EmptyGrid{d1},UCGrid{d}}, points::AbstractVector{Point{d}}) where d where d1
    a = SVector{d,Float64}(dr.grid_endo.min)
    b = SVector{d,Float64}(dr.grid_endo.max)
    n = SVector{d,Int64}(dr.grid_endo.n)
    C = dr.itp[1]
    return eval_UC_spline(a, b, n, C, points)
end

#
function evaluate(dr::CubicDR{EmptyGrid{d1},UCGrid{d}}, y::Point{d}) where d where d1
    return evaluate(dr, [y])[1]
end

function evaluate(dr::CubicDR{EmptyGrid{d1},UCGrid{d}}, y::Vector{Float64}) where d where d1
    p = evaluate(dr, Point{d}(y...))
    return [p...]
end

evaluate(dr::CubicDR{EmptyGrid{d1},UCGrid{d}}, i, y) where d where d1 = evaluate(dr, y)

####
#### 2 CartesianGrid continous arguments d.r.
####

# function CubicDR(exo_grid::EmptyGrid, endo_grid::CartesianGrid{d}, i::Union{Val{nx}, Type{Val{nx}}}) where d where nx
#     dims = endo_grid.n+2
#     c = [zeros(Value{nx},dims...)]
#     CubicDR{EmptyGrid, CartesianGrid, nx, d}(exo_grid, endo_grid, c)
# end

function CubicDR(exo_grid::UCGrid{d1}, endo_grid::UCGrid{d2}, i::Union{Val{nx}, Type{Val{nx}}}) where d1 where d2 where nx
    grid = ProductGrid(exo_grid, endo_grid)
    dims = cat(exo_grid.n .+ 2, endo_grid.n .+ 2; dims=1)
    c = [zeros(Value{nx},dims...)]
    CubicDR{UCGrid{d1}, UCGrid{d2}, nx, d1+d2}(grid, exo_grid, endo_grid, c)
end
#
function set_values!(dr::CubicDR{UCGrid{d1},UCGrid{d2},n_x,d}, V::MSM{Value{n_x}}) where n_x where d where d1 where d2
    set_values!(dr, V.views)
end

function set_values!(dr::CubicDR{UCGrid{d1},UCGrid{d2},n_x,d}, V::AbstractVector{<:AbstractVector{Value{n_x}}}) where n_x where d where d1 where d2
    exog = dr.grid_exo
    endog = dr.grid_endo

    data = reshape( cat(V...; dims=1), endog.n..., exog.n...)
    d_m = ndims(exog)
    d_s = ndims(endog)
    perms = cat(((d_s+1):(d_s+d_m))..., (1:d_s)...; dims=1)
    RV = permutedims(data, perms)
    C = dr.itp[1]
    dims = size(C)
    inds = [2:(i-1) for i in dims]
    C[inds...] = RV
    prefilter!(C)
end


function evaluate(dr::CubicDR{UCGrid{d1},UCGrid{d2}, n_x, d}, z::AbstractVector{Point{d}}) where n_x where d where d1 where d2
    a = cat(dr.grid_exo.min, dr.grid_endo.min; dims=1)
    b = cat(dr.grid_exo.max, dr.grid_endo.max; dims=1)
    n = cat(dr.grid_exo.n, dr.grid_endo.n; dims=1)
    cc = dr.itp[1]
    res = eval_UC_spline(a, b, n, cc, z)
    return res
end

function evaluate(dr::CubicDR{UCGrid{d1},UCGrid{d2}}, x::AbstractVector{Point{d1}}, y::AbstractVector{Point{d2}}) where d1 where d2
    N = length(x)
    z = [ [x[i]; y[i]] for i=1:N]
    evaluate(dr, z)
end

function evaluate(dr::CubicDR{UCGrid{d1},UCGrid{d2}}, x::Point{d1}, y::AbstractVector{Point{d2}}) where d1 where d2
    N = length(y)
    z = [ [x; y[i]] for i=1:N]
    evaluate(dr, z)
end

function evaluate(dr::CubicDR{UCGrid{d1},UCGrid{d2}}, i::Int64, y::Union{Point{d2},AbstractVector{Point{d2}}}) where d1 where d2
    x = node(Point, dr.grid_exo, i)
    evaluate(dr, x, y)
end

# TODO replace by generic call?
function evaluate(dr::CubicDR{UCGrid{d1},UCGrid{d2}}, x::Point{d1}, y::Point{d2}) where d1 where d2
    dr([x],[y])[1]
end

####
#### UnstructuredGrid × CartesianGrid 2 continous arguments d.r.
####

function CubicDR(exo_grid::UnstructuredGrid{d1}, endo_grid::UCGrid{d2}, i::Union{Val{nx}, Type{Val{nx}}}) where d1 where d2 where nx
    c = [zeros(Value{nx},(endo_grid.n.+2)...) for _=1:n_nodes(exo_grid)]
    grid = ProductGrid(exo_grid, endo_grid)
    CubicDR{UnstructuredGrid{d1}, UCGrid{d2}, nx, d2}(grid, exo_grid, endo_grid, c)
end

function set_values!(dr::CubicDR{UnstructuredGrid{d1},UCGrid{d2}, n_x, d2}, V::MSM{Value{n_x}}) where d1 where d2 where n_x
    set_values!(dr, V.views)
end

function set_values!(dr::CubicDR{UnstructuredGrid{d1},UCGrid{d2}, n_x, d2}, values::AbstractVector{<:AbstractVector{Value{n_x}}}) where d1 where d2 where n_x
    orders = dr.grid_endo.n
    inds = [2:(o+1) for o in orders]
    for (i,C) in enumerate(dr.itp)
        C[inds...] = reshape(values[i],orders...)
        prefilter!(C)
    end
end

function evaluate(dr::CubicDR{<:UnstructuredGrid,<:UCGrid}, i::Int, z::AbstractVector{Point{d}}) where d
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[i]
    splines.eval_UC_spline(a, b, n, cc, z)
end

# TODO replace by generic call?
function evaluate(dr::CubicDR{<:UnstructuredGrid,<:UCGrid}, i::Int, z::SVector{d, U}) where d where U
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[i]
    splines.eval_UC_spline(a, b, n, cc, z)
end


# Experimental

function evaluate(dr::CubicDR{<:UnstructuredGrid,<:UCGrid}, ::Val{(0,2)}, i::Int, z::Vector{SVector{d, U}}) where d where U
    it = [evaluate(dr, Val((0,2)), i, u) for u in z]
    return (
        [e[1] for e in it],
        [e[2] for e in it]
    )
end


function evaluate(dr::CubicDR{<:UnstructuredGrid,<:UCGrid}, ::Val{(0,2)}, i::Int, z::SVector{d, U}) where d where U
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[i]
    f = u->splines.eval_UC_spline(a, b, n, cc, u)
    return (
        f(z),
        ForwardDiff.jacobian(f, z)
    )
    
end
