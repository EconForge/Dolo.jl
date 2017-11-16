struct CubicDR{S<:Grid, T<:CartesianGrid, nx, d} <: AbstractDecisionRule{S, T, nx}
    grid_exo::S
    grid_endo::T
    itp::Vector{Array{Value{nx}, d}}
end

(dr::CubicDR)(args...) = evaluate(dr, args...)

#####
##### 1-argument decision rule
#####

function CubicDR(exo_grid::EmptyGrid, endo_grid::CartesianGrid{d}, i::Union{Val{nx}, Type{Val{nx}}}) where d where nx
    dims = endo_grid.n+2
    c = [zeros(Value{nx},dims...)]
    CubicDR{EmptyGrid, CartesianGrid{d}, nx, d}(exo_grid, endo_grid, c)
end

function set_values!(dr::CubicDR{EmptyGrid,CartesianGrid{d},n_x,d},  V::Vector{<:Array{Value{n_x}}}) where n_x where d
    n = dr.grid_endo.n
    C = dr.itp[1]
    ind = [2:(n[i]+1) for i=1:length(n)]
    # this segfaults, don't know why
    # C[ind...] = reshape(V[1], n)
    C[ind...] = V[1]
    prefilter!(C)
end

function evaluate(dr::CubicDR{EmptyGrid,CartesianGrid{d}}, points::Vector{Point{d}}) where d
    a = SVector{d,Float64}(dr.grid_endo.min)
    b = SVector{d,Float64}(dr.grid_endo.max)
    n = SVector{d,Int64}(dr.grid_endo.n)
    C = dr.itp[1]
    return eval_UC_spline(a, b, n, C, points)
end

####
#### 2 CartesianGrid continous arguments d.r.
####

# function CubicDR(exo_grid::EmptyGrid, endo_grid::CartesianGrid{d}, i::Union{Val{nx}, Type{Val{nx}}}) where d where nx
#     dims = endo_grid.n+2
#     c = [zeros(Value{nx},dims...)]
#     CubicDR{EmptyGrid, CartesianGrid, nx, d}(exo_grid, endo_grid, c)
# end

function CubicDR(exo_grid::CartesianGrid{d1}, endo_grid::CartesianGrid{d2}, i::Union{Val{nx}, Type{Val{nx}}}) where d1 where d2 where nx
    dims = cat(1, exo_grid.n+2, endo_grid.n+2)
    c = [zeros(Value{nx},dims...)]
    CubicDR{CartesianGrid{d1}, CartesianGrid{d2}, nx, d1+d2}(exo_grid, endo_grid, c)
end
#
function set_values!(dr::CubicDR{CartesianGrid{d1},CartesianGrid{d2},n_x,d}, V::Vector{<:Array{Value{n_x}}}) where n_x where d where d1 where d2
    exog = dr.grid_exo
    endog = dr.grid_endo

    data = reshape( cat(1, V...), endog.n..., exog.n...)
    d_m = ndims(exog)
    d_s = ndims(endog)
    perms = cat(1, ((d_s+1):(d_s+d_m))..., (1:d_s)...)
    RV = permutedims(data, perms)
    C = dr.itp[1]
    dims = size(C)
    inds = [2:(i-1) for i in dims]
    C[inds...] = RV
    prefilter!(C)
end


function evaluate(dr::CubicDR{CartesianGrid{d1},CartesianGrid{d2}, n_x, d}, z::ListOfPoints{d}) where n_x where d where d1 where d2
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    n = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    cc = dr.itp[1]
    res = eval_UC_spline(a, b, n, cc, z)
    return res
end

function evaluate(dr::CubicDR{CartesianGrid{d1},CartesianGrid{d2}}, x::ListOfPoints{d1}, y::ListOfPoints{d2}) where d1 where d2
    N = length(x)
    z = [ [x[i]; y[i]] for i=1:N]
    evaluate(dr, z)
end

function evaluate(dr::CubicDR{CartesianGrid{d1},CartesianGrid{d2}}, x::Point{d1}, y::ListOfPoints{d2}) where d1 where d2
    N = length(y)
    z = [ [x; y[i]] for i=1:N]
    evaluate(dr, z)
end

function evaluate(dr::CubicDR{CartesianGrid{d1},CartesianGrid{d2}}, i::Int64, y::Union{Point{d2},ListOfPoints{d2}}) where d1 where d2
    x = node(Point, dr.grid_exo, i)
    evaluate(dr, x, y)
end

# TODO replace by generic call
function evaluate(dr::CubicDR{CartesianGrid{d1},CartesianGrid{d2}}, x::Point{d1}, y::Point{d2}) where d1 where d2
    dr([x],[y])[1]
end

####
#### UnstructuredGrid × CartesianGrid 2 continous arguments d.r.
####

function CubicDR(exo_grid::UnstructuredGrid{d1}, endo_grid::CartesianGrid{d2}, i::Union{Val{nx}, Type{Val{nx}}}) where d1 where d2 where nx
    c = [zeros(Value{nx},(endo_grid.n+2)...) for _=1:n_nodes(exo_grid)]
    CubicDR{UnstructuredGrid{d1}, CartesianGrid{d2}, nx, d2}(exo_grid, endo_grid, c)
end

function set_values!(dr::CubicDR{UnstructuredGrid{d1},CartesianGrid{d2}, n_x, d2}, values::Vector{<:Array{Value{n_x}}}) where d1 where d2 where n_x
    orders = dr.grid_endo.n
    inds = [2:(o+1) for o in orders]
    for (i,C) in enumerate(dr.itp)
        C[inds...] = reshape(values[i],orders...)
        prefilter!(C)
    end
end

function evaluate(dr::CubicDR{<:UnstructuredGrid,<:CartesianGrid}, i::Int, z::Vector{Point{d}}) where d
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[i]
    splines.eval_UC_spline(a, b, n, cc, z)
end

# TODO replace by generic call
function evaluate(dr::CubicDR{<:UnstructuredGrid,<:CartesianGrid}, i::Int, z::Point{d}) where d
    evaluate(dr, i, [z])[1]
end
