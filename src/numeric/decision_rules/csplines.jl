@compat const CubicSplineDR{S<:Grid,T<:Grid} = DecisionRule{S,T,<:Vector{<:Array{Float64}}}

#####
##### 1-argument decision rule
#####

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, n_x::Int)
    orders = grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...)]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid}, values::Array{Array{Float64,2},1})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    orders = dr.grid_endo.n
    dr.itp = [filter_mcoeffs(a, b, orders, v) for v in values]
end

function evaluate(dr::CubicSplineDR{<:EmptyGrid,<:CartesianGrid}, z::AbstractMatrix)
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

####
#### 2 CartesianGrid continous arguments d.r.
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

function set_values!(dr::CubicSplineDR{<:CartesianGrid, <:CartesianGrid}, values::Array{Array{Float64,2},1})
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

function evaluate(dr::CubicSplineDR{<:CartesianGrid,<:CartesianGrid}, z::AbstractMatrix)
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    n = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    cc = dr.itp[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end


(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(z::AbstractMatrix) = evaluate(dr, z)
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(z::AbstractVector) = dr(z')[:]
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(x::AbstractVector, y::AbstractVector) = dr(cat(1, x, y))
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(x::AbstractMatrix, y::AbstractMatrix) = dr([x y])
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(x::AbstractVector, y::AbstractMatrix) = dr([repmat(x', size(y, 1), 1) y])
(dr::DecisionRule{<:CartesianGrid, <:CartesianGrid})(i::Int, y::Union{AbstractVector,AbstractMatrix}) = dr(node(dr.grid_exo, i), y)

####
#### UnstructuredGrid × CartesianGrid 2 continous arguments d.r.
####

function DecisionRule(grid_exo::UnstructuredGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...) for i in 1:n_nodes(grid_exo)]
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs)
end

function set_values!(dr::CubicSplineDR{<:UnstructuredGrid,<:CartesianGrid}, values::Array{Matrix{Float64},1})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    orders = dr.grid_endo.n
    dr.itp = [filter_mcoeffs(a, b, orders, vals) for vals in values]
end

function evaluate(dr::CubicSplineDR{<:UnstructuredGrid,<:CartesianGrid}, i::Int, z::AbstractMatrix)
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.itp[i]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end
