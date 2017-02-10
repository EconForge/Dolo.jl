import splines

abstract AbstractDecisionRule
abstract AbtsractCubicSplineDR <: AbstractDecisionRule

# temporary

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

# Constant Decision Rule
type ConstantDecisionRule <: AbstractDecisionRule
    values::Array{Float64,1}
end
(dr::ConstantDecisionRule)(x::Array{Float64,1}) = dr.values
(dr::ConstantDecisionRule)(x::Array{Float64,2}) = repmat(dr.values', size(x, 1), 1)
(dr::ConstantDecisionRule)(i::Int, x::Union{Vector,Matrix}) =  dr(x)
(dr::ConstantDecisionRule)(i::Int, j::Int, x::Union{Vector, Matrix}) = dr(x)


## Cubic spline decision rules

function set_values!(dr::AbtsractCubicSplineDR, values::Array{Array{Float64,2},1})
    grid = dr.grid
    a = grid.min
    b = grid.max
    orders = grid.n
    dr.coefficients = [filter_mcoeffs(a, b, orders, v) for v in values]
end

function set_values!(dr::AbtsractCubicSplineDR, values::Array{Float64,2})
    set_values!(dr, [values])
end

# Decision Rule on Markov Process

type MCDecisionRule <: AbtsractCubicSplineDR
    process::DiscreteProcess
    grid::CartesianGrid
    n_x::Int # number of values
    coefficients::Array{Array{Float64}}
end

function DecisionRule(mc::DiscreteMarkovProcess, grid::CartesianGrid, n_x::Int)
    a = grid.min
    b = grid.max
    orders = grid.n
    n_ms = size(mc.values, 1)
    coeffs = [zeros(n_x, (orders+2)...) for v in 1:n_ms]
    return MCDecisionRule(mc, grid, n_x,  coeffs)
end

function DecisionRule(mc::DiscreteMarkovProcess, grid::CartesianGrid, values::Array{Array{Float64,2}})
    n_x = size(values[1], 2)
    dr = DecisionRule(mc, grid, n_x)
    set_values!(dr, values)
    return dr
end

function (dr::MCDecisionRule)(i::Int, x::Array{Float64,2})
    a = dr.grid.min
    b = dr.grid.max
    n = dr.grid.n
    cc = dr.coefficients[i]
    res = splines.eval_UC_multi_spline(a, b, n, cc, x)'
    return res
end

(dr::MCDecisionRule)(i::Int, x::Array{Float64,1}) =  dr(i, x')[:]
(dr::MCDecisionRule)(i::Int, j::Int, x::Union{Vector, Matrix}) = dr(j, x)


# Decision Rule on IID Process

type IDecisionRule <: AbtsractCubicSplineDR
    process::IIDExogenous
    grid::CartesianGrid
    n_x::Int # number of values
    coefficients::Vector{Array{Float64}}
end

function DecisionRule(proc::MvNormal, grid::CartesianGrid, n_x::Int)
    a = grid.min
    b = grid.max
    orders = grid.n
    coeffs = [zeros(n_x, (orders+2)...)]
    return IDecisionRule(proc, grid, n_x,  coeffs)
end

function DecisionRule(proc::MvNormal, grid::CartesianGrid, values::Array{Array{Float64,2},1})
    n_x = size(values[1], 2)
    dr = DecisionRule(proc, grid, n_x)
    set_values!(dr, values)
    return dr
end

function DecisionRule(proc::MvNormal, grid::CartesianGrid, values::Array{Float64,2})
    return DecisionRule(proc, grid, [values])
end

function (dr::IDecisionRule)(x::Array{Float64,2})
    a = dr.grid.min
    b = dr.grid.max
    n = dr.grid.n
    cc = dr.coefficients[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, x)'
    return res
end

(dr::IDecisionRule)(x::Array{Float64,1}) = dr(x')[:]
(dr::IDecisionRule)(i::Int, x::Union{Vector, Matrix}) = dr(x)
(dr::IDecisionRule)(i::Int, j::Int, x::Union{Vector, Matrix}) = dr(x)


### drs on VAR1

# Decision Rule on ContinuousProcess Process


type CPDecisionRule <: AbtsractCubicSplineDR
    process::DiscretizedProcess
    grid_exo::CartesianGrid
    grid_endo::CartesianGrid
    grid::CartesianGrid
    n_x::Int # number of values
    coefficients::Array{Array{Float64}}
end

function DecisionRule(dp::DiscretizedProcess, grid_endo::CartesianGrid, n_x::Int)
    # this assumes the discretized process
    # was discretized on a cartesian grid
    a_exo = minimum(dp.nodes, 1)[:]
    b_exo = maximum(dp.nodes, 1)[:]
    o_exo = Int64[length(unique(dp.nodes[:, i])) for i=1:size(dp.nodes, 2)]
    grid_exo = CartesianGrid(a_exo, b_exo, o_exo)
    grid = CartesianGrid(
            cat(1, grid_exo.min, grid_endo.min),
            cat(1, grid_exo.max, grid_endo.max),
            cat(1, grid_exo.n, grid_endo.n)
        )
    a = grid.min
    b = grid.max
    orders = grid.n
    # n_ms = size(dp.values, 1)
    coeffs = [zeros(n_x, (orders+2)...)]
    return CPDecisionRule(dp, grid_exo, grid_endo, grid, n_x, coeffs)
end

function DecisionRule(dp::DiscretizedProcess, grid::CartesianGrid, values::Array{Array{Float64,2}})
    n_x = size(values[1], 2)
    dr = DecisionRule(dp, grid, n_x)
    set_values!(dr, values)
    return dr
end

function (dr::CPDecisionRule)(z::Array{Float64,2})
    a = dr.grid.min
    b = dr.grid.max
    n = dr.grid.n
    cc = dr.coefficients[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

(dr::CPDecisionRule)(z::Vector) = dr(z')[:]


function (dr::CPDecisionRule)(i::Int64, x::Array{Float64,2})
    a = dr.grid.min
    b = dr.grid.max
    n = dr.grid.n
    cc = dr.coefficients[1]
    N = size(x, 1)
    e = dr.process.nodes[i, :]
    z = cat(2, repmat(e', N, 1), x)
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end



function (dr::CPDecisionRule)(i::Int64, j::Int64, x::Array{Float64,2})
    a = dr.grid.min
    b = dr.grid.max
    n = dr.grid.n
    cc = dr.coefficients[1]
    N = size(x, 1)
    e = dr.process.integration_nodes[i][j, :]
    z = cat(2, repmat(e', N, 1), x)
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end



(dr::CPDecisionRule)(i::Int64, x::Vector) = dr(i, x')[:]
(dr::CPDecisionRule)(i::Int64, j::Int64, x::Vector) = dr(i, j,x')[:]
