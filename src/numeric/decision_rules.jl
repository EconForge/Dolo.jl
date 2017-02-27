
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


abstract AbstractDecisionRule{S,T}
abstract AbstractCachedDecisionRule{S,T} <: AbstractDecisionRule{S,T}

type DecisionRule{S,T} <: AbstractDecisionRule{S,T}
    grid_exo::S
    grid_endo::T
    n_x::Int # number of values
    coefficients::Array{Array{Float64}}
end

type CachedDecisionRule{S,T} <: AbstractCachedDecisionRule{S,T}
    grid_exo::S
    grid_endo::T
    n_x::Int # number of values
    coefficients::Array{Array{Float64}}
    process::AbstractDiscretizedProcess
end

function CachedDecisionRule(process::DiscretizedIIDProcess, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    grid_exo = EmptyGrid()
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs, process)
end

function CachedDecisionRule(process::DiscretizedProcess, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assume we have cartesian grid
    grid_exo = process.grid
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs, process)
end

function CachedDecisionRule(process::DiscreteMarkovProcess, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    grid_exo = UnstructuredGrid(process.nodes)
    return DecisionRule(grid_exo, grid_endo, n_x, coeffs, process)
end


#
# function DecisionRule(process::DiscretizedIIDProcess, grid_endo::CartesianGrid, n_x::Int)
#     # hmm kind of silently assuming we have cartesian grid
#     grid_exo = EmptyGrid()
#     return DecisionRule{EmptyGrid, CartesianGrid}(grid_exo, grid_endo, n_x, coeffs, process)
# end
#
#
# function DecisionRule(process::DiscretizedIIDProcess, grid_endo::CartesianGrid, values::Array{Array{Float64,2}})
#     n_x = size(values[1], 2)
#     dr = DecisionRule(grid_exo, grid_endo, n_x)
#     set_values!(dr, values)
#     return dr
# end

function DecisionRule(grid_exo::EmptyGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
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


function set_values!(dr::DecisionRule{EmptyGrid, CartesianGrid}, values::Array{Array{Float64,2},1})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    orders = dr.grid_endo.n
    dr.coefficients = [filter_mcoeffs(a, b, orders, v) for v in values]
end

function set_values!(dr::DecisionRule{EmptyGrid, CartesianGrid}, values::Array{Float64,2})
    set_values!(dr, [values])
end

function (dr::DecisionRule{EmptyGrid, CartesianGrid})(z::Matrix{Float64})
    a = dr.grid_endo.min
    b = dr.grid_endo.max
    n = dr.grid_endo.n
    cc = dr.coefficients[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

(dr::DecisionRule{EmptyGrid, CartesianGrid})(z::Vector{Float64}) = dr(z')[:]
(dr::DecisionRule{EmptyGrid, CartesianGrid})(i::Int64, x::Union{Vector{Float64},Matrix{Float64}}) = dr(x)
(dr::DecisionRule{EmptyGrid, CartesianGrid})(i::Int64, j::Int64, x::Union{Vector{Float64},Matrix{Float64}}) = dr(x)


#
function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, n_x::Int)
    # hmm kind of silently assuming we have cartesian grid
    orders = grid_exo.n + grid_endo.n
    coeffs = [zeros(n_x, (orders+2)...)]
    return (DecisionRule{CartesianGrid, CartesianGrid})(grid_exo, grid_endo, n_x, coeffs)
end
#
function DecisionRule(grid_exo::CartesianGrid, grid_endo::CartesianGrid, values::Array{Array{Float64,2}})
    n_x = size(values[1], 2)
    dr = DecisionRule(grid_exo, grid_endo, n_x)
    set_values!(dr, values)
    return dr
end


function set_values!(dr::DecisionRule{CartesianGrid, CartesianGrid}, values::Matrix{Float64})
    a = cat(1,dr.grid_exo.min,dr.grid_endo.min)
    b = cat(1,dr.grid_exo.max,dr.grid_endo.max)
    orders = cat(1,dr.grid_exo.n,dr.grid_endo.n)
    dr.coefficients = [filter_mcoeffs(a, b, orders, values)]
end

function set_values!(dr::DecisionRule{CartesianGrid, CartesianGrid}, values::Array{Matrix{Float64},1})
    set_values!(dr, values[1])
end

function (dr::DecisionRule{CartesianGrid, CartesianGrid})(z::Matrix{Float64})
    a = cat(1, dr.grid_exo.min, dr.grid_endo.min)
    b = cat(1, dr.grid_exo.max, dr.grid_endo.max)
    n = cat(1, dr.grid_exo.n, dr.grid_endo.n)
    cc = dr.coefficients[1]
    res = splines.eval_UC_multi_spline(a, b, n, cc, z)'
    return res
end

(dr::DecisionRule{CartesianGrid, CartesianGrid})(z::Vector{Float64}) = dr(z')[:]
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::Vector{Float64},y::Vector{Float64}) = dr(cat(1,x,y))
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::Matrix{Float64},y::Matrix{Float64}) = dr([x y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(x::Vector{Float64},y::Matrix{Float64}) = dr([repmat(x',size(y,1),1) y])
(dr::DecisionRule{CartesianGrid, CartesianGrid})(i::Int64,y::Union{Vector{Float64},Matrix{Float64}}) = dr(node(dr.grid_exo,i),y)

# (dr::DecisionRule{CartesianGrid, CartesianGrid})(i::Int64, x::Union{Vector{Float64},Matrix{Float64}}) = dr(x)
# # (dr::DecisionRule{EmptyGrid, CartesianGrid})(i::Int64, x::Matrix{Float64}) = dr(x)
# (dr::DecisionRule{CartesianGrid, CartesianGrid})(i::Int64, j::Int64, x::Union{Vector{Float64},Matrix{Float64}}) = dr(x)
