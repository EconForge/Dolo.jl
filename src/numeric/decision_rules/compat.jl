## Common:


### backward compat

function DecisionRule(exo_grid, endo_grid, values::Vector{Matrix{Float64}})
    N,n_x = size(values[1])
    vals = [to_LOP(v) for v in values]
    DecisionRule(exo_grid, endo_grid, vals)
end

function DecisionRule(exo_grid, endo_grid, values::Vector{Matrix{Float64}}, tt)
    N,n_x = size(values[1])
    vals = [reinterpret(Value{n_x}, v', (N,)) for v in values]
    DecisionRule(exo_grid, endo_grid, vals, tt)
end

# EmptyGrid x CartesianGrid

function evaluate(dr::AbstractDecisionRule{EmptyGrid, CartesianGrid{d}, n_x}, z::AbstractMatrix{Float64}) where n_x where d
    N = size(z,1)
    assert(size(z,2)==d)
    points = to_LOP(Point{d}, z)
    out = evaluate(dr,points)
    return reinterpret(Float64, out, (n_x,N))'
end

evaluate(dr::AbstractDecisionRule{EmptyGrid, <:CartesianGrid}, z::Vector) = vec(evaluate(dr,z'))
evaluate(dr::AbstractDecisionRule{EmptyGrid, <:CartesianGrid}, i::Int, y::Union{Vector,AbstractMatrix}) = evaluate(dr,y) # kind of nonsensical
evaluate(dr::AbstractDecisionRule{EmptyGrid, <:CartesianGrid}, x::Union{Vector,AbstractMatrix}, y::Union{Vector,AbstractMatrix}) = evaluate(dr,y)

function set_values!(dr::AbstractDecisionRule{EmptyGrid, <:CartesianGrid, n_x}, vals::Vector{Matrix{Float64}}) where n_x
    N = size(vals[1],1)
    dims = tuple(dr.grid_endo.n...)
    V = [reinterpret(Value{n_x}, vals[1]', dims)]
    set_values!(dr, V)
end

# CartesianGrid x CartesianGrid

function evaluate(dr::AbstractDecisionRule{CartesianGrid{d1}, CartesianGrid{d2}, n_x}, x::AbstractMatrix, y::AbstractMatrix) where d1 where d2 where n_x
    N = size(x,1)
    assert(size(y,1)==N)
    xx = to_LOP(Point{d1}, x)
    yy = to_LOP(Point{d2}, y)
    res = evaluate(dr,xx,yy)
    return reinterpret(Float64, res, (n_x, N))'
end
#
evaluate(dr::AbstractDecisionRule{<:CartesianGrid, <:CartesianGrid}, x::Vector{Float64}, y::Matrix{Float64}) = evaluate(dr, repmat(x', size(y,1)), y)
evaluate(dr::AbstractDecisionRule{<:CartesianGrid, <:CartesianGrid}, x::Vector{Float64}, y::Vector{Float64}) = evaluate(dr, vector_to_matrix(x), vector_to_matrix(y))
evaluate(dr::AbstractDecisionRule{CartesianGrid{d1}, CartesianGrid{d2}}, z::AbstractMatrix) where d1 where d2 = evaluate(dr, [z[:,1:d1] z[:,1:(d1+1:end)]])
evaluate(dr::AbstractDecisionRule{<:CartesianGrid, <:CartesianGrid}, z::Vector{Float64}) = vec(evaluate(dr,vector_to_matrix(z)))
evaluate(dr::AbstractDecisionRule{CartesianGrid{d1}, CartesianGrid{d2}}, i::Int, y::Union{Vector{Float64},AbstractMatrix{Float64}}) where d1 where d2  = evaluate(dr,node(dr.grid_exo,i), y)


function set_values!(dr::AbstractDecisionRule{<:CartesianGrid, <:CartesianGrid, n_x}, vals::Vector{Matrix{Float64}}) where n_x
    N = size(vals[1],1)
    dims = tuple(dr.grid_endo.n...)
    V = [reinterpret(Value{n_x}, v', dims) for v in vals]
    set_values!(dr, V)
end

# UnstructuredGrid x CartesianGrid


function evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, CartesianGrid{d2},n_x}, i::Int, z::AbstractMatrix{Float64}) where n_x where d1 where d2
    N = size(z,1)
    assert(size(z,2)==d2)
    points = to_LOP(Point{d2}, z)
    out = evaluate(dr,i,points)
    reinterpret(Float64, out, (n_x, N))'
end

evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, CartesianGrid{d2}, n_x}, i::Int, x::Vector{Float64}) where n_x where d1 where d2 = vec(evaluate(dr,i,vector_to_matrix(x)))
function evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, CartesianGrid{d2}, n_x}, inds::Vector{Int}, x::AbstractMatrix{Float64}) where n_x where d1 where d2
    l = [evaluate(dr,inds[n],x[n,:])' for n=1:size(x,1)]
    vcat(l...)
end

function evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, CartesianGrid{d2}, n_x}, inds::Vector{Int}, x::Vector{Float64}) where n_x where d1 where d2
    vcat([evaluate(dr,i,x)' for i in inds]...)
end

function set_values!(dr::AbstractDecisionRule{<:UnstructuredGrid, <:CartesianGrid, n_x}, vals::Vector{Matrix{Float64}}) where n_x
    N = size(vals[1],1)
    dims = tuple(dr.grid_endo.n...)
    V = [to_LOP(Value{n_x}, v) for v in vals]
    set_values!(dr, V)
end
