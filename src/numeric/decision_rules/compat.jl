## Common:


### backward compat

function DecisionRule(exo_grid, endo_grid, values::Vector{Matrix{Float64}})
    N,n_x = size(values[1])
    vals = [to_LOP(v) for v in values]
    DecisionRule(exo_grid, endo_grid, vals)
end

function DecisionRule(exo_grid, endo_grid, values::Vector{Matrix{Float64}}, tt)
    N,n_x = size(values[1])
    vals = [reshape(reinterpret(Value{n_x}, vec(v')), (N,)) for v in values]
    DecisionRule(exo_grid, endo_grid, vals, tt)
end

# EmptyGrid x UCGrid

function evaluate(dr::AbstractDecisionRule{EmptyGrid{d1}, UCGrid{d}, n_x}, z::AbstractMatrix{Float64}) where n_x where d
    N = size(z,1)
    @assert size(z,2)==d
    points = to_LOP(Point{d}, z)
    out = evaluate(dr,copy(points))
    return reshape(reinterpret(Float64, vec(out)), (n_x,N))'
end

evaluate(dr::AbstractDecisionRule{EmptyGrid{d1}, <:UCGrid}, z::Vector{Float64})  where d1 = vec(evaluate(dr,z'))
evaluate(dr::AbstractDecisionRule{EmptyGrid{d1}, <:UCGrid}, i::Int, y::Union{Vector,AbstractMatrix})  where d1 = evaluate(dr,y) # kind of nonsensical
evaluate(dr::AbstractDecisionRule{EmptyGrid{d1}, <:UCGrid}, x::Union{Vector,AbstractMatrix}, y::Union{Vector,AbstractMatrix})  where d1 = evaluate(dr,y)
UCGrid
function set_values!(dr::AbstractDecisionRule{EmptyGrid{d1}, <:UCGrid, n_x}, vals::Vector{Matrix{Float64}}) where n_x where d1
    N = size(vals[1],1)
    dims = tuple(dr.grid_endo.n...)
    V = [reshape(reinterpret(Value{n_x}, vec(vals[1]')), dims)]
    set_values!(dr, V)
end

# UCGrid x UCGrid

function evaluate(dr::AbstractDecisionRule{UCGrid{d1}, UCGrid{d2}, n_x}, x::AbstractMatrix, y::AbstractMatrix) where d1 where d2 where n_x
    N = size(x,1)
    @assert size(y,1)==N
    xx = to_LOP(Point{d1}, x)
    yy = to_LOP(Point{d2}, y)
    res = evaluate(dr,copy(xx),copy(yy))
    return reshape(reinterpret(Float64,vec(res)), (n_x, N))'
end
#
evaluate(dr::AbstractDecisionRule{<:UCGrid, <:UCGrid}, x::Vector{Float64}, y::Matrix{Float64}) = evaluate(dr, repeat(x', size(y,1)), y)
evaluate(dr::AbstractDecisionRule{<:UCGrid, <:UCGrid}, x::Vector{Float64}, y::Vector{Float64}) = evaluate(dr, vector_to_matrix(x), vector_to_matrix(y))
evaluate(dr::AbstractDecisionRule{UCGrid{d1}, UCGrid{d2}}, z::AbstractMatrix) where d1 where d2 = evaluate(dr, [z[:,1:d1] z[:,1:(d1+1:end)]])
evaluate(dr::AbstractDecisionRule{<:UCGrid, <:UCGrid}, z::Vector{Float64}) = vec(evaluate(dr,vector_to_matrix(z)))
evaluate(dr::AbstractDecisionRule{UCGrid{d1}, UCGrid{d2}}, i::Int, y::Union{Vector{Float64},AbstractMatrix{Float64}}) where d1 where d2  = evaluate(dr,node(dr.grid_exo,i), y)


function set_values!(dr::AbstractDecisionRule{<:UCGrid, <:UCGrid, n_x}, vals::Vector{Matrix{Float64}}) where n_x
    N = size(vals[1],1)
    dims = tuple(dr.grid_endo.n...)
    V = [reshape(reinterpret(Value{n_x}, vec(v')), dims) for v in vals]
    set_values!(dr, V)
end

# UnstructuredGrid x UCGrid


function evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, UCGrid{d2},n_x}, i::Int, z::AbstractMatrix{Float64}) where n_x where d1 where d2
    N = size(z,1)
    @assert size(z,2)==d2
    points = to_LOP(Point{d2}, z)
    # TODO: remove copy()
    out = evaluate(dr,i,copy(points))
    reshape(reinterpret(Float64, vec(out)), (n_x, N))'
end

evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, UCGrid{d2}, n_x}, i::Int, x::Vector{Float64}) where n_x where d1 where d2 = vec(evaluate(dr,i,vector_to_matrix(x)))
function evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, UCGrid{d2}, n_x}, inds::Vector{Int}, x::AbstractMatrix{Float64}) where n_x where d1 where d2
    l = [evaluate(dr,inds[n],x[n,:])' for n=1:size(x,1)]
    vcat(l...)
end

function evaluate(dr::AbstractDecisionRule{UnstructuredGrid{d1}, UCGrid{d2}, n_x}, inds::Vector{Int}, x::Vector{Float64}) where n_x where d1 where d2
    vcat([evaluate(dr,i,x)' for i in inds]...)
end

function set_values!(dr::AbstractDecisionRule{<:UnstructuredGrid, <:UCGrid, n_x}, vals::Vector{Matrix{Float64}}) where n_x
    N = size(vals[1],1)
    dims = tuple(dr.grid_endo.n...)
    V = [to_LOP(Value{n_x}, v) for v in vals]
    set_values!(dr, V)
end
