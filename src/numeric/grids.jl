@compat abstract type Grid{d} end
#

# # backward backward compatibility
nodes(::Type{<:Union{ListOfPoints,ListOfPoints{d}}}, grid::Grid{d}) where d = nodes(grid)
nodes(::Type{<:Matrix}, grid::Grid) = from_LOP(nodes(grid))

node(::Type{<:Union{Point,Point{d}}}, grid::Grid{d}, i::Int) where d = node(grid,i)
node(::Type{<:Vector}, grid::Grid, i::Int) = Vector(node(grid,i))



import Base

Base.ndims(grid::Grid{d}) where d = d

function Base.show(io::IO, grid::Grid)
    print(io, typeof(grid))
end


immutable EmptyGrid <: Grid{0}
    # this grid does not exist ;-)
end

nodes(grid::EmptyGrid) = nothing
n_nodes(grid::EmptyGrid) = 0
node(grid::EmptyGrid, i::Int) = nothing # fail if i!=1 ?

##########################
# Grid made of one point #
##########################

immutable PointGrid{d} <: Grid{d}
    point::SVector{d,Float64}
end

function (::Type{<:PointGrid})(point::Vector{Float64})
    d = length(point)
    PointGrid{d}(SVector{d,Float64}(point...))
end

nodes(grid::PointGrid) = [grid.point]
n_nodes(grid::PointGrid) = 1
node(grid::PointGrid, i::Int) = grid.point # fail if i!=1 ?


#####################
# Unstructured Grid #
#####################

immutable UnstructuredGrid{d} <: Grid{d}
    nodes::ListOfPoints{d}
end

# Old-convention
function (::Type{UnstructuredGrid{d}})(nodes::Matrix{Float64}) where d
    N = size(nodes,1)
    @assert d == size(nodes,2)
    UnstructuredGrid{d}(reinterpret(Point{d}, nodes, (N,)))
end

nodes(grid::UnstructuredGrid) = grid.nodes
n_nodes(grid::UnstructuredGrid) = length(grid.nodes)
node(grid::UnstructuredGrid, i::Int) = grid.nodes[i] # fail if i!=1 ?
node(::Type{<:Point}, grid::UnstructuredGrid, i::Int) = node(grid,i)

function Product(a::UnstructuredGrid{d1}, b::UnstructuredGrid{d2}) where d1 where d2
    A = [Base.product(a.nodes, b.nodes)...]
    N = length(A)
    d = d1 + d2
    nodes = reinterpret(Point{d},A,(N,))
    return UnstructuredGrid{d}(nodes)
end

#################
# CartesianGrid #
#################

function mlinspace(min, max, n)
    # this now returns an iterator
    nodes = map(linspace, min, max, n)
    return Base.product(nodes...)
end

immutable CartesianGrid{d} <: Grid{d}
    min::Point{d}
    max::Point{d}
    n::SVector{d,Int}
    nodes::ListOfPoints{d}
end

function (::Type{<:CartesianGrid})(min::SVector{d,Float64}, max::SVector{d,Float64}, n::SVector{d,Int64}) where d
    A = [mlinspace(min, max, n)...]
    N = prod(n)
    mm = reinterpret(Point{d},A,(N,))
    return CartesianGrid{d}(min, max, n, mm)
end

(::Type{<:CartesianGrid})(min::Vector{Float64},max::Vector{Float64},n::Vector{Int64}) = CartesianGrid(SVector(min...), SVector(max...), SVector(n...))

nodes(grid::CartesianGrid{d}) where d = grid.nodes
n_nodes(grid::CartesianGrid{d}) where d = length(grid.nodes)
node(grid::CartesianGrid{d}, i::Int) where d = grid.nodes[i]

nodes(::Type{<:ListOfPoints}, grid::CartesianGrid{d}) where d  = grid.nodes
node(::Type{<:Point}, grid::CartesianGrid{d}, i::Int64) where d = node(grid,i)

function Product(a::CartesianGrid{d1}, b::CartesianGrid{d2}) where d1 where d2
  return Dolo.CartesianGrid{d1+d2}( [a.min; b.min], [a.max; b.max], [a.n; b.n])
end


################
# Smolyak Grid #
################

immutable SmolyakGrid{d} <: Grid{d}
    smol_params::BM.SmolyakParams{Float64,Vector{Int}}
    nodes::Matrix{Float64}
    B_nodes::Matrix{Float64}
end

function SmolyakGrid{d}(min::Point{d}, max::Point{d}, mu::SVector{d,Int64}) where d
    sp = BM.SmolyakParams(d, Vector(mu), Vector(min), Vector(max))
    nodes = BM.nodes(sp)
    B_nodes = BM.evalbase(sp, nodes)
    return SmolyakGrid{d}(sp, nodes, B_nodes)
end

SmolyakGrid(min::Point{d}, max::Point{d}, mu::SVector{d,Int64}) where d = SmolyakGrid{d}(min,max,mu)
SmolyakGrid(min::Point{d}, max::Point{d}, mu::Int64) where d = SmolyakGrid{d}(min, max, SVector([mu for i=1:d]...))
SmolyakGrid{d}(min::Point{d}, max::Point{d}, mu::Int64) where d = SmolyakGrid{d}(min, max, SVector([mu for i=1:d]...))

function SmolyakGrid(min::Vector{Float64},max::Vector{Float64},mu::Union{Vector{Int64},Int64})
    d = length(min)
    mmu = mu isa Int? ffill(mu,d) : mu
    @assert d == length(max) == length(mmu)
    SmolyakGrid{d}(SVector(min...), SVector(max...), SVector(mu...))
end
#
function SmolyakGrid{d}(min::Vector{Float64},max::Vector{Float64},mu::Union{Vector{Int64},Int64}) where d
    mmu = mu isa Int? fill(mu,d) : mu
    @assert d == length(min) == length(max) == length(mmu)
    SmolyakGrid{d}(SVector(min...), SVector(max...), SVector(mu...))
end


nodes(grid::SmolyakGrid{d}) where d  = reinterpret(Point{d}, grid.nodes',( size(grid.nodes,1),))
n_nodes(grid::SmolyakGrid) = size(grid.nodes,1)
node(grid::SmolyakGrid{d}, i::Int) where d = Point{d}(grid.nodes[i,:]...)

###############
# Random Grid #
###############

immutable RandomGrid{d} <: Grid{d}
    min::SVector{d,Float64}
    max::SVector{d,Float64}
    n::Int64
    nodes::ListOfPoints{d}
end

function  (::Type{<:Union{RandomGrid,RandomGrid{d}}})(min::Point{d}, max::Point{d}, n::Int) where d
    nodes = reinterpret(Point{d}, rand(d, n), (n,))  # on [0, 1]
    for n=1:length(nodes)
        nodes[n] = nodes[n] .* (max-min) + min
    end
    RandomGrid(min, max, n, nodes)
end

function (::Type{RandomGrid})(min::Vector{Float64}, max::Vector{Float64}, n::Int)
    d = length(min)
    @assert d == length(max)
    dim_err = DimensionMismatch("min was length $d, max must be also")
    length(max) == d || dim_err
    all(max .> min) || error("max must be greater than min")
    RandomGrid{d}(SVector(min...),SVector(max...),n)
end

function (::Type{RandomGrid{d}})(min::Vector{Float64}, max::Vector{Float64}, n::Int) where d
    @assert d == length(min)
    RandomGrid(min,max,d)
end

nodes(grid::RandomGrid) = grid.nodes
nodes(::Type{<:ListOfPoints}, grid::RandomGrid) = nodes(grid)

n_nodes(grid::RandomGrid) = length(grid.nodes)

node(grid::RandomGrid,i::Int) = grid.nodes[i]
