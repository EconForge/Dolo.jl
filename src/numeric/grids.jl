abstract type AbstractGrid end

abstract type Grid{d} end

struct ProductGrid{T,S}
    exo::T
    endo::S
end



×(grid1, grid2) = ProductGrid(grid1, grid2)
#

# # backward backward compatibility
nodes(::Type{<:Union{ListOfPoints,ListOfPoints{d}}}, grid::Grid{d}) where d = nodes(grid)
nodes(::Type{<:Matrix}, grid::Grid) = copy(from_LOP(nodes(grid)))

node(::Type{<:Union{Point,Point{d}}}, grid::Grid{d}, i::Int) where d = node(grid,i)
node(::Type{<:Vector}, grid::Grid, i::Int) = Vector(node(grid,i))

import Base


struct PGrid{T1, T2}
    exo::T1
    endo::T2
end

Base.ndims(grid::Grid{d}) where d = d

function Base.show(io::IO, grid::Grid)
    print(io, typeof(grid))
end


struct EmptyGrid <: Grid{0}
    # this grid does not exist ;-)
end

nodes(grid::EmptyGrid) = nothing
n_nodes(grid::EmptyGrid) = 0
node(grid::EmptyGrid, i::Int) = nothing # fail if i!=1 ?

##########################
# Grid made of one point #
##########################

struct PointGrid{d} <: Grid{d}
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

struct UnstructuredGrid{d} <: Grid{d}
    nodes::ListOfPoints{d}
end

# Old-convention
function UnstructuredGrid{d}(nodes::Matrix{Float64}) where d
    N = size(nodes,1)
    @assert d == size(nodes,2)
    UnstructuredGrid{d}(reshape(reinterpret(Point{d}, vec(copy(nodes'))), (N,)))
end

nodes(grid::UnstructuredGrid) = grid.nodes
n_nodes(grid::UnstructuredGrid) = length(grid.nodes)
node(grid::UnstructuredGrid, i::Int) = grid.nodes[i] # fail if i!=1 ?
node(::Type{<:Point}, grid::UnstructuredGrid, i::Int) = node(grid,i)

function Product(a::UnstructuredGrid{d1}, b::UnstructuredGrid{d2}) where d1 where d2
    A = [Base.product(a.nodes, b.nodes)...]
    N = length(A)
    d = d1 + d2
    nodes = reshape(reinterpret(Point{d},(A)),(N,))
    return UnstructuredGrid{d}(nodes)
end

#################
# UCGrid #
#################

function mlinspace(min, max, n)
    # this now returns an iterator
    nodes = map((x,y,z)->range(x,stop=y,length=z), min, max, n)
    return Base.product(nodes...)
end

struct UCGrid{d} <: Grid{d}
    min::Point{d}
    max::Point{d}
    n::SVector{d,Int}
    nodes::ListOfPoints{d}
end

ndims(grid::UCGrid{d}) where d = d

function (::Type{<:UCGrid})(min::SVector{d,Float64}, max::SVector{d,Float64}, n::SVector{d,Int64}) where d
    A = [mlinspace(min, max, n)...]
    N = prod(n)
    mm = reshape(reinterpret(Point{d},vec(A)),(N,))
    return UCGrid{d}(min, max, n, mm)
end

(::Type{<:UCGrid})(min::Vector{Float64},max::Vector{Float64},n::Vector{Int64}) = UCGrid(SVector(min...), SVector(max...), SVector(n...))
(::Type{<:UCGrid})(min::Vector{Float64},max::Vector{Float64},n::Int64) = UCGrid(min, max, fill(n,length(min)))

scales(grid::UCGrid{d}) where d = Tuple{Vararg{Vector{Float64},d}}([range(grid.min[i], grid.max[i];length=grid.n[i]) for i=1:d])

nodes(grid::UCGrid{d}) where d = grid.nodes
n_nodes(grid::UCGrid{d}) where d = length(grid.nodes)
node(grid::UCGrid{d}, i::Int) where d = grid.nodes[i]

nodes(::Type{<:ListOfPoints}, grid::UCGrid{d}) where d  = grid.nodes
node(::Type{<:Point}, grid::UCGrid{d}, i::Int64) where d = node(grid,i)

function Product(a::UCGrid{d1}, b::UCGrid{d2}) where d1 where d2
  return Dolo.UCGrid{d1+d2}( [a.min; b.min], [a.max; b.max], [a.n; b.n])
end


################
# Smolyak Grid #
################

struct SmolyakGrid{d} <: Grid{d}
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
    mmu = mu isa Int ? ffill(mu,d) : mu
    @assert d == length(max) == length(mmu)
    SmolyakGrid{d}(SVector(min...), SVector(max...), SVector(mu...))
end
#
function SmolyakGrid{d}(min::Vector{Float64},max::Vector{Float64},mu::Union{Vector{Int64},Int64}) where d
    mmu = mu isa Int ? fill(mu,d) : mu
    @assert d == length(min) == length(max) == length(mmu)
    SmolyakGrid{d}(SVector(min...), SVector(max...), SVector(mu...))
end


nodes(grid::SmolyakGrid{d}) where d  = reshape(reinterpret(Point{d}, vec(copy(grid.nodes'))),( size(grid.nodes,1),))
n_nodes(grid::SmolyakGrid) = size(grid.nodes,1)
node(grid::SmolyakGrid{d}, i::Int) where d = Point{d}(grid.nodes[i,:]...)

###############
# Random Grid #
###############

struct RandomGrid{d} <: Grid{d}
    min::SVector{d,Float64}
    max::SVector{d,Float64}
    n::Int64
    nodes::ListOfPoints{d}
end

function  (::Type{<:Union{RandomGrid,RandomGrid{d}}})(min::Point{d}, max::Point{d}, n::Int) where d
    nodes = reshape(reinterpret(Point{d}, vec(rand(d, n))), (n,))  # on [0, 1]
    for nn=1:length(n)
        nodes[nn] = nodes[nn] .* (max-min) + min
    end
    RandomGrid(min, max, n, copy(nodes))
end

# function RandomGrid(min::Vector{Float64}, max::Vector{Float64}, n::Int)
#     d = length(min)
#     RandomGrid{d}(SVector(min...),SVector(max...),n)
# end

function RandomGrid{d}(min::Vector{Float64}, max::Vector{Float64}, n::Int) where d
    @assert d == length(min)
    all(max .> min) || error("max must be greater than min")
    RandomGrid{d}(SVector{d}(min...),SVector{d}(max...),n)
end

nodes(grid::RandomGrid) = grid.nodes
nodes(::Type{<:ListOfPoints}, grid::RandomGrid) = nodes(grid)

n_nodes(grid::RandomGrid) = length(grid.nodes)

node(grid::RandomGrid,i::Int) = grid.nodes[i]


### Cartesian

struct Cartesian <: AbstractGrid
    a::Vector{Float64}
    b::Vector{Float64}
    orders::Vector{Int}
end


### Domains




abstract type AbstractDomain end


struct EmptyDomain <: AbstractDomain 
    states::Vector{Symbol}
end



struct CartesianDomain<: AbstractDomain
    states
    min::Vector{Float64}
    max::Vector{Float64}
end

const Domain = CartesianDomain

ndims(dom::CartesianDomain) = length(dom.min)


function discretize(dom::CartesianDomain; n=Union{Int, Vector{Int}})
    if typeof(n)<:Int
        nv = fill(n, ndims(dom))
    else
        nv = n
    end
    min = dom.min
    max = dom.max
    return CartesianGrid(min, max, n)
    
end
