abstract type AbstractGrid end

abstract type Grid{d} end

import Base: ndims, length

Base.length(g::Grid{d}) where d = max(1,n_nodes(g))  # we should deprecate n_nodes

Base.ndims(grid::T) where T<:Grid{d} where d = d

function Base.show(io::IO, grid::Grid)
    print(io, typeof(grid))
end


# # # backward backward compatibility
nodes(::Type{<:Union{ListOfPoints,ListOfPoints{d}}}, grid::Grid{d}) where d = nodes(grid)
nodes(::Type{<:Matrix}, grid::Grid) = copy(from_LOP(nodes(grid)))

node(::Type{<:Union{Point,Point{d}}}, grid::Grid{d}, i::Int) where d = node(grid,i)
node(::Type{<:Vector}, grid::Grid, i::Int) = Vector(node(grid,i)...)


import Base

### Product Grids

# special product grid, with endo/exo distinction

struct ProductGrid{T,S}
    exo::T
    endo::S
end

⊗(grid1::Grid{d1}, grid2::Grid{d2}) where d1 where d2 = ProductGrid(grid1, grid2)

Base.show(io::IO, pg::ProductGrid) = print(io, pg.exo, "⊗",  pg.endo)
Base.length(g::ProductGrid{g1,g2}) where g1 where g2 = length(g.exo)*length(g.endo)

struct PGrid{d, T<:Tuple} <: Grid{d}
    grids::T
end

function PGrid(grids...)
    d = sum( ndims(g) for g in grids)
    t = typeof(grids)
    return PGrid{d, t}(grids)
end

function iter(pg::T) where T <:PGrid
    (concat(e...) for e in Iterators.product( (g for g in pg.grids)... ))
end

Base.length(pg:: PGrid{d, T} ) where d where T = prod( (length(g) for g in pg.grids) )
Base.eltype(pg:: PGrid{d, T} ) where d where T = SVector{d, Float64}

function Base.show(io::IO,x::PGrid{d,T})  where d where T
    N = length(x.grids)
    print(io, x.grids[1])
    if N==1
        return
    end
    for n=2:N
        print(" ⨯ ")
        print(io, x.grids[n])
    end
end

function enumerate_product(pg::T) where T <:PGrid
    iter1 = Iterators.product( (1:length(g) for g in pg.grids)... )
    iter2 = (concat(e...) for e in Iterators.product( (g for g in pg.grids)... ))
    # return iter1, iter2
    zip(iter1, iter2)
end



### EmptyGrid


## TODO: this should probably be of dimension d
struct EmptyGrid{d} <: Grid{d}
    # this grid does not exist ;-)
end
Base.show(io::IO, g::EmptyGrid{d}) where d = print(io, "EmptyGrid{$d}")

nodes(grid::EmptyGrid) = nothing
n_nodes(grid::EmptyGrid) = 0 ##### Reconsider: ???
function node( grid::EmptyGrid{d}, i::Int) where d
    return fill(NaN, SVector{d, Float64})
end




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

Base.show(io::IO, gr::UnstructuredGrid{d}) where d = print(io, "UNSgrid{$d}")

# Old-convention
function UnstructuredGrid{d}(nodes::Matrix{Float64}) where d
    N = size(nodes,1)
    @assert d == size(nodes,2)
    UnstructuredGrid{d}(reshape(reinterpret(Point{d}, vec(copy(nodes'))), (N,)))
end

nodes(grid::UnstructuredGrid) = grid.nodes
n_nodes(grid::UnstructuredGrid) = length(grid.nodes)
node(grid::UnstructuredGrid, i::Int) = grid.nodes[i] # fail if i!=1 ?

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

Base.iterate(ucg::UCGrid{d}, args...) where d = iterate(ucg.nodes, args...)
Base.length(ucg::UCGrid{d}) where d = length(ucg.nodes)
Base.eltype(ucg::UCGrid{d}) where d = SVector{d, Float64}


function Base.show(io::IO, g::UCGrid{d}) where d
    print("UCGrid{", d, "}")
end


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
  return Dolo.UCGrid{d1+d2}(
      SVector(a.min..., b.min...),
      SVector(a.max..., b.max...),
      SVector(a.n..., b.n...)
  )
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


concat(v...) = length(v)==1 ? v[1] : SVector(v[1]..., concat(v[2:end]...)... )
