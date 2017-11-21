@compat abstract type Grid{d} end
#

# # backward backward compatibility
nodes(::Type{<:Matrix}, grid::Grid{d}) where  d  = from_LOP(nodes(grid))

# nodes(::Type{Array{Float64,2}}, grid::Grid) = from_LOP(nodes(grid))

# node(::typeof(Point), grid::Grid{d}, i::Int) where d= SVector{d,Float64}(node(grid, i)...)
# inode(::typeof(Point), grid::Grid{d}, i::Int) where d = SVector{d,Float64}(inode(grid, i)...)



# # # type is not predictible yet
nodes(::Type{<:ListOfPoints}, grid::Grid{d}) where  d = (nodes(grid))
# node(::typeof(Point), grid::Grid{d}, i::Int) where d= SVector{d,Float64}(node(grid, i)...)
# inode(::typeof(Point), grid::Grid{d}, i::Int) where d = SVector{d,Float64}(inode(grid, i)...)
#


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
    @assert d == size(nodes,d)
    UnstructuredGrid{d}(reinterpret(Point{d}, nodes, (N,)))
end

nodes(grid::UnstructuredGrid) = grid.nodes
n_nodes(grid::UnstructuredGrid) = length(grid.nodes)
node(grid::UnstructuredGrid, i::Int) = grid.nodes[i] # fail if i!=1 ?

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



immutable SmolyakGrid{N} <: Grid{N}
    smol_params::BM.SmolyakParams{Float64,Vector{Int}}
    nodes::Matrix{Float64}
    B_nodes::Matrix{Float64}

    function (::Type{SmolyakGrid{N}}){N}(min::Vector{Float64}, max::Vector{Float64}, mu::Vector{Int})
        d = length(min)
        @assert d == N
        dim_err = DimensionMismatch("min was length $d, max and mu must match")
        length(max) == d || throw(dim_err)
        length(mu) == d || throw(dim_err)

        sp = BM.SmolyakParams(d, mu, min, max)
        nodes = BM.nodes(sp)
        B_nodes = BM.evalbase(sp, nodes)
        return new{N}(sp, nodes, B_nodes)
    end
    function (::Type{SmolyakGrid{N}}){N}(min::Array{Float64,1}, max::Array{Float64,1}, mu::Int)
        return SmolyakGrid{N}(min, max, fill(mu, N))
    end
end


immutable RandomGrid{N} <: Grid{N}
    min::Vector{Float64}
    max::Vector{Float64}
    n::Int
    nodes::Matrix{Float64}
    function (::Type{RandomGrid{N}}){N}(min, max, n)
        d = length(min)
        @assert d == N
        dim_err = DimensionMismatch("min was length $d, max must be also")
        length(max) == d || throw(dim_err)
        all(max .> min) || error("max must be greater than min")

        nodes = rand(n, d)  # on [0, 1]
        nodes .= nodes .* (max - min)' .+ min' # scale/shift to [min, max]

        return new{N}(min, max, n, nodes)
    end
end
