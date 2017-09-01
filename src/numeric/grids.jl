@compat abstract type Grid{N} end

function Base.show(io::IO, grid::Grid)
    print(io, typeof(grid))
end

function mlinspace(min, max, n)
    nodes = map(linspace, min, max, n)
    return QE.gridmake(nodes...)
end

immutable EmptyGrid{N} <: Grid{N}
    # this grid does not exist ;-)
end

nodes(grid::EmptyGrid) = nothing
n_nodes(grid::EmptyGrid) = 0
node(grid::EmptyGrid, i::Int) = nothing # fail if i!=1 ?


immutable PointGrid{N} <: Grid{N}
    point::Vector{Float64}

    function (::Type{PointGrid{N}}){N}(point::Vector{Float64})
        @assert N == length(point)
        new{N}(point)
    end
end
nodes(grid::PointGrid) = grid.point'
n_nodes(grid::PointGrid) = 1
node(grid::PointGrid, i::Int) = point # fail if i!=1 ?


immutable UnstructuredGrid{N} <: Grid{N}
    nodes::Matrix{Float64}

    function (::Type{UnstructuredGrid{N}}){N}(nodes::Matrix{Float64})
        @assert N == size(nodes, 2)
        new{N}(nodes)
    end
end
nodes(grid::UnstructuredGrid) = grid.nodes
n_nodes(grid::UnstructuredGrid) = size(grid.nodes, 1)
node(grid::UnstructuredGrid, i::Int) = grid.nodes[i, :] # fail if i!=1 ?


immutable CartesianGrid{N} <: Grid{N}
    min::Vector{Float64}
    max::Vector{Float64}
    n::Vector{Int}
    nodes::Matrix{Float64}
    function (::Type{CartesianGrid{N}}){N}(min, max, n)
        @assert N == length(min) == length(max) == length(n)
        nodes = mlinspace(min, max, n)
        return new{N}(min, max, n, nodes)
    end
end
nodes(grid::Grid) = grid.nodes
n_nodes(grid::Grid) = size(grid.nodes, 1)
node(grid::Grid, i::Int) = grid.nodes[i, :]

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
