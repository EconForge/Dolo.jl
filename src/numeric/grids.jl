@compat abstract type Grid end

function Base.show(io::IO, grid::Grid)
    print(typeof(grid))
end

function mlinspace(min, max, n)
    nodes = map(linspace, min, max, n)
    return QE.gridmake(nodes...)
end

immutable EmptyGrid <: Grid
    # this grid does not exist ;-)
end

nodes(grid::EmptyGrid) = nothing
n_nodes(grid::EmptyGrid) = 0
node(grid::EmptyGrid, i::Int) = nothing # fail if i!=1 ?


immutable PointGrid <: Grid
    point::Vector{Float64}
end
nodes(grid::PointGrid) = grid.point'
n_nodes(grid::PointGrid) = 1
node(grid::PointGrid, i::Int) = point # fail if i!=1 ?


immutable UnstructuredGrid <: Grid
    nodes::Matrix{Float64}
end
nodes(grid::UnstructuredGrid) = grid.nodes
n_nodes(grid::UnstructuredGrid) = size(grid.nodes, 1)
node(grid::UnstructuredGrid, i::Int) = grid.nodes[i, :] # fail if i!=1 ?


immutable CartesianGrid <: Grid
    min::Vector{Float64}
    max::Vector{Float64}
    n::Vector{Int}
    nodes::Matrix{Float64}
    function CartesianGrid(min, max, n)
        nodes = mlinspace(min, max, n)
        return new(min, max, n, nodes)
    end
end
nodes(grid::Grid) = grid.nodes
n_nodes(grid::Grid) = size(grid.nodes, 1)
node(grid::Grid, i::Int) = grid.nodes[i, :]

immutable SmolyakGrid <: Grid
    smol_params::BM.SmolyakParams{Float64,Vector{Int}}
    nodes::Matrix{Float64}
    B_nodes::Matrix{Float64}

    function SmolyakGrid(min::Vector{Float64}, max::Vector{Float64}, mu::Vector{Int})
        d = length(min)
        dim_err = DimensionMismatch("min was length $d, max and mu must match")
        length(max) == d || throw(dim_err)
        length(mu) == d || throw(dim_err)

        sp = BM.SmolyakParams(d, mu, min, max)
        nodes = BM.nodes(sp)
        B_nodes = BM.evalbase(sp, nodes)
        return new(sp, nodes, B_nodes)
    end
end

function SmolyakGrid(min::Array{Float64,1}, max::Array{Float64,1}, mu::Int)
    return SmolyakGrid(min, max, fill(mu, length(min)))
end
