abstract Grid

function mlinspace(min, max, n)
    nodes = [linspace(e...) for e in zip(min, max, n)]
    return QE.gridmake(nodes...)
end

type EmptyGrid <: Grid
    # this grid does not exist ;-)
end

nodes(grid::EmptyGrid) = nothing
n_nodes(grid::EmptyGrid) = 0
node(grid::EmptyGrid, i::Int64) = nothing # fail if i!=1 ?


type PointGrid <: Grid
    point::Vector{Float64}
end
nodes(grid::PointGrid) = grid.point[:,:]
n_nodes(grid::PointGrid) = 1
node(grid::PointGrid, i::Int64) = point # fail if i!=1 ?


type UnstructuredGrid <: Grid
    nodes::Matrix{Float64}
end
nodes(grid::UnstructuredGrid) = grid.nodes
n_nodes(grid::UnstructuredGrid) = size(grid.nodes,1)
node(grid::UnstructuredGrid, i::Int64) = grid.nodes[i,:] # fail if i!=1 ?


type CartesianGrid <: Grid
    min::Vector{Float64}
    max::Vector{Float64}
    n::Vector{Int64}
    nodes::Matrix{Float64}
    function CartesianGrid(min, max, n)
        nodes = mlinspace(min, max, n)
        return new(min,max,n,nodes)
    end
end
nodes(grid::Grid) = grid.nodes
n_nodes(grid::Grid) = size(grid.nodes,1)
node(grid::Grid, i::Int64) = grid.nodes[i,:]

type SmolyakGrid <: Grid
    min::Vector{Float64}
    max::Vector{Float64}
    mu::Vector{Int64}
    nodes::Matrix{Float64}
end
