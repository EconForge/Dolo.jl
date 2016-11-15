abstract Grid

function mlinspace(min, max, n)
    nodes = [linspace(e...) for e in zip(min, max, n)]
    return QE.gridmake(nodes...)
end

type CartesianGrid <: Grid
    min::Array{Float64,1}
    max::Array{Float64,1}
    n::Array{Int64,1}
end


type SmolyakGrid <: Grid
    min::Array{Real,1}
    max::Array{Real,1}
    mu::Array{Int,1}
end

function nodes(grid::CartesianGrid)
    return mlinspace(grid.min, grid.max, grid.n)
end
