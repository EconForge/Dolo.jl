module PlayGround

    import Base: product
    using StaticArrays

    abstract type Grid{n_x} end

    struct UCGrid{n_x} <: Grid{n_x}
        a::SVector{n_x, Float64}
        b::SVector{n_x, Float64}
        n::SVector{n_x, Int64}
        nodes:: Vector{SVector{n_x, Float64}}
    end

    import Base
    Base.iterate(ucg::UCGrid{n_x}, args...) where n_x = iterate(ucg.nodes, args...)
    Base.length(ucg::UCGrid{n_x}) where n_x = length(ucg.nodes)
    Base.eltype(ucg::UCGrid{n_x}) where n_x = SVector{n_x, Float64}

    function UCGrid{d}(a::SVector{d, Float64},
                         b::SVector{d, Float64},
                         n::SVector{d, Int64}
    ) where d
        intervals = tuple( [range(a[i], b[i]; length=n[i]) for i=1:d]... )
        # print( [product(intervals)...])
        # nn = [(e) for e in product(intervals...)]
        # return nn
        cnodes = [SVector(e...) for e in product(intervals...)]
        # TODO: can we avoid the cnodes[:] allocation
        return UCGrid{d}(a,b,n,cnodes[:])
    end

    function UCGrid(
        a::AbstractVector{Float64},
        b::AbstractVector{Float64},
        n::AbstractVector{Int64})

        n_x = length(n)
        @assert (length(a)==n_x) & (length(b)==n_x)

        aa = SVector(a...)
        bb = SVector(b...)
        nn = SVector(n...)

        return UCGrid{n_x}(aa, bb, nn)
    end

end

import Main.PlayGround

grid_1 = Main.PlayGround.UCGrid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 2, 2])
grid_2 = Main.PlayGround.UCGrid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 2, 2])

for (k,v) in IterTools.product(grid_1, grid_2)

    print(k,v)

end