module PlayGround

    import Base: product
    using StaticArrays

    abstract type Grid{d} end

    struct UCGrid{d} <: Grid{d}
        a::SVector{d, Float64}
        b::SVector{d, Float64}
        n::SVector{d, Int64}
        nodes:: Vector{SVector{d, Float64}}
    end

    ndims(g::Grid{d}) where d = d

    import Base
    Base.iterate(ucg::UCGrid{d}, args...) where d = iterate(ucg.nodes, args...)
    Base.length(ucg::UCGrid{d}) where d = length(ucg.nodes)
    Base.eltype(ucg::UCGrid{d}) where d = SVector{d, Float64}

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

        d = length(n)
        @assert (length(a)==d) & (length(b)==d)

        aa = SVector(a...)
        bb = SVector(b...)
        nn = SVector(n...)

        return UCGrid{d}(aa, bb, nn)
    end



    # Concatenate SVectors
end

import Main.PlayGround

grid_1 = Main.PlayGround.UCGrid([0.0, 0.0], [1.0, 1.0], [2, 2])
grid_2 = Main.PlayGround.UCGrid([0.0, 0.0], [1.0, 1.0], [2, 2])

pg = Main.PlayGround.PGrid(grid_1, grid_2)


show(pg)

for k in Main.PlayGround.iter(pg)
    println(k)
end

for k in Main.PlayGround.enumerate_product(pg)
    println(k)
end



