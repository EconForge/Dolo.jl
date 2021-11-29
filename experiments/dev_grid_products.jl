# using Dolo
# import Dolo: Cartesian

module Temp

    using IterTools
    using StaticArrays
    const Point{d} = SVector{d, Float64}

    abstract type AbstractGrid{d} end

    struct CGrid{d} <: AbstractGrid{d}
        min::Point{d}
        max::Point{d}
        n::SVector{d, Int64}
    end

    n_nodes(cg::CGrid{d}) where d = prod(cg.n)

    function iter(cg::CGrid{d}) where d
        # ranges = (range(e[1], e[2]; length=e[3]) for e in zip(cg.min, cg.max, cg.n))
        ranges = (range(cg.min[i], cg.max[i];length=cg.n[i]) for i=1:d)
        return ( SVector{d}(v...) for v in Iterators.product(ranges...) )
    end


    
    CGrid{d}(min::Vector, max::Vector, n::Vector) where d = CGrid(Point{d}(min...), Point{d}(max...), SVector(n...))

    struct PGrid{d, G1, G2} <: AbstractGrid{d}
        g1:: G1
        g2:: G2
    end

    # that one looks efficient
    function iter(pg::PGrid{d, G1, G2}) where d where G1 where G2
        iters = (iter(pg.g1), iter(pg.g2))
        return ( ((k[1][1], k[2][1]), SVector(k[1][2]..., k[2][2]...)) for (k) in Iterators.product(iters...))
    end

    PGrid(g1::AbstractGrid{d1}, g2::AbstractGrid{d2}) where d1 where d2 = PGrid{d1+d2, typeof(g2), typeof(g2)}(g1, g2)
    n_nodes(pg::PGrid) = n_nodes(pg.g1)*n_nodes(pg.g2)

    function test_iter_pg(pg::PGrid)
        k = 0.0
        for i in iter(pg)
            v = i[2]
            k += sum(v)
        end
        return k
    end


    struct PPGrid{T<:Tuple}
        grids::T
    end

    return (
        v
        for v in  Iterators.product(
            tuple( (enumerate(iter(gr)) for gr in ppg.grids)...  )...
        )
    )

    PPGrid(grids...) = PPGrid(grids)

    n_nodes(pg::PPGrid) = prod( (n_nodes(gr) for gr in pg.grids) )

    concat(x1::SVector{d1,Float64}, x2::SVector{d2,Float64}) where d1 where d2 = SVector{d1+d2, Float64}(x1...,x2...)
    concat(x1::SVector{d1,Float64}, x2::SVector{d2,Float64}, x3::SVector{d3,Float64}) where d1 where d2 where d3 = SVector{d1+d2+d3, Float64}(x1...,x2..., x3...)

    function iter(ppg::PPGrid)
        nodes = tuple( ([e for e in iter(gr)] for gr in ppg.grids)... )
        return ( concat(v...) for v in  Iterators.product(nodes...) )
    end
    
    import Base
    function Base.enumerate(ppg::PPGrid)
         return (
            ( tuple( (e[1] for e in v)...) , concat( (e[2] for e in v)...))
            for v in Iterators.product( (enumerate(iter(gr)) for gr in ppg.grids)... )
        )
    end
    
    # Base.getindex(ppg, inds...) = concat()
    
    # uses less memory but slower
    function iter2(ppg::PPGrid)
        return ( concat(v...) for v in Iterators.product( (iter(gr) for gr in ppg.grids) ... ) )
    end



    function test_iter(cg::CGrid{d}) where d
        l = 0.0
        for v in iter(cg)
            l += sum(v)
        end
        return l
    end


    function test_iter_ppg(ppg::PPGrid)
        l = 0.0
        for v in iter(ppg)
            l += sum(v)
        end
        return l
    end


    function test_iter2_ppg(ppg::PPGrid)
        l = 0.0
        for v in iter2(ppg)
            l += sum(v)
        end
        return l
    end


    # PGrid(g1::AbstractGrid{d1}, g2::AbstractGrid{d2}) where d1 where d2 = PGrid{d1+d2, typeof(g2), typeof(g2)}(g1, g2)


end

import Main.Temp as M

g = M.CGrid{2}([0.0, 0.0], [1.0, 1.0], [5,5])

pg = M.PGrid(g,g)


pggg = M.PPGrid(g,g)


@time M.test_iter(g)
@time M.test_iter_pg(pg)
@time M.test_iter_ppg(pggg)
@time M.test_iter2_ppg(pggg)


@time M.test_iter(g)
@time M.test_iter_pg(pg)
@time M.test_iter_ppg(pggg)
@time M.test_iter2_ppg(pggg)


@code_warntype M.test_iter(g)

@code_warntype M.test_iter_ppg(pggg)
@code_warntype M.test_iter2_ppg(pggg)



function test()
    t = 0.0
    for k in M.enumerate(pggg)
        t += sum(k[2])
        
    end
    return t
end