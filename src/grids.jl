abstract type AGrid{d} end
abstract type ASGrid{d} <: AGrid{d} end

import Base: eltype, iterate, size

eltype(cg::AGrid{d}) where d = SVector{d, Float64}
ndims(cg::AGrid{d}) where d = d

struct CGrid{d} <: AGrid{d}
    ranges::NTuple{d, Tuple{Float64, Float64, Int64}}
end

size(cg::CGrid{d}) where d = tuple((e[3] for e in cg.ranges)... )


getindex(g::CGrid{d}, inds::Vararg{Int64,d}) where d = SVector{d}(
    ( 
        ( g.ranges[i][1] + (g.ranges[i][2]-g.ranges[i][1])*( (inds[i]-1)/(g.ranges[i][3]-1)) )
        for i=1:d
    )
)

getindex(g::CGrid{d}, ci::CartesianIndex) where d = g[ci.I...]

# enum(g::CGrid{1})= (QP((i,), g[i]) for i in 1:length(g))
# enum(g::CGrid{1})= ((;loc=(i,), val=g[i]) for i in 1:length(g))
enum(g::CGrid{d}) where d = (QP(c, g[c...]) for c in Iterators.product( tuple( ((1:r[3]) for r in g.ranges)... )... ) )

getindex(g::CGrid{1}, ::Colon) = [SVector(i) for i in range(g.ranges[1]...)]

@inline to__linear_index(g::CGrid{d}, ind::Vararg{Int64, d}) where d = LinearIndices(
    tuple( (1:r[3] for r in g.ranges)... )
)[ind...]

@inline from_linear(g::CGrid{d}, n) where d = CartesianIndices(
    tuple( (1:r[3] for r in g.ranges)... )
)[n].I




struct SGrid{n,d} <: ASGrid{d}
    points::SVector{n,SVector{d, Float64}}
end

function SGrid(Q::Matrix)
    n,d = size(Q)
    return SGrid{n,d}([SVector(Q[i,:]...) for i=1:size(Q,1)])
end

function SGrid(v::Vector)
    return SGrid(SVector(v...))
end


struct ProductGrid{G1, G2, d} <: AGrid{d}
    g1::G1
    g2::G2
    # points::Vector{SVector{d, Float64}}
end

const PGrid = ProductGrid


getindex(g::SGrid{d}, ::Colon) where d = g.points
getindex(g::SGrid{d}, i::Int) where d = g.points[i]
getindex(g::PGrid, c::CartesianIndex) = g[c[1],c[2]]
getindex(g::PGrid, ::Colon) = [g...]

cover(m,v::SVector{d,T}) where d where T = SVector{d,T}(
    m...,
    (v[i] for i=length(m)+1:length(v))...
)

PGrid(g1::SGrid{n, d1}, g2::CGrid{d2}) where n where d1 where d2 = PGrid{typeof(g1), CGrid{d2}, d1+d2}(g1, g2)
cross(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid(g1,g2)

# another way to define multi-dimension cartesian grids
PGrid(g1::CGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid{typeof(g1), CGrid{d2}, d1+d2}(g1, g2)

import Base: getindex

from_linear(g::PGrid{G1, G2, d}, n) where G1 where G2 where d = let x=divrem(n-1, length(g.g1)); (x[2]+1, x[1]+1) end

getindex(g::PGrid{G1, G2, d}, n::Int) where G1 where G2 where d = getindex(g, from_linear(g, n)...)

function getindex(g::PGrid{G1, G2, d}, i::Int64, j::Int64) where G1<:SGrid{d1} where G2<:CGrid{d2} where d where d1 where d2
    SVector{d,Float64}(g.g1[i]..., g.g2[j]...)
end


getindex(g::PGrid{G1, G2, d}, i::Int64, ::Colon) where G1 where G2 where d = g.g2[:] # TODO: should error if i out of bounds
getindex(g::PGrid{G1, G2, d}, ::Colon, i::Int64) where G1 where G2 where d = g.g1[:]

@inline to__linear_index(g::CGrid{2}, ind::Tuple{Int64, Int64}) = let 
    i,j = ind
    p = g.ranges[1][3]
    return i + p*(j-1)
end
   
@inline to__linear_index(g::PGrid, ind::Tuple{Int64, Int64}) =  ind[1] + length(g.g1)*(ind[2]-1)


show(io::IO, g::SGrid{d1, d2}) where d1 where d2 = print(io, "SGrid{$(d1)}")

show(io::IO, g::CGrid{d}) where d = let 
    s = join( tuple((r[3] for r in g.ranges)...), "×")
    print(io,"CGrid{$( s )}")
end


show(io::IO, g::PGrid) = println(io, "$(g.g1)×$(g.g2)")

import Base: iterate
import Base: length
import Base: getindex
import Base: setindex!


length(pg::PGrid{G1, G2, d}) where G1 where G2 where d = length(pg.g1)*length(pg.g2)
length(sg::SGrid{d}) where d = length(sg.points)
length(cg::CGrid{d}) where d = prod(e[3] for e in cg.ranges)

### Fast iterators for 1d case

import Base: iterate

getindex(s::CGrid{d}, n::Int) where d = let
    t = tuple( (e[3] for e in s.ranges)... )
    inds = CartesianIndices( t )[n].I
    s[inds...]
end

iterate(s::SGrid) = (s.points[1], 2)
iterate(s::SGrid, state) = state<=length(s) ? (s.points[state], state+1) : nothing

eltype(cg::CGrid{d}) where d = SVector{d, Float64}

function iterate(cg::CGrid{d}) where d
    
    tt = tuple( (1:(e[3]) for e in cg.ranges)...)
    # pit = Iterators.ProductIterator{Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}(tt)
    pit = Iterators.ProductIterator(tt)
    c1, state = iterate(pit)

    return cg[c1...], (pit,state)

end

function iterate(cg::CGrid{d}, sstate) where d
    
    pit,state = sstate

    res = iterate(pit, state)
    if res isa Nothing
        return nothing
    else
        c, newstate = res
        return cg[c...], (pit,newstate)
    end

end

function isdone(cg::CGrid{d}, sstate) where d
    pit, state= sstate
    return isdone(pit, state)
end



## TODO: check the following


function Base.iterate(g::PGrid{G1, G2, d}) where G1 where G2 where d
    x = g.g1[1]
    y = g.g2[1]
    return (SVector{d, Float64}(x...,y...),(y,1,1))
end

function Base.iterate(g::PGrid{G1,G2,d},state) where G1 where G2 where d
    y,i,j=state
    if i<length(g.g1)
        i += 1
        x = g.g1[i]
        return (SVector{d,Float64}(x..., y...), (y,i,j))
    else
        if j==length(g.g2)
            return nothing
        else
            j += 1
            i = 1
            x = g.g1[i]
            y = g.g2[j]
            return (SVector{d,Float64}(x..., y...), (y,i,j))
        end
    end
end

enum(grid::PGrid) = (
    QP((c[1],c[2]), grid[c...]) for c in Iterators.product(1:length(grid.g1), 1:length(grid.g2))
)

