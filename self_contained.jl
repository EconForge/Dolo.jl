using Test
using StaticArrays
import Base: getindex, length

struct CGrid{d}
    ranges::NTuple{d, Tuple{Float64, Float64, Int64}}
end


getindex(g::CGrid{d}, inds::Vararg{Int64,d}) where d = SVector{d}(
    ( 
        ( g.ranges[i][1] + (g.ranges[i][2]-g.ranges[i][1])*( (inds[i]-1)/(g.ranges[i][3]-1)) )
        for i=1:d
    )
)


length(g::CGrid{1})= g.ranges[1][3]
enum(g::CGrid{1})= ( (;loc=(i,), val=g[i]) for i=1:length(g))
enum(g::CGrid{d}) where d = ((;loc=c, val=g[c...]) for c in Iterators.product( tuple( ((1:r[3]) for r in g.ranges)... )... ) )


d = 3

NN = 20

vars =  tuple( (Symbol("d$i") for i in 1:d)... )
sizes = tuple( (NN+i+1 for i in 1:d)...)
args = tuple( ( (0.0, 1.0*i, sizes[i]) for i in 1:d)... )
cg = CGrid( args )

cg.ranges

eeenum(cg::CGrid{d}) where d =  (cg[c...] for c in Iterators.product( tuple( ((1:r[3]) for r in cg.ranges)... )... ) )


pit = Iterators.ProductIterator{Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}((1:22, 1:23, 1:24))

v1, s2 = iterate(pit)
v2, s3 = iterate(pit, s2)

import Base: iterate
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

eltype(cg::CGrid{d}) where d = SVector{d, Float64}

import Base: length, eltype, isdone
length(cg::CGrid{d}) where d = prod( (e[3] for e in cg.ranges))

function isdone(cg::CGrid{d}, sstate) where d
    pit, state= sstate
    return isdone(pit, state)
end

function check(cg)
    sum(sum( cg ))
end

@code_warntype check(cg)
check(cg)
@allocations check(cg) == 1




# function fun2(cg::CGrid{d}) where d
#     if sum( sum( (sum(loc...)*val) for (;loc, val) in enum(cg) ) ) < -Inf
#         print("")
#     end
# end

function fun3(cg::CGrid{d}) where d
    if sum( sum( (cg[loc...].*val) for (;loc, val) in enum(cg) ) ) < -Inf
        print("")
    end
end


# @code_warntype fun(cg)
eltype(enum(cg))

# fun2(cg)
# println( @allocated fun2(cg) )


fun3(cg)
println( @allocated fun3(cg) )
