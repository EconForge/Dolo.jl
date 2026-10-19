import Adapt: adapt, adapt_structure

adapt_structure(to, s::SplineInterpolator{G,C,k}) where G where C where k = let
    θθ = adapt(to, s.θ)
    CC = typeof(θθ)
    SplineInterpolator{G,CC,k}(s.grid, θθ)
end

function adapt_structure(to, v::GVector{G,V}) where G where V
    data = adapt(to, v.data)
    GVector{G, typeof(data)}(
        v.grid,
        data
    )
end

function adapt_structure(to, f::DFun{Dom, Gar, Itp, vars}) where Dom where Gar where Itp where vars
    itp = adapt(to, f.itp)
    values = adapt(to, f.values)
    DFun{Dom, typeof(values), typeof(itp), vars}(
        f.domain,
        values,
        itp
    )
end

function adapt_structure(to, L::LL{G,D,F}) where G where D where F
    LL(L.grid, adapt_blocks(to, L.D), adapt(to, L.φ))
end

# blocks moved off the CPU are stored component by component (see `BlockVector`)
function adapt_blocks(to, D::Vector{<:BlockTuple})
    on_cpu = adapt(to, zeros(Int32, 1)) isa Array
    return on_cpu ? adapt(to, D) : adapt(to, BlockVector(D))
end
adapt_blocks(to, D) = adapt(to, D)

adapt_structure(to, B::BlockVector{Tf,p,q,pq,d,ds,K}) where {Tf,p,q,pq,d,ds,K} =
    BlockVector{Tf,p,q,d,ds,K}(adapt(to, B.F), adapt(to, B.I), adapt(to, B.X), adapt(to, B.V))

maxabs(u::Number, v::Number) = abs(max(u,v))
