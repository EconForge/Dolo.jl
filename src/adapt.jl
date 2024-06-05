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

function adapt_structure(to, v::GArray{G,V}) where G where V
    data = adapt(to, v.data)
    GArray(
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
    LL(L.grid, adapt(to, L.D), adapt(to, L.φ))
end

maxabs(u::Number, v::Number) = abs(max(u,v))
