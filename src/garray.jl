## GArrays

struct GArray{G,U}
    grid::G
    data::U
end

const GVector{G,T} = GArray{G,T}
const GDist{G,Tf} = GArray{G, Vector{Tf}}

GDist(g::G, v::Vector{Tf}) where G where Tf = GArray{G,Vector{Tf}}(g, v)


norm(v::GArray) = maximum(u->maximum(abs, u), v.data)


distance(u::SVector, v::SVector) = maximum(abs, u-v)
# distance(t::Tuple{<:SVector, <:SVector}) = distance(t[1],t[2])

# distance(u::GArray, v::GArray) = maximum( (u,v)->distance(u,v), zip(u.data, v.data))

distance(u::GVector, v::GVector) = maximum( t->distance(t[1],t[2]), zip(u.data, v.data))

import Base: size
size(a::GArray{PGrid{G1, G2, d}, T}) where G1 where G2 where d where T = (length(a.grid.grids[1]), length(a.grid.grids[2]))


# TODO: check
getindex(a::GArray{PGrid{G1, G2, d}, T}, i::Int, j::Int) where G1 where G2 where d where T = a.data[ i + length(a.grid.grids[1])*(j-1)]
getindex(a::GArray{PGrid{G1, G2, d}, T}, i::Int, ::Colon) where G1 where G2 where d where T = [a[i,j] for j=1:length(a.grid.grids[2])]
getindex(a::GArray{PGrid{G1, G2, d}, T}, ::Colon, j::Int) where G1 where G2 where d where T = [a[i,j] for i=1:length(a.grid.grids[1])]

getindex(a::GArray{PGrid{G1, G2, d}, T}, ::Colon) where G1 where G2 where d where T = a.data


# TODO: check
function setindex!(a::GArray{PGrid{G1, G2, d}, T}, v, i::Int, j::Int) where G1 where G2 where d where T 
    (a.data[ i + length(a.grid.grids[1])*(j-1)] = v)
    nothing
end

getprecision(g::GArray{G,T}) where G where T = eltype(g).types[1].types[1]


eltype(g::GArray{G,T}) where G where T = eltype(T)

# warning: these functions don't copy any data
function ravel(g::GArray) 
    Tf = eltype(g.data[1])
    reinterpret(Tf, g.data)
end
unravel(g::GArray, x) = GArray(
    g.grid,
    reinterpret(eltype(g), x)
)

iterate(g::GArray) = iterate(g.data)
iterate(g::GArray, i) = iterate(g.data, i)

length(g::GArray) = length(g.data)

getindex(g::GArray{G,T}, i::Int64) where G where T = g.data[i]
setindex!(g::GArray, x, i) = (g.data[i] = x)

function setindex!(
        g::GArray{<:Dolo.CGrid{2},Vector{T}},
        v::T,
        i::Int64, 
        j::Int64) where T
    setindex!(g.data, v, to__linear_index(g.grid, i, j))
end

getindex(g::GArray{G,T}, inds::Vararg{Int64, d}) where G<:AGrid{d} where d where T = g.data[to__linear_index(g.grid, inds)]
getindex(g::GArray{G,T}, inds::Int64) where G<:AGrid{d} where d where T = g.data[inds]

# interpolating indexing
function (xa::GArray{PGrid{G1, G2, d}, T})(i::Int64, p::SVector{d2, U}) where G1<:SGrid where G2<:CGrid{d2} where d where d2 where T where U
    g1 = xa.grid.grids[1]
    g2 = xa.grid.grids[2]
    dims = tuple(length(g1), (e[3] for e in g2.ranges)... )
    # ranges = tuple( (range(e...) for e in g2.ranges)... )
    # v = view(reshape(xa.data, dims),i,:)
    # v = view(reshape(xa.data, dims),i,:)
    v = reshape(view(xa.data, :), dims) ### Weird but should not allocate
    # v = view(xa.data, :) #1d only
    res = interp(g2.ranges, view(v,i,:), p...)
    res
end

(xa::GArray{PGrid{G1, G2, d}, T})(i::Int64, j::Int64) where G1 where G2 where d where T  = xa[i,j]
# (xa::GArray{PGrid{G1, G2, d}, T})(S::Tuple{Tuple{Int64}, U}) where G1 where G2 where U where d where T = xa(S[1][1],S[2][2])
(xa::GArray{PGrid{G1, G2, d}, T})(S::Tuple{Tuple{Int64, Int64}, U}) where G1 where G2 where U where d where T = xa[S[1]...]

function (xa::GArray{PGrid{G1, G2, d}, T})(S::Tuple{Tuple{Int64}, U}) where G1 where G2 where U<:SVector where d where T
    #### TODO: replace
    n_x =  ndims(xa.grid.grids[2])
    V = S[2]
    n = length(V)
    s = SVector((V[i] for i=n-n_x+1:n)...)
    xa(S[1][1],s)
end



# enum(g::GArray) = enumerate(g)

import Base: *, \, +, -, /

*(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.*B.data)
\(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.\B.data)
/(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data./B.data)
+(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.+B.data)
-(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.-B.data)

*(A::GArray{G,T}, x::Number) where G where T  = GArray(A.grid, A.data .* x)
*(x::Number, A::GArray{G,T}) where G where T = GArray(A.grid, x .* A.data)


*(A::GArray{G,Vector{T}}, x::SVector{q, Tf}) where G where T <:SMatrix{p, q, Tf, n}  where p where q where n where Tf = GArray(A.grid, [M*x for M in A.data])


*(A::GArray{G,T}, B::AbstractArray{Tf}) where G where T <:SMatrix{p, q, Tf, n}  where p where q where n where Tf = 
    ravel(
        GArray(
            A.grid,
            A.data .* reinterpret(SVector{q, Tf}, B)
        )
    )


import Base: convert

function Base.convert(::Type{Matrix}, A::GArray{G,Vector{T}}) where G where T <:SMatrix{p, q, Tf, k}  where p where q where k where Tf
    N = length(A.data)
    n0 = N*p
    n1 = N*q
    M = zeros(n0, n1)
    for n=1:N
        i = p*(n-1)+1
        j = q*(n-1)+1
        M[i:i+p-1, j:j+q-1] = A.data[n]
    end
    return M
end



### TODO

duplicate(g::GArray{G,T}) where G where T = GArray(g.grid, deepcopy(g.data))