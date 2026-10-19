# Linear operator L acting on vectors v defined on the grid:
#
#     (L.v)[n] = Σ_k F_x[n,k] dv(S[n,k])
#
# where dv interpolates v. It is stored as one tuple of blocks (;F_x, S) per grid point.
# The derivative of the residuals w.r.t. the decision rule has this form (see `dF_2`).

import Base: *, \, -
using LinearMaps
using BlockDiagonals
using LinearAlgebra: UniformScaling, lu
import LinearAlgebra: mul!, ldiv!


struct LL{G,D_,F}
    grid::G
    D::D_     # D[n] = ((;F_x, S), ...)
    φ::F      # interpolator for the vectors L is applied to
end

const LF = LL

struct I_L{_L}
    L::_L
end

-(::UniformScaling, L::LL) = I_L(L)
*(A::I_L, b) = b-A.L*b

convert(::Type{LinearMap}, IL::I_L) = I - convert(LinearMap,IL.L)


## Storage of the blocks

# On a GPU, the blocks are stored component by component (`BlockVector`), so that
# neighbouring threads (grid points) read neighbouring memory addresses. With one tuple
# of blocks per grid point, they would read addresses several kilobytes apart.

"""
    BlockVector

Vector of tuples of `K` blocks `(;F_x, S)`, with `F_x::SMatrix{p,q}` and `S::QP` a point
of a `SGrid × CGrid` grid, stored as arrays indexed by `[n, k, component]`.
Reading or writing an element works as for a vector of tuples.
"""
struct BlockVector{Tf,p,q,pq,d,ds,K,AF,AI} <: AbstractVector{NTuple{K, NamedTuple{(:F_x, :S), Tuple{SMatrix{p,q,Tf,pq}, QP{Tuple{Int,SVector{d,Tf}}, SVector{ds,Tf}}}}}}
    F::AF     # F[n,k,:] = F_x[:]
    I::AI     # I[n,k] = S.loc[1] (exogenous node)
    X::AF     # X[n,k,:] = S.loc[2] (endogenous coordinates)
    V::AF     # V[n,k,:] = S.val (all coordinates)
end

BlockVector{Tf,p,q,d,ds,K}(F::AF, I::AI, X::AF, V::AF) where {Tf,p,q,d,ds,K,AF,AI} =
    BlockVector{Tf,p,q,p*q,d,ds,K,AF,AI}(F, I, X, V)

const BlockTuple{K,Tf,p,q,pq,d,ds} = NTuple{K, NamedTuple{(:F_x, :S), Tuple{SMatrix{p,q,Tf,pq}, QP{Tuple{Int,SVector{d,Tf}}, SVector{ds,Tf}}}}}

function BlockVector(D::AbstractVector{<:BlockTuple{K,Tf,p,q,pq,d,ds}}) where {K,Tf,p,q,pq,d,ds}
    N = length(D)
    B = BlockVector{Tf,p,q,d,ds,K}(zeros(Tf, N, K, p*q), zeros(Int32, N, K), zeros(Tf, N, K, d), zeros(Tf, N, K, ds))
    for n in 1:N
        B[n] = D[n]
    end
    return B
end

Base.size(B::BlockVector) = (size(B.I, 1),)
Base.IndexStyle(::Type{<:BlockVector}) = IndexLinear()

"Block `k` at grid point `n`."
@inline function getblock(B::BlockVector{Tf,p,q,pq,d,ds}, n::Int, k::Int) where {Tf,p,q,pq,d,ds}
    F_x = SMatrix{p,q,Tf,pq}(ntuple(j -> @inbounds(B.F[n,k,j]), Val(pq)))
    x = SVector{d,Tf}(ntuple(j -> @inbounds(B.X[n,k,j]), Val(d)))
    v = SVector{ds,Tf}(ntuple(j -> @inbounds(B.V[n,k,j]), Val(ds)))
    return (; F_x, S = QP((Int(@inbounds B.I[n,k]), x), v))
end

@inline function setblock!(B::BlockVector{Tf,p,q,pq,d,ds}, b, n::Int, k::Int) where {Tf,p,q,pq,d,ds}
    for j in 1:pq
        @inbounds B.F[n,k,j] = b.F_x[j]
    end
    @inbounds B.I[n,k] = b.S.loc[1]
    for j in 1:d
        @inbounds B.X[n,k,j] = b.S.loc[2][j]
    end
    for j in 1:ds
        @inbounds B.V[n,k,j] = b.S.val[j]
    end
    return B
end

"Set the matrix `F_x` of block `k` at grid point `n`."
@inline function setF!(B::BlockVector{Tf,p,q,pq}, F_x, n::Int, k::Int) where {Tf,p,q,pq}
    for j in 1:pq
        @inbounds B.F[n,k,j] = F_x[j]
    end
end

@inline Base.getindex(B::BlockVector{Tf,p,q,pq,d,ds,K}, n::Int) where {Tf,p,q,pq,d,ds,K} =
    ntuple(k -> getblock(B, n, k), Val(K))

@inline function Base.setindex!(B::BlockVector{Tf,p,q,pq,d,ds,K}, v, n::Int) where {Tf,p,q,pq,d,ds,K}
    for k in 1:K
        setblock!(B, v[k], n, k)
    end
    return B
end

Base.copyto!(A::BlockVector, B::BlockVector) = (copyto!(A.F, B.F); copyto!(A.I, B.I); copyto!(A.X, B.X); copyto!(A.V, B.V); A)
Base.similar(B::BlockVector{Tf,p,q,pq,d,ds,K}) where {Tf,p,q,pq,d,ds,K} = BlockVector{Tf,p,q,d,ds,K}(similar(B.F), similar(B.I), similar(B.X), similar(B.V))
Base.Array(B::BlockVector) = collect(adapt(Array, B))

# for vectors of tuples of blocks
@inline getblock(D::AbstractVector, n, k) = @inbounds D[n][k]

"Number of blocks at each grid point."
nblocks(D::AbstractVector) = fieldcount(eltype(D))
nblocks(::BlockVector{Tf,p,q,pq,d,ds,K}) where {Tf,p,q,pq,d,ds,K} = K

"Size `(p, q)` of the blocks `F_x`."
blocksize(D::AbstractVector) = size(fieldtype(fieldtype(eltype(D), 1), :F_x))
blocksize(::BlockVector{Tf,p,q}) where {Tf,p,q} = (p, q)


## Kernels (at grid point n)

"`(L.v)[n]`, where `dφ` interpolates `v`."
# hot loop of Neumann/GMRES iterations: inlining and accumulating from zero are measurably faster
# blocks are read one at a time (with a BlockVector, D[n] would load all of them)
@inline function apply_kernel(n, D, dφ)
    (; F_x, S) = getblock(D, n, 1)
    t = F_x * dφ(S)
    for k in 2:nblocks(D)
        (; F_x, S) = getblock(D, n, k)
        t += F_x * dφ(S)
    end
    return t
end

"Blocks of `α L` at grid point `n`."
scale_kernel(n, D, α) = map(b -> (; F_x = α*b.F_x, S = b.S), D[n])

"Blocks of `J⁻¹ L` at grid point `n`, for a block-diagonal `J`."
# J[n] is factorized once for all blocks
function ldiv_kernel(n, D, J)
    F = lu(J[n])
    map(b -> (; F_x = F\b.F_x, S = b.S), D[n])
end


# same, one block at a time, for a BlockVector
function scale_blocks!(n, D, α)
    for k in 1:nblocks(D)
        setF!(D, α*getblock(D, n, k).F_x, n, k)
    end
end

function ldiv_blocks!(n, D, J)
    F = lu(J[n])
    for k in 1:nblocks(D)
        setF!(D, F\getblock(D, n, k).F_x, n, k)
    end
end


## In place operations

"`L ← α L`"
mul!(L::LL, α::Number, engine=Serial()) = (pointwise!(scale_kernel, L.D, engine, L.D, α); L)
mul!(L::LL{<:Any,<:BlockVector}, α::Number, engine=Serial()) = (foreach_point!(scale_blocks!, length(L.D), engine, L.D, α); L)

"`L ← J⁻¹ L`, where `J` is block diagonal"
ldiv!(J::GArray, L::LL, engine=Serial()) = (pointwise!(ldiv_kernel, L.D, engine, L.D, J.data); L)
ldiv!(J::GArray, L::LL{<:Any,<:BlockVector}, engine=Serial()) = (foreach_point!(ldiv_blocks!, length(L.D), engine, L.D, J.data); L)

"`r ← L.v`"
function mul!(r::GArray, L::LL, v::GArray, engine=Serial())
    fit!(L.φ, v)
    pointwise!(apply_kernel, r.data, engine, L.D, L.φ)
    return r
end

# same, on flat vectors
function mul!(r::AbstractVector{T}, L::LL, v::AbstractVector{T}, engine=Serial()) where T<:Union{Float32, Float64}
    p, q = blocksize(L.D)
    mul!(
        GArray(L.grid, reinterpret(SVector{p,T}, r)),
        L,
        GArray(L.grid, reinterpret(SVector{q,T}, v)),
        engine
    )
    return r
end


## Allocating versions

*(L::LL, v::GArray) = mul!(duplicate(v), L, v)
*(L::LL, v::AbstractVector{T}) where T<:Union{Float32, Float64} = mul!(similar(v), L, v)

\(J::GArray, L::LL) = ldiv!(J, deepcopy(L))


## (I-L)⁻¹ by Neumann series

function neumann(L::LL, r0; K=1000, τ_η=1e-10, engine=Serial())
    dx = duplicate(r0)
    neumann!(dx, L, r0, engine; K=K, τ_η=τ_η)
    return dx
end

# computes dx = (I-L)^(-1) r0 = Σ_k L^k r0, returns the number of terms
function neumann!(dx, L::LL, r0, engine=Serial(); mem=(;du=duplicate(r0), dv=duplicate(r0)), K=1000, τ_η=1e-10, tmr=NO_TIMER)

    (; du, dv) = mem

    dx.data .= r0.data
    du.data .= r0.data
    dv.data .= r0.data

    k = 0
    while k < K
        k += 1
        @timeit tmr "eval L" mul!(du, L, dv, engine)
        dx.data .+= du.data
        η = norm(du)
        if η<τ_η
            break
        end
        du,dv=dv,du
    end

    return k

end


## Conversions

convert(::Type{BlockDiagonal}, ga::GArray{G, Vector{T}}) where T where G = BlockDiagonal(ga.data)

convert(::Type{Matrix}, L::LL) = convert(Matrix, convert(LinearMap,L))

function convert(::Type{LinearMap}, L::LL, engine=Serial())
    p, q = blocksize(L.D)
    N = length(L.grid)
    LinearMap(
        (v,u)-> mul!(v, L, u, engine),
        p*N,
        q*N;
        ismutating=true
    )
end
