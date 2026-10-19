# Linear operator L acting on vectors v defined on the grid:
#
#     (L.v)[n] = Σ_k F_x[n,k] dv(S[n,k])
#
# where dv interpolates v. It is stored as one tuple of blocks (;F_x, S) per grid point.
# The derivative of the residuals w.r.t. the decision rule has this form (see `dF_2`).

import Base: *, \, -
using LinearMaps
using BlockDiagonals
using LinearAlgebra: UniformScaling
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


## Kernels (at grid point n)

"`(L.v)[n]`, where `dφ` interpolates `v`."
# hot loop of Neumann/GMRES iterations: inlining and accumulating from zero are measurably faster
@inline function apply_kernel(n, D, dφ)
    t = zero(D[n][1].F_x[:, 1])
    for k in 1:length(D[n])
        (; F_x, S) = D[n][k]
        t += F_x * dφ(S)
    end
    return t
end

"Blocks of `α L` at grid point `n`."
scale_kernel(n, D, α) = map(b -> (; F_x = α*b.F_x, S = b.S), D[n])

"Blocks of `J⁻¹ L` at grid point `n`, for a block-diagonal `J`."
ldiv_kernel(n, D, J) = map(b -> (; F_x = J[n]\b.F_x, S = b.S), D[n])


## In place operations

"`L ← α L`"
mul!(L::LL, α::Number, engine=Serial()) = (pointwise!(scale_kernel, L.D, engine, L.D, α); L)

"`L ← J⁻¹ L`, where `J` is block diagonal"
ldiv!(J::GArray, L::LL, engine=Serial()) = (pointwise!(ldiv_kernel, L.D, engine, L.D, J.data); L)

"`r ← L.v`"
function mul!(r::GArray, L::LL, v::GArray, engine=Serial())
    fit!(L.φ, v)
    pointwise!(apply_kernel, r.data, engine, L.D, L.φ)
    return r
end

# same, on flat vectors
function mul!(r::AbstractVector{T}, L::LL, v::AbstractVector{T}, engine=Serial()) where T<:Union{Float32, Float64}
    p, q = size(L.D[1][1].F_x)
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
    p, q = size(L.D[1][1].F_x)
    N = length(L.grid)
    LinearMap(
        (v,u)-> mul!(v, L, u, engine),
        p*N,
        q*N;
        ismutating=true
    )
end
