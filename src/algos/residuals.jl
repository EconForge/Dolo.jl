# Residuals of the arbitrage equations, and their derivatives.
#
# For controls x on the grid, and a decision rule φ for next period, the residual
# at grid point s is
#
#     F(s, x, φ) = complementarities( Σ_S w(S) f(s, x, S, φ(S)) )
#
# Each quantity is defined three times, from the same pointwise implementation:
# - pointwise, at one grid point `s` (canonical and non-allocating)
# - as a kernel, at grid point `n`, for `pointwise!`
# - vectorized: in place on the whole grid with a given engine, or allocating.


## Pointwise

"Expectation of the arbitrage equations at grid point `s`."
function expected_arbitrage(dmodel::ADModel, s::QP, x::SVector{d,T}, φ) where d where T
    r = zero(SVector{d,T})
    for (w,S) in τ(dmodel, s, x)
        r += w*arbitrage(dmodel, s, x, S, φ(S))
    end
    return r
end

"Residuals at grid point `s`."
F(dmodel::ADModel, s::QP, x::SVector, φ) =
    complementarities(dmodel.model, s, x, expected_arbitrage(dmodel, s, x, φ))

"""
Decision rule `φ` replaced by its first order expansion at the point where it is evaluated:
`φ(S₀) + φ'(S₀)(S-S₀)`, where S₀ is the value of S. When S is a vector of dual numbers,
it has the same value and the same derivative as φ(S), but φ is only differentiated w.r.t.
the endogenous states (with as many partials), not through the dual numbers of S.
"""
struct Linearized{F}
    φ::F
end

struct LinearizedTag end

function (l::Linearized)(S::QP)
    i, z = S.loc
    z0 = ForwardDiff.value.(z)
    T = eltype(z0)
    d = length(z0)
    Tag = ForwardDiff.Tag{LinearizedTag, T}
    zd = SVector(ntuple(k -> ForwardDiff.Dual{Tag}(z0[k], ntuple(j -> T(j == k), Val(d))), Val(d)))
    r = l.φ(i, zd)
    X0 = ForwardDiff.value.(r)
    G = SMatrix{length(r), d}(ntuple(n -> ForwardDiff.partials(r[(n-1) % length(r) + 1], (n-1) ÷ length(r) + 1), Val(length(r)*d)))
    return X0 + G*(z - z0)
end

"Jacobian of the residuals at grid point `s` w.r.t. the controls `x`."
dF_1(dmodel::ADModel, s::QP, x::SVector, φ) =
    ForwardDiff.jacobian(u->F(dmodel, s, u, φ), x)

# with exogenous Markov states, φ is a spline in the endogenous states: it is linearized
dF_1(dmodel::ADModel, s::QP, x::SVector, φ::DFun{<:Any, <:GArray{<:PGrid{<:SGrid, <:CGrid}}}) =
    ForwardDiff.jacobian(u->F(dmodel, s, u, Linearized(φ)), x)

"""
Derivative of the residuals at grid point `s` w.r.t. the decision rule `φ`:
a tuple of blocks `(;F_x, S)`, one for each future state `S`, such that `dF = Σ F_x dφ(S)`.
"""
function dF_2(dmodel::ADModel, s::QP, x::SVector, φ, ::Val{true})
    r_F = ForwardDiff.jacobian(
        r->complementarities(dmodel.model, s, x, r),
        expected_arbitrage(dmodel, s, x, φ)
    )
    return tuple((
        (; F_x = w*r_F*ForwardDiff.jacobian(u->arbitrage(dmodel, s, x, S, u), φ(S)), S = S)
        for (w,S) in τ(dmodel, s, x)
    )...)
end

function dF_2(dmodel::ADModel, s::QP, x::SVector, φ, ::Val{false})
    return tuple((
        (; F_x = w*ForwardDiff.jacobian(u->arbitrage(dmodel, s, x, S, u), φ(S)), S = S)
        for (w,S) in τ(dmodel, s, x)
    )...)
end

dF_2(dmodel::ADModel, s::QP, x::SVector, φ) = dF_2(dmodel, s, x, φ, Val(has_complementarities(dmodel.model)))


## Kernels (at grid point n)

F_kernel(n, dmodel, x, φ) = F(dmodel, point(dmodel.grid, n), x[n], φ)
dF_1_kernel(n, dmodel, x, φ) = dF_1(dmodel, point(dmodel.grid, n), x[n], φ)
dF_2_kernel(n, dmodel, x, φ, v::Val) = dF_2(dmodel, point(dmodel.grid, n), x[n], φ, v)
dF_2_kernel(n, dmodel, x, φ) = dF_2_kernel(n, dmodel, x, φ, Val(has_complementarities(dmodel.model)))


## Vectorized, in place

"`r ← F(x, φ)`"
F!(r::GArray, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (pointwise!(F_kernel, r.data, engine, dmodel, x, φ); r)

"`J ← ∂F/∂x (x, φ)`, block diagonal"
dF_1!(J::GArray, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (pointwise!(dF_1_kernel, J.data, engine, dmodel, x, φ); J)

"`(r, J) ← F(x, φ)` and `∂F/∂x (x, φ)` simultaneously in a single pass"
struct FdF1Tag end

function F_and_dF_1_pointwise(dmodel::ADModel, s::QP, xn::SVector{d, T}, φ) where {d, T}
    Tag = ForwardDiff.Tag{FdF1Tag, T}
    xd = SVector(ntuple(k -> ForwardDiff.Dual{Tag}(xn[k], ntuple(j -> T(j == k), Val(d))), Val(d)))
    lin_φ = Linearized(φ)
    res_dual = F(dmodel, s, xd, lin_φ)
    r = SVector(ntuple(i -> ForwardDiff.value(res_dual[i]), Val(d)))
    J = SMatrix{d, d, T, d*d}(ntuple(idx -> begin
        row = (idx - 1) % d + 1
        col = (idx - 1) ÷ d + 1
        ForwardDiff.partials(res_dual[row], col)
    end, Val(d*d)))
    return r, J
end

function F_and_dF_1_kernel!(n, r, J, dmodel, x, φ)
    s = point(dmodel.grid, n)
    rn, Jn = F_and_dF_1_pointwise(dmodel, s, x[n], φ)
    r[n] = rn
    J[n] = Jn
end

function F_and_dF_1!(r::GArray, J::GArray, dmodel::ADModel, x::GArray, φ, engine=Serial())
    foreach_point!(F_and_dF_1_kernel!, length(x), engine, r.data, J.data, dmodel, x, φ)
    return r, J
end

"`L ← ∂F/∂φ (x, φ)`"
dF_2!(L::LL, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (pointwise!(dF_2_kernel, L.D, engine, dmodel, x, φ, Val(has_complementarities(dmodel.model))); L)

# with a BlockVector, each block is written as soon as it is computed
function dF_2_blocks!(n, D, dmodel, x, φ, ::Val{true})
    s = point(dmodel.grid, n)
    xn = x[n]
    r_F = ForwardDiff.jacobian(
        r->complementarities(dmodel.model, s, xn, r),
        expected_arbitrage(dmodel, s, xn, φ)
    )
    for (k, (w,S)) in enumerate(τ(dmodel, s, xn))
        F_x = w*r_F*ForwardDiff.jacobian(u->arbitrage(dmodel, s, xn, S, u), φ(S))
        setblock!(D, (; F_x, S), n, k)
    end
end

function dF_2_blocks!(n, D, dmodel, x, φ, ::Val{false})
    s = point(dmodel.grid, n)
    xn = x[n]
    for (k, (w,S)) in enumerate(τ(dmodel, s, xn))
        F_x = w*ForwardDiff.jacobian(u->arbitrage(dmodel, s, xn, S, u), φ(S))
        setblock!(D, (; F_x, S), n, k)
    end
end

dF_2_blocks!(n, D, dmodel, x, φ) = dF_2_blocks!(n, D, dmodel, x, φ, Val(has_complementarities(dmodel.model)))

dF_2!(L::LL{<:Any,<:BlockVector}, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (foreach_point!(dF_2_blocks!, length(L.D), engine, L.D, dmodel, x, φ, Val(has_complementarities(dmodel.model))); L)


## Vectorized, allocating

F(dmodel::ADModel, x::GArray, φ) =
    GArray(dmodel.grid, [F_kernel(n, dmodel, x, φ) for n in 1:length(x)])

dF_1(dmodel::ADModel, x::GArray, φ) =
    GArray(dmodel.grid, [dF_1_kernel(n, dmodel, x, φ) for n in 1:length(x)])

dF_2(dmodel::ADModel, x::GArray, φ) =
    LL(dmodel.grid, [dF_2_kernel(n, dmodel, x, φ) for n in 1:length(x)], deepcopy(φ))
