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
function dF_2(dmodel::ADModel, s::QP, x::SVector, φ)
    # complementarities are applied to the expectation, not to each term
    r_F = ForwardDiff.jacobian(
        r->complementarities(dmodel.model, s, x, r),
        expected_arbitrage(dmodel, s, x, φ)
    )
    return tuple((
        (; F_x = w*r_F*ForwardDiff.jacobian(u->arbitrage(dmodel, s, x, S, u), φ(S)), S = S)
        for (w,S) in τ(dmodel, s, x)
    )...)
end


## Kernels (at grid point n)

F_kernel(n, dmodel, x, φ) = F(dmodel, point(dmodel.grid, n), x[n], φ)
dF_1_kernel(n, dmodel, x, φ) = dF_1(dmodel, point(dmodel.grid, n), x[n], φ)
dF_2_kernel(n, dmodel, x, φ) = dF_2(dmodel, point(dmodel.grid, n), x[n], φ)


## Vectorized, in place

"`r ← F(x, φ)`"
F!(r::GArray, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (pointwise!(F_kernel, r.data, engine, dmodel, x, φ); r)

"`J ← ∂F/∂x (x, φ)`, block diagonal"
dF_1!(J::GArray, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (pointwise!(dF_1_kernel, J.data, engine, dmodel, x, φ); J)

"`L ← ∂F/∂φ (x, φ)`"
dF_2!(L::LL, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (pointwise!(dF_2_kernel, L.D, engine, dmodel, x, φ); L)

# with a BlockVector, each block is written as soon as it is computed
function dF_2_blocks!(n, D, dmodel, x, φ)
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

dF_2!(L::LL{<:Any,<:BlockVector}, dmodel::ADModel, x::GArray, φ, engine=Serial()) =
    (foreach_point!(dF_2_blocks!, length(L.D), engine, L.D, dmodel, x, φ); L)


## Vectorized, allocating

F(dmodel::ADModel, x::GArray, φ) =
    GArray(dmodel.grid, [F_kernel(n, dmodel, x, φ) for n in 1:length(x)])

dF_1(dmodel::ADModel, x::GArray, φ) =
    GArray(dmodel.grid, [dF_1_kernel(n, dmodel, x, φ) for n in 1:length(x)])

dF_2(dmodel::ADModel, x::GArray, φ) =
    LL(dmodel.grid, [dF_2_kernel(n, dmodel, x, φ) for n in 1:length(x)], deepcopy(φ))
