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

"Jacobian of the residuals at grid point `s` w.r.t. the controls `x`."
dF_1(dmodel::ADModel, s::QP, x::SVector, φ) =
    ForwardDiff.jacobian(u->F(dmodel, s, u, φ), x)

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


## Vectorized, allocating

F(dmodel::ADModel, x::GArray, φ) =
    GArray(dmodel.grid, [F_kernel(n, dmodel, x, φ) for n in 1:length(x)])

dF_1(dmodel::ADModel, x::GArray, φ) =
    GArray(dmodel.grid, [dF_1_kernel(n, dmodel, x, φ) for n in 1:length(x)])

dF_2(dmodel::ADModel, x::GArray, φ) =
    LL(dmodel.grid, [dF_2_kernel(n, dmodel, x, φ) for n in 1:length(x)], deepcopy(φ))
