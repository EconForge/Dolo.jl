using StaticArrays
import Dolo: transition, arbitrage, complementarities, hypeof
import Dolo: ⫫

# Financial autarky (Coeurdacier et al., 2020): two countries, no asset trade.
# Epstein-Zin preferences, convex capital adjustment costs, nonnegative investment.
# Same as integration_B without the bond (bruteforce/models/integration_A.yaml).

model = let

    name = :integration_A

    β = 0.96
    θ = 0.3
    δ = 0.08
    γ = 4.0
    ψ = 4.0
    ρ_A = 0.9
    ξ = 0.2
    a2 = δ^ξ
    a1 = δ - a2/(1-ξ)*δ^(1-ξ)
    σ_A_1 = 0.025
    σ_A_2 = 0.05
    ζ = 0.0     # correlation of productivity shocks

    kmin = 2.0
    kmax = 10.0

    # steady-state (identical countries)
    A_1 = 0.0
    A_2 = 0.0
    i_1 = ((1/β-(1-δ))/θ)^(1/(θ-1)) * δ
    i_2 = i_1
    k_1 = i_1/δ
    k_2 = k_1
    c_1 = k_1^θ - δ*k_1
    c_2 = c_1
    w_1 = c_1
    w_2 = c_2
    ew_1 = c_1
    ew_2 = c_2

    calibration = (;β, θ, δ, γ, ψ, ρ_A, ξ, a1, a2, σ_A_1, σ_A_2, ζ,
        A_1, A_2, k_1, k_2, i_1, i_2, ew_1, ew_2, w_1, w_2)

    states = Dolo.CartesianSpace(;
        :A_1 => (-3σ_A_1/sqrt(1-ρ_A^2), 3σ_A_1/sqrt(1-ρ_A^2)),
        :A_2 => (-3σ_A_2/sqrt(1-ρ_A^2), 3σ_A_2/sqrt(1-ρ_A^2)),
        :k_1 => (kmin, kmax),
        :k_2 => (kmin, kmax)
    )

    controls = Dolo.CartesianSpace(;
        :i_1 => (0.0, Inf),
        :i_2 => (0.0, Inf),
        :ew_1 => (0.0, Inf),
        :ew_2 => (0.0, Inf),
        :w_1 => (0.0, Inf),
        :w_2 => (0.0, Inf)
    )

    Σ = @SMatrix [σ_A_1^2 ζ*σ_A_1*σ_A_2; ζ*σ_A_1*σ_A_2 σ_A_2^2]
    process = Dolo.VAR1( (:A_1, :A_2), ρ_A, Σ )

    Dolo.YModel(name, states, controls, process, calibration)

end


function intermediate(model::hypeof(model), s::NamedTuple, x::NamedTuple)

    p = model.calibration

    y_1 = s.k_1^p.θ*exp(s.A_1)
    y_2 = s.k_2^p.θ*exp(s.A_2)
    c_1 = y_1 - x.i_1
    c_2 = y_2 - x.i_2
    Φ_1 = p.a1 + p.a2/(1-p.ξ)*(x.i_1/s.k_1)^(1-p.ξ)
    Φ_2 = p.a1 + p.a2/(1-p.ξ)*(x.i_2/s.k_2)^(1-p.ξ)
    Φ_1__p = p.a2*(x.i_1/s.k_1)^(-p.ξ)
    Φ_2__p = p.a2*(x.i_2/s.k_2)^(-p.ξ)

    return (;y_1, y_2, c_1, c_2, Φ_1, Φ_2, Φ_1__p, Φ_2__p)

end


function transition(model::hypeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)

    p = model.calibration
    y = intermediate(model, s, x)

    k_1 = ( (1-p.δ) + y.Φ_1 )*s.k_1
    k_2 = ( (1-p.δ) + y.Φ_2 )*s.k_2

    (;k_1, k_2)  ## only the endogenous states

end


function arbitrage(model::hypeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)

    p = model.calibration
    (;β, γ, ψ, δ, θ) = p

    y = intermediate(model, s, x)
    Y = intermediate(model, S, X)

    # stochastic discount factors
    m_1 = β*( X.w_1/x.ew_1 )^(ψ-γ) * (Y.c_1/y.c_1)^(-ψ)
    m_2 = β*( X.w_2/x.ew_2 )^(ψ-γ) * (Y.c_2/y.c_2)^(-ψ)

    # return on capital
    r_1 = θ*(Y.y_1/S.k_1)*y.Φ_1__p + y.Φ_1__p/Y.Φ_1__p*( (1-δ) + Y.Φ_1 - X.i_1/S.k_1*Y.Φ_1__p )
    r_2 = θ*(Y.y_2/S.k_2)*y.Φ_2__p + y.Φ_2__p/Y.Φ_2__p*( (1-δ) + Y.Φ_2 - X.i_2/S.k_2*Y.Φ_2__p )

    res_1 = -(m_1*r_1 - 1)
    res_2 = -(m_2*r_2 - 1)
    res_3 = X.w_1^(1-γ) - x.ew_1^(1-γ)
    res_4 = X.w_2^(1-γ) - x.ew_2^(1-γ)
    res_5 = β*x.ew_1^(1-ψ) - x.w_1^(1-ψ) + (1-β)*y.c_1^(1-ψ)
    res_6 = β*x.ew_2^(1-ψ) - x.w_2^(1-ψ) + (1-β)*y.c_2^(1-ψ)

    return (;res_1, res_2, res_3, res_4, res_5, res_6)

end


# investment is nonnegative: 0 <= i_1, 0 <= i_2
function complementarities(model::hypeof(model), s::NamedTuple, x::NamedTuple, Fv::SVector)
    return (Fv[1] ⫫ x.i_1, Fv[2] ⫫ x.i_2, Fv[3], Fv[4], Fv[5], Fv[6])
end


model
