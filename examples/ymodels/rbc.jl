using StaticArrays
import Dolo: transition, arbitrage

# RBC model with AR1 productivity driven by an iid normal shock.
# Consumption, wage and return on capital are treated as controls, defined by equations.
# Model "RBC" of the improved time iteration paper (bruteforce/models/rbc.yaml).

model = let

    name = :rbc

    β = 0.99
    σ = 2.0
    η = 1.0
    δ = 0.025
    α = 0.33
    ρ = 0.8
    σ_z = 0.01

    # steady-state
    a = 0.0
    n = 0.33
    rk = 1/β - 1 + δ
    k = n/(rk/α)^(1/(1-α))
    w = (1-α)*exp(a)*(k/n)^α
    i = δ*k
    c = exp(a)*k^α*n^(1-α) - i
    χ = w/c^σ/n^η

    calibration = (;β, σ, η, χ, δ, α, ρ, σ_z, k, a, i, n, c, w, rk)

    states = Dolo.CartesianSpace(;
        :k => (k*0.5, k*2.0),
        :a => (-0.1, 0.1)
    )

    controls = Dolo.CartesianSpace(;
        :i => (-Inf, Inf),
        :n => (-Inf, Inf),
        :c => (0.0, Inf),
        :w => (0.0, Inf),
        :rk => (-Inf, Inf)
    )

    Σ = @SMatrix [σ_z^2;;]
    process = Dolo.MvNormal( (:z,), Σ )

    Dolo.YModel(name, states, controls, process, calibration)

end


function transition(model::typeof(model), s::NamedTuple, x::NamedTuple, e::NamedTuple)

    (;δ, ρ) = model.calibration

    k = (1-δ)*s.k + x.i
    a = ρ*s.a + e.z

    (;k, a)

end


function arbitrage(model::typeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)

    (;β, σ, η, χ, δ, α) = model.calibration

    res_1 = 1 - β*(x.c/X.c)^σ*(1-δ+X.rk)
    res_2 = χ*x.n^η*x.c^σ - x.w
    res_3 = exp(s.a)*s.k^α*x.n^(1-α) - x.i - x.c
    res_4 = α*exp(s.a)*(x.n/s.k)^(1-α) - x.rk
    res_5 = (1-α)*exp(s.a)*(s.k/x.n)^α - x.w

    return (;res_1, res_2, res_3, res_4, res_5)

end


model
