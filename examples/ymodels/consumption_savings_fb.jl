using StaticArrays
import Dolo: transition, arbitrage, initial_guess, hypeof
import Dolo: ⫫

# Consumption-savings with AR1 income and a borrowing constraint 0 <= c <= w.
# The constraint is written explicitly: h is its multiplier, and h ⫫ (w-c)
# (Fischer-Burmeister) is one of the equations, instead of using complementarities.
# Model "CS" of the improved time iteration paper (bruteforce/models/consumption_savings_new.yaml).

model = let

    name = :consumption_savings_fb

    β = 0.96
    γ = 4.0
    r = 1.02
    ρ = 0.9
    σ_y = sqrt(0.001)

    y = 0.0
    w = 1.0
    c = 0.9*w
    h = 0.0

    calibration = (;β, γ, r, ρ, σ_y, y, w, c, h)

    states = Dolo.CartesianSpace(;
        :y => (-3σ_y/sqrt(1-ρ^2), 3σ_y/sqrt(1-ρ^2)),
        :w => (0.5, 4.0)
    )

    controls = Dolo.CartesianSpace(;
        :c => (0.0, Inf),
        :h => (0.0, Inf)
    )

    Σ = @SMatrix [σ_y^2;;]
    process = Dolo.VAR1( (:y,), ρ, Σ )

    Dolo.YModel(name, states, controls, process, calibration)

end


function transition(model::hypeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)
    (;r) = model.calibration
    w = exp(M.y) + (s.w - x.c)*r
    (;w)  ## only the endogenous state
end


function arbitrage(model::hypeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)
    (;β, γ, r) = model.calibration
    res_1 = x.h ⫫ (s.w - x.c)
    res_2 = β*( X.c/x.c )^(-γ)*r - 1 + x.h
    (;res_1, res_2)
end


function initial_guess(model::hypeof(model), s::NamedTuple)
    T = typeof(s.w)   # precision of the model
    c = min(s.w*T(0.95), 1 + T(0.04)*(s.w - 1))
    h = T(0.05)
    (;c, h)
end


model
