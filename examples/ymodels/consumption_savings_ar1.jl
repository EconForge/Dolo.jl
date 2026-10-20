using StaticArrays
import Dolo: transition, arbitrage, complementarities, hypeof
import Dolo: ⫫

# Consumption-savings with AR1 income and a borrowing constraint 0 <= c <= w
# (bruteforce/models/consumption_savings.yaml).

model = let

    name = :consumption_savings_ar1

    β = 0.96
    γ = 4.0
    r = 1.02
    ρ = 0.9
    σ_y = sqrt(0.001)

    y = 0.0
    w = 1.0
    c = 0.9*w

    calibration = (;β, γ, r, ρ, σ_y, y, w, c)

    states = Dolo.CartesianSpace(;
        :y => (-3σ_y/sqrt(1-ρ^2), 3σ_y/sqrt(1-ρ^2)),
        :w => (0.5, 20.0)
    )

    controls = Dolo.CartesianSpace(;
        :c => (0.0, Inf),
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
    eq = β*( X.c/x.c )^(-γ)*r - 1
    (eq,)
end


# 0 <= c <= w
function complementarities(model::hypeof(model), s::NamedTuple, x::NamedTuple, Fv::SVector)
    eq = Fv[1] ⫫ x.c
    eq = -eq ⫫ (s.w - x.c)
    eq = -eq
    (eq,)
end


model
