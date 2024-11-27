using Dolo.Build
import Dolo.Build: *

using Dolo: SGrid, CGrid, PGrid, GArray, YModel
import Dolo: ×, ⟂, ⫫
import Dolo: transition, arbitrage, recalibrate, initial_guess, projection, equilibrium, complementarities
import Dolo: X, reward
import Dolo: GridSpace, CartesianSpace

using Dolo: rouwenhorst


using StaticArrays
model = let 

    β = 0.96
    γ = 4.0
    σ = 0.1
    ρ = 0.0
    r = 1.02
    
    w = 1.0
    c = 0.9*w
    cbar = c

    ξ = 0.0
    y = 0.0

    a = 0
    mr = 0

    rew = 0.0


    m = (;y)
    s = (;w)
    x = (;c)
    p = (;β, γ, σ, ρ, r, cbar)

    calibration = merge(p,m,s,x)


    states = CartesianSpace(;
        w=(0.5, 5.0)
    )

    controls = CartesianSpace(;
        c=(0.0,Inf),
    )

    Σ = @SMatrix [0.01;;]
    exogenous = MvNormal( (:y,), Σ)


    name = :consumption_savings

    YModel(
        name,
        states,
        controls,
        exogenous,
        calibration
    )

end

@assert isbits(model)


function transition(mod::typeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)
    p = mod.calibration
    w = exp(M.y) + (s.w-x.c)*(p.r)
    return ( (;w) )
end


function arbitrage(mod::typeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)
    p = mod.calibration
    eq = p.β*( X.c/x.c )^(-p.γ)*(p.r) - 1 # - x.λ
    (eq,)
end

function complementarities(mod::typeof(model), s::NamedTuple, x::NamedTuple, Fv::SVector)
    eq = Fv[1] ⫫ x.c
    eq = -eq ⫫ s.w-x.c
    eq = -eq
    (eq,)
end

function initial_guess(mod::typeof(model), s::NamedTuple)
    p = mod.calibration
    # c = p.c
    c = s.w*0.95
    # c = min(s.w, 1.0 + (s.w-1.0)*0.02)
    return (;c)
end

function reward(mod::typeof(model), s::NamedTuple, x::NamedTuple)
    p = mod.calibration
    c = x.c
    return c^(1-p.γ) / (1-p.γ)
end

# function X(mod::typeof(model), s::NamedTuple)

#     p = mod.calibration.p
#     y = s.y
#     ([0.0001], [y])

# end

model