using StaticArrays
import Dolo: transition, arbitrage, complementarities
import Dolo: ⫫

# New Keynesian model with Calvo pricing and price dispersion (Delta),
# productivity and cost-push shocks following a Markov chain.
# (bruteforce/models/NK_dsge.yaml)

model = let

    name = :nk_dsge

    β = 0.995
    ψ = 1.0
    ε = 11.0
    ϕ = 1.5
    θ = 0.4
    χ = 1.0
    τ = 1/ε
    p = 0.5

    # steady-state
    a = 0.0
    cp = 0.0
    π = 1.0
    n = 1.0
    c = 1.0
    rw = χ*n^ψ*c
    mc = (1-τ)*rw/exp(a)
    F1 = (ε/(ε-1))*mc/(1-θ*β*π^ε)
    F2 = 1/(1-θ*β*π^(ε-1))
    Delta__L = 1.0

    calibration = (;β, ψ, ε, ϕ, θ, χ, τ, p, a, cp, Delta__L, c, n, F1, F2, π)

    P = @SMatrix [p p/2 p/2; p/2 0.5 p/2; p/2 p/2 0.5]
    Q = SVector( SVector(0.0, 0.0), SVector(-0.005, 0.0), SVector(0.005, 0.0) )

    process = Dolo.MarkovChain( (:a, :cp), P, Q )

    states = Dolo.ProductSpace(
        Dolo.GridSpace((:a, :cp), Q),
        Dolo.CartesianSpace(;
            :Delta__L => (1.0, 1.05),
        )
    )

    controls = Dolo.CartesianSpace(;
        :c => (0.0, Inf),
        :n => (0.0, Inf),
        :F1 => (-Inf, Inf),
        :F2 => (-Inf, Inf),
        :π => (1-0.000001, 2.0)
    )

    Dolo.YModel(name, states, controls, process, calibration)

end


function intermediate(model::typeof(model), s::NamedTuple, x::NamedTuple)
    p = model.calibration
    r = (1/p.β)*x.π^p.ϕ
    rw = p.χ*x.n^p.ψ*x.c
    mc = (1-p.τ*exp(s.cp))*rw/exp(s.a)
    Delta = p.θ*x.π^p.ε*s.Delta__L + (1-p.θ)*(x.F1/x.F2)^(-p.ε)
    (;r, rw, mc, Delta)
end


function transition(model::typeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)
    y = intermediate(model, s, x)
    (;Delta__L = y.Delta)  ## only the endogenous state
end


function arbitrage(model::typeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)
    (;β, ε, θ) = model.calibration
    y = intermediate(model, s, x)
    res_1 = 1/x.c - β*y.r/(X.c*X.π)
    res_2 = (ε/(ε-1))*y.mc + θ*β*X.π^ε*X.F1 - x.F1
    res_3 = 1 + θ*β*X.π^(ε-1)*X.F2 - x.F2
    res_4 = x.c - exp(s.a)/y.Delta*x.n
    res_5 = x.π - ((1-(x.F1/x.F2)^(1-ε)*(1-θ))/θ)^(1/(ε-1))
    (;res_1, res_2, res_3, res_4, res_5)
end


# 1-1e-6 <= π <= 2
function complementarities(model::typeof(model), s::NamedTuple, x::NamedTuple, Fv::SVector)
    eq = Fv[5] ⫫ (x.π - (1-0.000001))
    eq = -eq ⫫ (2.0 - x.π)
    eq = -eq
    (Fv[1], Fv[2], Fv[3], Fv[4], eq)
end


model
