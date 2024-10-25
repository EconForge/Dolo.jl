using Dolo.Build
import Dolo.Build: *

model = let 

    β = 0.96
    γ = 4.0
    σ = 0.1
    ρ = 0.8
    # r = 1.025
    # y = 1.0 # available income
    
    K = 40.0
    α = 0.36
    A = 1
    δ = 0.025

    r = α*(1/K)^(1-α) - δ
    w = (1-α)*K^α

    y = w

    c = 0.9*y

    λ = 0.1

    e = 0
    cbar = c

    N = 200


    m = (;w,r,e)
    s = (;y)
    x = (;c, λ)
    y = (;K)
    z = (;z=0.0)
    p = (;β, γ, σ, ρ, cbar, α, δ, N)

    calibration = merge(m,s,x,y,z,p)


    n_m = 3
    mc = rouwenhorst(3,ρ,σ)
    
    P = convert(SMatrix{n_m, n_m}, mc.P)
    ## decide whether this should be matrix or smatrix
    Q = SVector( (SVector(w,r,e) for e in mc.V)... )

    states = GridSpace(
        (:w,:r,:e),
        SVector( [Q[i] for i=1:size(Q,1)]...  )
    )×CartesianSpace(;
        y=(0.01, 100.0)
    )

    controls = CartesianSpace(;
        c=(0.0,Inf),
    )
    exogenous = Dolo.MarkovChain(
        (:w,:r,:e), P,Q
    )
    
    name = :ayiagari

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
    y = exp(M.e)*M.w + (s.y-x.c)*(1+M.r)
    return ( (;y) )
end


function arbitrage(mod::typeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)
    p = mod.calibration
    eq = 1 - p.β*( X.c/x.c )^(-p.γ)*(1+S.r) # - x.λ
    # @warn "The euler equation is satisfied only if c<w. If c=w, it can be strictly positive."
    # eq2 = x.λ ⟂ s.y-x.c
    (eq,)
end

function complementarities(mod::typeof(model), s::NamedTuple, x::NamedTuple, Fv::SVector)
    # eq = Fv[1] ⫫ s.y-x.c
    eq = Fv[1] ⟂ s.y-x.c
    (eq,)
end

# function initial_guess(model, m::SLArray, s::SLArray, p)
#     # c = min( 1.0 + 0.01*(s.y - 1.0), s.y)
#     c = 0.8*s.y
#     λ = 0.01 # max( 1.0 + 0.01*(s.y - 1.0), 0.01*(s.y-1))
#     return SLVector(;c, λ)
# end

function initial_guess(mod::typeof(model), s::NamedTuple)
    p = mod.calibration
    c = s.y*0.9
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