using Dolo

using StaticArrays
using LabelledArrays
using Dolo: SGrid, CGrid, PGrid, GArray, YModel
import Dolo: ×, ⟂, ⫫
import Dolo: transition, arbitrage, recalibrate, initial_guess, projection, equilibrium
import Dolo: X, reward
import Dolo: GridSpace, CartesianSpace

using Dolo: rouwenhorst

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
        y=[0.01, 100]
    )

    controls = CartesianSpace(;
        c=(0,Inf),
        λ=(0,Inf)
    )
    exogenous = Dolo.MarkovChain(
        (:w,:r,:e), P,Q
    )
    # grid 1 = SGrid(Q) × CGrid(((0.01,80.0,N),))
    
    name = :ayiagari

    YModel(
        name,
        states,
        controls,
        exogenous,
        calibration
    )

end

# small workaround to limit endless printing
# show(io::IO, tt::Type{typeof(model)}) = print(io, "DModel{#$(hash(tt))}")

@assert isbits(model)



# function recalibrate(mod::typeof(model); K=mod.calibration.y.K, N=mod.calibration.p.N)

#     # this lacks generality but should be future proof

#     (;m, s, x, y, z, p) = mod.calibration
    
#     (;ρ, σ, α, δ) = p

#     r = α*(1/K)^(1-α) - δ
#     w = (1-α)*K^α

#     m = merge(m, (;w=w, r=r))
#     y = merge(y, (;K=K))
#     p = merge(p, (;N=N))


#     N_ = convert(Int,N)
    
#     ## decide whether this should be matrix or smatrix
#     n_m = 3
#     mc = rouwenhorst(3,ρ,σ)
    
#     P = convert(SMatrix{n_m, n_m}, mc.p)
#     ## decide whether this should be matrix or smatrix
#     Q = SVector( (SVector(w,r,e) for e in mc.state_values)... )

#     grid = SSGrid(Q) × CGrid(((0.01,80.0,N_),))
    
#     nmodel = DModel(
#         (;m, s, x, y, z, p),
#         grid,
#         P
#     )

#     # important to check we are not changing the type
#     @assert (typeof(nmodel) == typeof(model))

#     return nmodel

# end


function transition(mod::typeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)
    p = model.calibration
    y = exp(M.e)*M.w + (s.y-x.c)*(1+M.r)
    return ( (;y) )
end


function arbitrage(mod::typeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)
    p = model.calibration
    eq = 1 - p.β*( X.c/x.c )^(-p.γ)*(1+S.r) - x.λ
    # @warn "The euler equation is satisfied only if c<w. If c=w, it can be strictly positive."
    # eq2 = x.λ ⟂ s.y-x.c
    eq2 = x.λ ⫫ s.y-x.c
    return ( (;eq, eq2) )
end

# function initial_guess(model, m::SLArray, s::SLArray, p)
#     # c = min( 1.0 + 0.01*(s.y - 1.0), s.y)
#     c = 0.8*s.y
#     λ = 0.01 # max( 1.0 + 0.01*(s.y - 1.0), 0.01*(s.y-1))
#     return SLVector(;c, λ)
# end

function initial_guess(model, s::NamedTuple)
    p = model.calibration
    # c = min( 1.0 + 0.01*(s.y - 1.0), s.y)
    # c = exp(m.e)*m.w *0.8
    c = s.y*0.9
    return (;c)
end

function reward(model::typeof(model), s::NamedTuple, x::NamedTuple)
    p = model.calibration
    c = x.c
    return c^(1-p.γ) / (1-p.γ)
end

function X(model, s::NamedTuple)

    p = model.calibration.p
    y = s.y
    ([0.0001], [y])

end

model