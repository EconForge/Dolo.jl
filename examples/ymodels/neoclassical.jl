
using Dolo.Build


## Define Model

model = let 

    α = 0.3
    β = 0.96
    γ = 2.0
    δ = 0.1

    k = ((1/β - (1-δ))/α)^(1/(α-1))
    i = δ*k
    z = 0.0

    p = (;α, β, γ, δ)
    
    m = ( (; z))
    s = ( (; k))
    x = ( (; i))

    calibration = merge(m,s,x,p)

    P = @SMatrix [0.9 0.1; 0.1 0.9]
    # Q = @SMatrix [-0.05; 0.05]
    Q = @SMatrix [-0.1; 0.1]
    exogenous = MarkovChain(
        (:z,),
        Matrix(P),
        Matrix(Q)

    )

    states = GridSpace(
        (:z,),
        SVector( [Q[i,:] for i=1:size(Q,1)]... )
    ) × CartesianSpace(
        k=[0.1, 5.0]
    )

    controls = CartesianSpace(
        i=(0.0, 10.0)
    )

    # exo = SGrid( [Q[i,:] for i=1:size(Q,1)] )
    # endo = CGrid( ((0.1, 5.0, 500),) )
    # grid = exo × endo

    YModel(
        :neoclassical,
        states,
        controls,
        exogenous,
        calibration
    )

end

function transition(model::typeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)
    p = model.calibration
    K = s.k*(1-p.δ) + x.i
    return ( (;k=K) )
end

function arbitrage(model::typeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)
    p = model.calibration
    c = exp(s.z)*s.k^p.α - x.i
    C = exp(S.z)*S.k^p.α - X.i
    r = p.β*(C/c)^(-p.γ)*(1-p.δ + p.α*exp(S.z)*S.k^(p.α-1)) - 1
    return ( (;r) )
end

function X(model, s)

    p = model.calibration.p

    y = exp(s.z)*s.k^p.α 
    ([0], [y*0.95])

end

function reward(model::typeof(model), s::NamedTuple, x::NamedTuple, p)
    c = exp(s.z)*s.k^p.α - x.i
    return c^(1-p.γ) / (1-p.γ)
end

model