using StaticArrays
import Dolo: arbitrage, transition

model = let 

    name = :rbc_iid

    # calibrate some parameters
    β = 0.99
    σ = 5
    η = 1
    δ = 0.025
    α = 0.33
    ρ = 0.8
    zbar = 0.0
    σ_z = 0.016
    n = 0.33
    z = zbar
    rk = 1/β - 1+δ
    k = n/(rk/α)^(1/(1-α))
    w = (1-α)*exp(z)*(k/n)^α
    y = exp(z)*k^α*n^(1-α)
    i = δ*k
    c = y - i
    χ =  w/c^σ/n^η

    sig_z = 0.016

    calibration = (;β, σ, η, δ, α, ρ, z, n, k, w, y, i, c, χ)


    states = Dolo.CartesianSpace(;
        :z => (-2*sig_z/(1-ρ^2)^0.5,  2*sig_z/(1-ρ^2)^0.5),
        :k => ( k*0.5,  k*1.5)
    )

    controls = Dolo.CartesianSpace(;
        :i => (0.0, 10.0),
        :n => (0.0, 1.5)        
    )
    
    Σ = @SMatrix [sig_z^2;]
    process = Dolo.MvNormal( (:ϵ,), Σ )

    Dolo.YModel(name, states, controls, process, calibration)

end;

function transition(model::typeof(model), s::NamedTuple, x::NamedTuple, e::NamedTuple)
    
    (;δ, ρ) = model.calibration
    
    Z = s.z * ρ + e.ϵ
    K = s.k * (1-δ) + x.i

    (;Z, K)

end


function intermediate(model::typeof(model),s::NamedTuple, x::NamedTuple)
    
    p = model.calibration

	y = exp(s.z)*(s.k^p.α)*(x.n^(1-p.α))
	w = (1-p.α)*y/x.n
	rk = p.α*y/s.k
	c = y - x.i
	return ( (; y, c, rk, w))

end


function arbitrage(model::typeof(model), s::NamedTuple, x::NamedTuple,  S::NamedTuple, X::NamedTuple)

    p = model.calibration

	y = intermediate(model, s, x)
	Y = intermediate(model, S, X)
	res_1 = p.χ*(x.n^p.η)*(y.c^p.σ) - y.w
	res_2 = (p.β*(y.c/Y.c)^p.σ)*(1 - p.δ + Y.rk) - 1
    
    return ( (;res_1, res_2) )

end

model