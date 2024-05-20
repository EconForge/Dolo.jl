
using StaticArrays
import Dolo: transition, arbitrage
import Dolo

model = let 

    name = :rbc_mc

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

    calibration = (;β, σ, η, δ, α, ρ, z, n, k, w, y, i, c, χ)


    P = @SMatrix [0.4 0.6; 0.6 0.4]
	# Q = @SMatrix [-0.01; 0.01]
    Q = SVector( SVector(-0.01), SVector(0.01) )

    # P = @SMatrix [1.0;]
    # Q = SVector( (SVector(0.0),) )
        
    process = Dolo.MarkovChain( (:z,), P, Q )

    states = Dolo.ProductSpace(
        Dolo.GridSpace((:z,), Q),
        Dolo.CartesianSpace(;
            :k => ( k*0.9, k*1.1),
        )
    )
    controls = Dolo.CartesianSpace(;
        :i => (0.0, 10.0),
        :n => (0.0, 1.5)
    )

    # calibration = (;α, β, γ, δ, ρ,χ, η = 2.0, σ = 2.0, i=0.1, n=0.8)


    Dolo.YModel(name, states, controls, process, calibration)

end


using LoopVectorization
dmodel = Dolo.discretize(model, Dict(:endo=>10000000))

(;φ, x0) = Dolo.time_iteration_workspace(dmodel)

function FF!(r0, dm, x0,  φ)

    MM = length(dm.grid.g1)
    a,b,NN = dm.grid.g2.ranges[1]
    g1 = dm.grid.g1
    g2 = dm.grid.g2

    p = dm.model.calibration
    
    # (;β, σ, η, δ, α, ρ, χ) = dm.model.calibration

    P = dm.model.exogenous.P

    itps = φ.itp

    for i_m=1:MM

        z = g1.points[i_m][1]

        for n_=1:NN

            i,n = x0[i_m, n_]

            k = g2[n_][1] # a+(b-a)*(n_-1)/(N-1)

            y = exp(z)*(k^p.α)*(n^(1-p.α))
            w = (1-p.α)*y/n

            c = y - i

            res1 = 0.0
            res2 = 0.0

            for i_M=1:MM

                prob = P[i_m,i_M]
                Z = g1.points[i_M][1]
                
                
                itp = itps[i_M]
                
                
                K = k*(1-p.δ) + i

                I,N = itp(K)
                Y = exp(Z)*(K^p.α)*(N^(1-p.α))
                C = Y - I

                res1 += prob*(p.χ*(n^p.η)*(c^p.σ) - w)
                res2 += prob*(p.β*(C/c)^(-p.σ)*(1-p.δ+exp(Z)*p.α*K^(p.α-1)*N^(1-p.α))-1)

            end

            r0[i_m,n_] = SVector(res1,res2)

        end
    end

end


function Dolo.transition(model::typeof(model), s::NamedTuple, x::NamedTuple, M::NamedTuple)
    
    (;δ, ρ) = model.calibration
    
    # Z = e.Z
    K = s.k * (1-δ) + x.i

    (;k=K,)  ## This is only the endogenous state

end



function intermediate(model::typeof(model),s::NamedTuple, x::NamedTuple)
    
    p = model.calibration

	y = exp(s.z)*(s.k^p.α)*(x.n^(1-p.α))
	w = (1-p.α)*y/x.n
	rk = p.α*y/s.k
	c = y - x.i
	return ( (; y, c, rk, w))

end


function arbitrage(model::typeof(model), s::NamedTuple, x::NamedTuple, S::NamedTuple, X::NamedTuple)

    p = model.calibration

	y = intermediate(model, s, x)
	Y = intermediate(model, S, X)
	res_1 = p.χ*(x.n^p.η)*(y.c^p.σ) - y.w
	res_2 = (p.β*(y.c/Y.c)^p.σ)*(1 - p.δ + Y.rk) - 1
    
    return ( (;res_1, res_2) )

end


r0 = deepcopy(x0*0)
r1 = deepcopy(x0*0)


using BenchmarkTools


# not faster !
@time FF!(r0, dmodel, x0, φ)
@time Dolo.F!(r1, dmodel, x0, φ)

@benchmark FF!(r0, dmodel, x0, φ)

@benchmark Dolo.F!(r1, dmodel, x0, φ)