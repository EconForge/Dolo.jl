using StaticArrays
import Dolo: transition, arbitrage
import Dolo

model = let 

    name = :rbc_mc

    # calibrate some parameters
    β = 0.99
    σ = 5.0
    η = 1.0
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


dmodel = Dolo.discretize(model)


model32 = Dolo.convert_precision(Float32,model)
dmodel32 = Dolo.discretize(model32)

wksp = Dolo.time_iteration_workspace(dmodel32)

itps = wksp.φ.itp

(;x0,r0, φ)  = wksp

Dolo.F(dmodel32, x0, φ)





Dolo.time_iteration(model)