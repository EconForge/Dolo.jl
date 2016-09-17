# ---------- #
# Grid types #
# ---------- #

abstract AbstractGrid

immutable Cartesian <: AbstractGrid
    a::Vector{Float64}
    b::Vector{Float64}
    orders::Vector{Int}
end

function Cartesian(stuff::Associative, calib::ModelCalibration)
    kind = get(stuff, :tag, nothing)
    if kind != :Cartesian
        error("Can't build Cartesian from dict with kind $(kind)")
    end

    a = eval_with(calib, stuff[:a])
    b = eval_with(calib, stuff[:b])
    orders = eval_with(calib, stuff[:orders])
    Cartesian(a, b, orders)
end

Cartesian(;a=[], b=[], orders=[]) = Cartesian(a, b, orders)

function _build_grid(data::Associative, calib::ModelCalibration)
    if data[:tag] == :Cartesian
        return Cartesian(data, calib)
    else
        m = "don't know how to handle grid of type $(data[:tag])"
        error(m)
    end
end

# ------------------ #
# Distribution types #
# ------------------ #

abstract AbstractDistribution

function _build_dist(data::Associative, calib::ModelCalibration)
    if data[:tag] == :Normal
        n = length(calib[:shocks])
        sigma = reshape(vcat(data[:sigma]...), n, n)
        return MvNormal(_to_Float64(sigma))
    else
        m = "don't know how to handle distribution of type $(data[:tag])"
        error(m)
    end
end

abstract AbstractExogenous
abstract DiscreteExogenous <: AbstractExogenous
abstract ContinuousExogenous <: AbstractExogenous
abstract IIDExogenous <: AbstractExogenous

# TOOD: generalize to non-scalar...

immutable Normal <: IIDExogenous
    sigma::Matrix{Float64}
end

Normal(sigma::Float64) = Normal(reshape([sigma],1,1))

immutable VAR1 <: ContinuousExogenous
    rho::Matrix{Float64}
    sigma::Matrix{Float64}
    N::Vector{Int}
end

AR1(rho,sigma,N) = VAR1(rho,sigma,N)
AR1(rho::Float64, sigma::Float64, N::Int) = VAR1(reshape([rho],1,1),reshape([sigma],1,1),reshape([N],1))

immutable MarkovChain{T} <: DiscreteExogenous
    transitions::Matrix{T}
    values::Matrix{T}
end

# ------------------------- #
# Discrete Transition types #
# ------------------------- #


function _build_exogenous_entry(data::Associative, calib::ModelCalibration)
    if data[:tag] == :MarkovChain
        # need to extract/clean up P and Q

        P = eval_with(calib, data[:P])
        n = length(P)
        state_values = Array(Float64, n, n)
        for i in 1:n
            state_values[i, :] = P[i]
        end

        # n = length(state_values)
        Q = eval_with(calib, data[:Q])
        Π = Array(Float64, n, n)
        for i in 1:n
            Π[i, :] = Q[i]
        end
        return MarkovChain(Π, state_values)
    elseif data[:tag] == :AR1
        # need to extract rho an dsigma
        rho = eval_with(calib, data[:rho])
        sigma = eval_with(calib, data[:sigma])
        N = eval_with(calib, get(data, :N, 10))  # TODO: should default be 10??
        return AR1(rho, sigma, N)
    elseif data[:tag] == :Normal
        # need to extract rho an dsigma
        sigma = eval_with(calib, data[:sigma])
        return Normal(sigma)
    end
    m = "don't know how to handle exogenous process of type $(data[:tag])"
    error(m)
end
