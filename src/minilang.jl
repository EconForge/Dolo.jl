# ---------- #
# Grid types #
# ---------- #

@compat abstract type AbstractGrid end

immutable Cartesian <: AbstractGrid
    a::Vector{Float64}
    b::Vector{Float64}
    orders::Vector{Int}
end

type Domain
    states::Vector{Symbol}
    min::Vector{Float64}
    max::Vector{Float64}
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

@compat abstract type AbstractDistribution end

function _build_dist(data::Associative, calib::ModelCalibration)
    if data[:tag] == :Normal
        n = length(calib[:shocks])
        Sigma = reshape(vcat(data[:Sigma]...), n, n)
        return Distributions.MvNormal(_to_Float64(Sigma))
    else
        m = "don't know how to handle distribution of type $(data[:tag])"
        error(m)
    end
end

# ------------------------- #
# Discrete Transition types #
# ------------------------- #

to_vector(tab::Number) = [tab]
to_matrix(tab::Number) = reshape([tab], 1, 1)
to_matrix(tab::Array) = hcat([Array{Float64}(e) for e in tab]...)
to_matrix(tab::Array{Array{Float64,1},1}) = cat(1, [e' for e in tab]...)

function _build_exogenous_entry(data::Associative, calib::ModelCalibration)

    if data[:tag] == :Product
        p1 = _build_exogenous_entry(data[:p1], calib)
        p2 = _build_exogenous_entry(data[:p2], calib)
        return ProductProcess(p1,p2)
    elseif data[:tag] == :MarkovChain
        # need to extract/clean up P and Q
        values = eval_with(calib, data[:values])
        states_values = to_matrix(values)
        transitions = eval_with(calib, data[:transitions])
        Π = to_matrix(transitions)
        return MarkovChain(Π, states_values)
    elseif data[:tag] == :VAR1
        # need to extract rho an dSigma
        rho = eval_with(calib, data[:rho])
        Sigma = eval_with(calib, data[:Sigma])
        N = eval_with(calib, get(data, :N, 10))  # TODO: should default be 10??
        # rho = to_matrix(rho)
        Sigma = to_matrix(Sigma)
        N = to_vector(N)
        return VAR1(rho, Sigma)
    elseif data[:tag] == :Normal
        # need to extract rho an dSigma
        Sigma = eval_with(calib, data[:Sigma])
        Sigma = to_matrix(Sigma)
        return Normal(Sigma)
    end
    m = "don't know how to handle exogenous process of type $(data[:tag])"
    error(m)
end
