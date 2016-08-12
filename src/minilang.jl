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

# ------------------------- #
# Discrete Transition types #
# ------------------------- #

function _build_discrete_transition(data::Associative, calib::ModelCalibration)
    if data[:tag] == :MarkovChain
        # need to extract/clean up P and Q
        state_values = map(Vector{Float64}, eval_with(calib, data[:P]))
        n = length(state_values)
        Π = Array(Float64, n, n)
        Q = eval_with(calib, data[:Q])
        for i in 1:n
            Π[i, :] = Q[i]
        end
        return MarkovChain(Π, state_values)
    else
        m = "don't know how to handle discrete_transition of type $(data[:tag])"
        error(m)
    end
end
