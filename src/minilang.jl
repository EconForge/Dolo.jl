# ---------- #
# Grid types #
# ---------- #

abstract AbstractGrid

immutable Cartesian <: AbstractGrid
    a::Vector{Float64}
    b::Vector{Float64}
    orders::Vector{Int}
end

function Cartesian(stuff::Associative)
    kind = get(stuff, :tag, nothing)
    if kind != :Cartesian
        error("Can't build Cartesian from dict with kind $(kind)")
    end

    Cartesian(stuff[:a], stuff[:b], stuff[:orders])
end

Cartesian(;a=[], b=[], orders=[]) = Cartesian(a, b, orders)

function _build_grid(data::Associative)
    if data[:tag] == :Cartesian
        return Cartesian(data)
    else
        m = "don't know how to handle grid of type $(data[:tag])"
        error(m)
    end
end

# ------------------ #
# Distribution types #
# ------------------ #

abstract AbstractDistribution

immutable Normal <: AbstractDistribution
    sigma::Matrix{Float64}
end

Normal(;sigma=zeros(0, 0)) = Normal(sigma)

function _build_dist(data::Associative, calib::ModelCalibration)
    if data[:tag] == :Normal
        n = length(calib[:shocks])
        sigma = reshape(vcat(data[:sigma]...), n, n)
        return Normal(_to_Float64(sigma))
    else
        m = "don't know how to handle distribution of type $(data[:tag])"
        error(m)
    end
end
