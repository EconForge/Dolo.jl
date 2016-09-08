immutable TaylorExpansion{order} <: AbstractDecisionRule
    s0::Vector{Float64}
    x0::Vector{Float64}
    x_1::Array{Float64,2}
    x_2::Array{Float64,3}
    x_3::Array{Float64,4}
end

TaylorExpansion(s0, x0, x_1) =
    TaylorExpansion{1}(s0, x0, x_1, zeros(0, 0, 0), zeros(0, 0, 0, 0))

TaylorExpansion(s0, x0, x_1, x_2) =
    TaylorExpansion{2}(s0, x0, x_1, x_2, zeros(0, 0, 0, 0))

TaylorExpansion(s0, x0, x_1, x_2, x_3) =
    TaylorExpansion{3}(s0, x0, x_1, x_2, x_3)

function _check_call(ts::TaylorExpansion, points::AbstractVector,
                     out::AbstractVector)
    ns = length(ts.s0)
    nx = length(ts.x0)

    if size(points, 1) != ns
        msg = "Points should be of length $(ns), found $(size(points, 1))"
        throw(DimensionMismatch(msg))
    end

    if size(out, 1) != nx
        msg = "out should be of length $(nx), found $(size(out, 1))"
        throw(DimensionMismatch(msg))
    end

    if size(out, 2) != size(points, 2)
        msg = "Out should have same number of columns as points"
        throw(DimensionMismatch(msg))
    end
end

# points should be a vector of observations of all state variables at one time
# period
@compat function (ts::TaylorExpansion{1})(points::AbstractVector,
                                          out::AbstractVector=similar(ts.x0),
                                          chk::Bool=true)
    chk && _check_call(ts, points, out)
    s0, x0, x_1 = ts.s0, ts.x0, ts.x_1
    ns = length(s0)
    nx = length(x0)

    copy!(out, x0)

    @inbounds for i in 1:ns
        @simd for n in 1:nx
            out[n] += x_1[n, i]*(points[i] - s0[i])
        end
    end
    out
end

@compat function (ts::TaylorExpansion{2})(points::AbstractVector,
                                          out::AbstractVector=similar(ts.x0),
                                          chk::Bool=true)
    chk && _check_call(ts, points, out)
    s0, x0, x_1, x_2 = ts.s0, ts.x0, ts.x_1, ts.x_2
    ns = length(s0)
    nx = length(x0)

    copy!(out, x0)

    @inbounds for j in 1:ns
        foo_j = points[j] - s0[j]
        for i in 1:ns
            foo_i = points[i] - s0[i]
            @simd for n in 1:nx
                out[n] += x_1[n, i]*foo_i
                out[n] += x_2[n, i, j]*foo_i*(foo_j)/2.0
            end
        end
    end
    out
end

@compat function (ts::TaylorExpansion{3})(points::AbstractVector,
                                          out::AbstractVector=similar(ts.x0),
                                          chk::Bool=true)
    chk && _check_call(ts, points, out)
    s0, x0, x_1, x_2, x_3 = ts.s0, ts.x0, ts.x_1, ts.x_2, ts.x_3
    ns = size(points, 1)
    nx = size(x0, 1)

    copy!(out, x0)

    @inbounds for k in 1:ns
        Δk = points[k] - s0[k]
        for j in 1:ns
            Δj = points[j] - s0[j]
            for i in 1:ns
                Δi = points[i] - s0[i]
                @simd for n in 1:nx
                    out[n] += x_1[n, i]*Δi
                    out[n] += x_2[n, i, j]*Δi*(Δj)/2.0
                    out[n] += x_3[n, i, j, k]*(Δi)*(Δj)*(Δk)/6.0
                end
            end
        end
    end
    out
end

# Each row is an observation of all the state variables
@compat function (ts::TaylorExpansion)(points::AbstractMatrix,
                                       out::AbstractMatrix=Array(Float64, size(points, 1), length(ts.x0)))
    ns = length(ts.s0)
    nx = length(ts.x0)

    if size(points, 2) != ns
        msg = "Points should have $(ns) columns, found $(size(points, 2))"
        throw(DimensionMismatch(msg))
    end

    if size(out, 2) != nx
        msg = "out should have $(nx) columns, found $(size(out, 2))"
        throw(DimensionMismatch(msg))
    end

    if size(out, 1) != size(points, 1)
        msg = "Out should have same number of rows as points"
        throw(DimensionMismatch(msg))
    end

    for n in 1:size(points, 1)
        ts(view(points, n, :), view(out, n, :), false)
    end
    out
end
