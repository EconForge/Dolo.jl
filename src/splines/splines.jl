module splines

export filter_coeffs, interpolant_cspline, filter_coeffs
export eval_UC_spline, eval_UC_spline!
export prefilter!

include("csplines.jl")
include("splines_filter.jl")
include("cubic_prefilter.jl")
include("interp.jl")


function interpolant_cspline(a, b, orders, V)

    coefs = filter_coeffs(a, b, orders, V)

    function fun(s::Array{Tf,2}) where Tf
        return eval_UC_spline(a, b, orders, coefs, s)
    end

    function fun(p::Tf...) where Tf
        return fun([p...]')
    end

    return fun

end

abstract type Linear end
abstract type Cubic end

struct MLinear<:Linear end
struct MCubic<:Cubic end

struct CubicInterpolator{G,C} <: Cubic
    grid::G
    θ::C
end

struct SplineInterpolator{G,C,k}
    grid::G
    θ::C
end




    function prefilter(ranges::NTuple{d,Tuple{Tf,Tf,i}}, V::AbstractArray{T, d}, ::Val{3}) where d where i<:Int where T where Tf
        θ = zeros(eltype(V), ((e[3]+2) for e in ranges)...)
        ind = tuple( (2:(e[3]+1) for e in ranges )...)
        θ[ind...] = V
        splines.prefilter!(θ)
        return θ
    end


    function prefilter!(θ::AbstractArray{T, d}, grid::NTuple{d,Tuple{Tf,Tf,i}}, V::AbstractArray{T, d}, ::Val{3}) where d where i<:Int where T where Tf
        splines.prefilter!(θ)
    end

    function prefilter(ranges::NTuple{d,Tuple{Tf,Tf,i}}, V::AbstractArray{T, d}, ::Val{1}) where d where i<:Int where T where Tf
        θ = copy(V)
        return θ
    end


    function prefilter!(θ::AbstractArray{T, d}, grid::NTuple{d,Tuple{Tf,Tf,i}}, V::AbstractArray{T, d}, ::Val{1}) where d where i<:Int where T where Tf
        θ .= V
    end

    # function CubicInterpolator(grid; values=nothing)

    #     println("Hejo")

    #     n = [e[3] for e in grid.ranges]
    #     θ = zeros(eltype(values), (i+2 for i in n)...)
    #     ci = CubicInterpolator{typeof(grid), typeof(θ)}(grid, θ)
    #     if !isnothing(values)
    #         ind = tuple( (2:(e[3]+1) for e in grid.ranges )...)
    #         ci.θ[ind...] .= values
    #         splines.prefilter!(ci.θ)
    #     end
    #     return ci
    
    # end


    function fit!(spl::SplineInterpolator{G,C,k}, values) where G where C where k

        # grid = spl.grid
        # n = [e[3] for e in spl.grid]
        nn = tuple( ((e[3]) for e in spl.grid )...)
        rhs = reshape(view(values,:),nn...)
        if k==3
            fill!(spl.θ, zero(eltype(spl.θ)))
            ind = tuple( (2:(e[3]+1) for e in spl.grid )...)
            spl.θ[ind...] .= rhs
            prefilter!(spl.θ, Val(:KA))
        elseif k==1
            spl.θ .= rhs
        end
    
    end

    function SplineInterpolator(ranges::Tuple; values=nothing, k=3)


        n = [e[3] for e in ranges]
        dims = tuple((i+k-1 for i in n)...)
        θ_ = zeros(eltype(values), dims...)
        nn = length(θ_)

        θ = θ_ # TODO: use MArray here
        # return θ
        # θ = MVector{dims, eltype(values), nn}(θ_...) # TODO: check whether we always want that
        ci = SplineInterpolator{typeof(ranges), typeof(θ), k}(ranges, θ)
        if !isnothing(values)
            fit!(ci, values)
        end
        
        return ci

    end

    function (spl::CubicInterpolator{G,C})(x::SVector) where G where C
        a = tuple( (e[1] for e in spl.grid)...)
        b = tuple( (e[2] for e in spl.grid)...)
        n = tuple( (e[3] for e in spl.grid)...)
        splines.eval_UC_spline(a,b,n, spl.θ, x)
    end

    function (spl::SplineInterpolator{G,C,3})(x::SVector) where G where C
        a = tuple( (e[1] for e in spl.grid)...)
        b = tuple( (e[2] for e in spl.grid)...)
        n = tuple( (e[3] for e in spl.grid)...)
        splines.eval_UC_spline(a,b,n, spl.θ, x)
    end


    function (spl::SplineInterpolator{G,C,1})(x) where G where C
        # ranges = spl.grid.ranges
        # data = spl.θ
        # dims = tuple( (e[3] for e in ranges)... )
        # v = reshape(view(data, :), dims) 
        interp(spl.grid, spl.θ, x...)
    end

end
