module splines

export filter_coeffs, interpolant_cspline, filter_coeffs
export eval_UC_spline, eval_UC_spline!
export prefilter!

include("csplines.jl")
include("splines_filter.jl")
include("cubic_prefilter.jl")


function interpolant_cspline(a, b, orders, V)

    coefs = filter_coeffs(a, b, orders, V)

    function fun(s::Array{Float64,2})
        return eval_UC_spline(a, b, orders, coefs, s)
    end

    function fun(p::Float64...)
        return fun([p...]')
    end

    return fun

end



end
