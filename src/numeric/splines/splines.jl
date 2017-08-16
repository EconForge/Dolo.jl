module splines

export filter_coeffs, interpolant_cspline, filter_coeffs
export eval_UC_spline, eval_UC_spline_G, eval_UC_multi_spline, eval_UC_spline!

include("csplines.jl")
include("splines_filter.jl")

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
