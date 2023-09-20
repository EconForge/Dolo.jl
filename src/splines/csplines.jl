include("macros.jl")

using StaticArrays


# old style call
function eval_UC_spline(a, b, orders, C::AbstractArray{Float64}, S::Matrix{Float64})
    d,N = size(S)

    n_x = size(C,1)
    dims = tuple(size(C)[2:end]...)
    out = zeros(n_x,N)
    CC = reshape(reinterpret(SVector{n_x,Float64},vec(C)),dims)
    SS = reshape(reinterpret(SVector{d,Float64},vec(S)),(N,))
    V  = reshape(reinterpret(SVector{n_x,Float64},vec(out)),(N,))
    eval_UC_spline!(a, b, orders, CC, SS, V)
    return out
end


function eval_UC_spline(a, b, orders, C::AbstractArray{T, d}, S::Vector{SVector{d,Float64}}) where T where d
    N = length(S)
    vals = zeros(T,N)
    eval_UC_spline!(a, b, orders, C, S, vals)
    return vals
end


@generated function eval_UC_spline!(a, b, orders, C::AbstractArray{T, d}, S::AbstractVector{SVector{d, Float64}}, V::AbstractVector{T}) where d where T
    # d = C.parameters[2]-1 # first dimension of C indexes the splines
    # the speed penalty of extrapolating points when iterating over point
    # seems very small so this is the default
    fun = (create_function(d,"natural"))
    return fun.args[2].args[2]
end

@generated function eval_UC_spline_(a, b, orders, C::AbstractArray{T, d}, S::SVector{d,U}) where d where T where U
    fun = create_function(d,"natural"; vectorize=false)
    return fun.args[2].args[2]
end

eval_UC_spline(a, b, orders, C::AbstractArray{T, d}, S::SVector{d,U}) where d where U where T = eval_UC_spline_(a, b, orders, C, S)
