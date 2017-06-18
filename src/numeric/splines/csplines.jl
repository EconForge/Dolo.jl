include("macros.jl")

Ad = [
   [-1.0/6.0  3.0/6.0 -3.0/6.0 1.0/6.0];
   [ 3.0/6.0 -6.0/6.0  0.0/6.0 4.0/6.0];
   [-3.0/6.0  3.0/6.0  3.0/6.0 1.0/6.0];
   [ 1.0/6.0  0.0/6.0  0.0/6.0 0.0/6.0]
]

dAd = [
   [ 0.0 -0.5  1.0 -0.5];
   [ 0.0  1.5 -2.0  0.0];
   [ 0.0 -1.5  1.0  0.5];
   [ 0.0  0.5  0.0  0.0]
]

d2Ad = [
   [ 0.0 0.0 -1.0  1.0];
   [ 0.0 0.0  3.0 -2.0];
   [ 0.0 0.0 -3.0  1.0];
   [ 0.0 0.0  1.0  0.0]
]

function eval_UC_spline(smin, smax, orders, C, S)

    d = size(S,2)
    N = size(S,1)

    vals = zeros(N)
    eval_UC_spline!(smin, smax, orders, C, S, vals)
    return vals

end

function eval_UC_spline_G(a, b, orders, C, S)

    d = size(S,2)
    N = size(S,1)

    vals = zeros(N)
    grad = zeros(N,d)

    eval_UC_spline_G!(a, b, orders, C, S, vals, grad)

    return (vals, grad)

end

function eval_UC_multi_spline(a, b, orders, C, S)

    d = size(S,2)
    N = size(S,1)
    K = size(C,1) # number of splines to evaluate
    vals = zeros(K,N)
    eval_UC_multi_spline!(a, b, orders, C, S, vals)
    return vals

end

# problem with this approach: the functions don't get cached.

# fun = (create_function(1,"natural"))
# fun.args[2].args[2]
# fun.args[2].args[2]
#
# for d = 1:4
#     eval(create_function(d,"natural"))
# end
#
# for d = 1:4
#     eval(create_function_with_gradient(d,"natural"))
# end

# with this approach functions don't get cached either (yet)
# but it's cleaner

@generated function eval_UC_spline!(a, b, orders, C, S, V)
    d = C.parameters[2]
    # the speed penalty of extrapolating points when iterating over point
    # seems very small so this is the default
    fun = (create_function(d,"natural"))
    return fun.args[2].args[2]
end


@generated function eval_UC_spline_G!(a, b, orders, C, S, V, dV)
    d = C.parameters[2]
    fun = (create_function_with_gradient(d,"natural"))
    return fun.args[2].args[2]
end

@generated function eval_UC_multi_spline!(a, b, orders, C, S, V)
    d = C.parameters[2]-1 # first dimension of C indexes the splines
    # the speed penalty of extrapolating points when iterating over point
    # seems very small so this is the default
    fun = (create_function_multi_spline(d,"natural"))
    return fun.args[2].args[2]
end
