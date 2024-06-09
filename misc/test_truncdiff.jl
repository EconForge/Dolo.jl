using Dolo: splines


ranges = ((0.0, 1.0, 10), (0.0, 1.0, 20))
methods(splines.interp)

using StaticArrays

x0 = SVector( 10,20)


values = rand(typeof(x0), 10,20)

values2 = rand(typeof(x0), 12,22)

using ForwardDiff


fun = u->splines.interp(ranges, values, u)

ForwardDiff.jacobian(
    fun,
    x0
)

a = tuple((e[1] for e in ranges)...)
b = tuple((e[2] for e in ranges)...)
orders = tuple((e[3] for e in ranges)...)

fun2 = u->splines.eval_UC_spline(a,b,orders, values2, u)
fun2(x0)




ForwardDiff.jacobian(
    fun2,
    x0
)

# import ForwardDiff: unsafe_trunc

# function ForwardDiff.unsafe_trunc(TI, number::ForwardDiff.Dual{Tag, Tf}) where Tag where Tf<:Union{Float32, Float64}
#     ForwardDiff.Dual{Tag}(
#         unsafe_trunc(TI, ForwardDiff.value(number)),
#         0*ForwardDiff.partials(number)
#     )

# end