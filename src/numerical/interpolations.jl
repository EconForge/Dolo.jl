using Interpolations: BSplineInterpolation, AbstractInterpolationWrapper

"""
AbstractDoloInterpoland interface:

- constructor that takes `a, b, n, values`
- `set_values!` or `update_coefficients!` method
- `Base.call`
- `grid` method

I have two sets of dimensions to worry about:

1. The number of variables in the basis space. For dolo this should almost
always be the number of state variables in the model.
2. The number of variables in the image. This is the number of functions to
interpolate over. In dolo this is almost always the number of control variables
or value variables (for vfi) in the model.

TODO: think about this more, but here's an idea

Maybe what we should do is split the interface into 2 parts:

1. `grid{T<:AbstractDoloInterpoland}(::Type{T}, a, b, n)``
2. constructor, set_values, call/eval mutating, call/eval non-mutating

The reason for the split is that often before we can actually build the
interpoland we need to have a grid of values on which to evaluate the function
to be interpolated. For example, if we want an approximation of `f(a, b)`,
usualy we need to have a grid of `a × b`, we then feed that into `f`, and build
an interpoland using that discretized version of `f`.

"""
abstract AbstractDoloInterpoland{N}

linspace_grid(a, b, n) = [linspace(A, B, N) for (A, B, N) in zip(a, b, n)]
mlinspace(a, b, n) = gridmake(linspace_grid(a, b, n)...)::Matrix{Float64}

typealias _SplineInterp Union{AbstractInterpolationWrapper{
                                  TypeVar(:T), TypeVar(:N),
                                  TypeVar(:TS, BSplineInterpolation)
                              },
                              BSplineInterpolation}

# simple wrapper type around Interpolations.jl so we can implement our own
# methods on top of it and store other data we need
type SplineInterpoland{N,TI<:_SplineInterp} <: AbstractDoloInterpoland{N}
    a::Vec{N,Float64}  # TODO: parse a,b,n in as Vec instead of Vector
    b::Vec{N,Float64}
    n::Vec{N,Int}
    itp::TI
end

# interface part 1
grid(::Type{SplineInterpoland}, a, b, n) = linspace_grid(a, b, n)

# this is the main interface we want.
function CubicSplines{N}(a::Vec{N,Float64}, b::Vec{N,Float64},
                         n::Vec{N,Int}, values, bc::Interpolations.Flag)
    itp = interpolate(values, BSpline(Cubic(bc)), OnGrid())
    lsgrid = linspace_grid(a, b, n)
    sitp = scale(itp, lsgrid...)  # TODO: type unstable right here
    SplineInterpoland(a, b, n, sitp)
end

CubicSplines{N,T<:Interpolations.Flag}(a::Vec{N,Float64}, b::Vec{N,Float64},
                                       n::Vec{N,Int}, values,
                                       bc::Type{T}=Natural) =
    CubicSplines(a, b, n, values, T())
