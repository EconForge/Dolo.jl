# bring in some extra unexported names here to make code below readable
using Interpolations: BSplineInterpolation, AbstractInterpolationWrapper, Flag

"""
The `AbstractDoloInterpoland` is made of three parts:

1. Required fields:
    - `a::FixedSizeArrays.Vec{N,Float64}`: a `Vec` of lower borders for the
    domain in each of the `N` dimensions
    - `b::FixedSizeArrays.Vec{N,Float64}`: a `Vec` of upper borders for the
    domain in each of the `N` dimensions
    - `n::FixedSizeArrays.Vec{N,Float64}`: a `Vec` for the order of
    approximation or number of points (depending on interpolation type) for
    each of the `N` dimensions
2. Methods for before the values on the grid are known:
    - `grid{T<:AbstractDoloInterpoland}(::Type{T}, a, b, n)`: produces the grid
    for the type of interpolation given lower bounds a, upper bounds b, and
    orders n.
3. Methods for after the values on the grid are known
    - A constructor that takes `a, b, n, values, other_args...`: where `values`
    is the value of the function on the grid and `other_args` are any other
    arguments needed for the interpolation type to be fully specified. Each
    argument in `other_args` should have a default value so that some form of
    all interpolation types can be construced with just `a, b, n, values`
    - Methods for `set_values!(dr::AbstractDoloInterpoland,  new_values)`.
    This function will update the values on the grid. If possible, the routine
    mutate `dr`, but regardless will _always_ return an updated version of
    the interpolation type. Users should do `dr = set_values!(dr, new_values)`
    just in case it is not possible to directly mutate the origial `dr`.
    - A method `evalaute(dr, xs::Real...)` that evaluate `dr` at a single point

The methods below are provided by the interface, conditional on the above
methods being provided. All methods below can be implemented for specific
subtypes if there are obvious efficiency gains to be had over the default
implementation. For methods where an "equivalent to" comparison is provided,
the actual implementation may differ from the stated equivalence.

- `grid{T<:AbstractDoloInterpoland}(::T)`
- `evaluate(dr, xs::AbstractVector...)`: evaluates `dr` at the `i`th  point in
each of the vectors. This is equivalent to `map(x->evaluate(dr, x...), xs...)`
- `evaluate!(dr, out::AbstractVector, xs::AbstractVector...)`: evaluates `dr`
at the `i`th  point in each of the vectors and stores the result in `out`.
This is equivalent to `map!(x->evaluate(dr, x...), out, xs...)`
- `evaluate(dr, xs::AbstractMatrix)`: evaluates `dr` at `size(xs, 1)` points
by assuming the `i`th column of `xs` constitues the points for the `i`th
dimension of the basis of `dr`. Equivalent to
`mapslices(x->evaluate(dr, x...), 1)`
- `evaluate(dr, out::AbstractVector, xs::AbstractMatrix)`: evaluates `dr` at
`size(xs, 1)` points by assuming the `i`th column of `xs` constitues the points
for the `i`th dimension of the basis of `dr`. Equivalent to
`out[:] = mapslices(x->evaluate(dr, x...), 1)`
- `update_coefs!(dr, new_vals)`: another name for `set_values!(dr, new_vals)`

"""
abstract AbstractDoloInterpoland{N}

typealias _SplineInterp Union{AbstractInterpolationWrapper{
                                  TypeVar(:T), TypeVar(:N),
                                  TypeVar(:TS, BSplineInterpolation)
                              },
                              BSplineInterpolation}

# simple wrapper type around Interpolations.jl so we can implement our own
# methods on top of it and store other data we need
type SplineInterpoland{N,TI<:_SplineInterp} <: AbstractDoloInterpoland{N}
    a::Vec{N,Float64}
    b::Vec{N,Float64}
    n::Vec{N,Int}
    itp::TI
end

## interface part 1.
# NOTE type param allows this to work for all type params on SplineInterpoland
grid{T<:SplineInterpoland}(::Type{T}, a, b, n) = linspace_grid(a, b, n)
mgrid{T<:SplineInterpoland}(::Type{T}, a, b, n) = mlinspace(a, b, n)

## interface part 2
# constructors
function CubicSplines{N}(a::Vec{N,Float64}, b::Vec{N,Float64}, n::Vec{N,Int},
                         values, bc::Flag)
    itp = interpolate(values, BSpline(Cubic(bc)), OnGrid())
    lsgrid = linspace_grid(a, b, n)
    sitp = scale(itp, lsgrid...)  # TODO: type unstable right here
    SplineInterpoland(a, b, n, sitp)
end

CubicSplines{N,T<:Flag}(a::Vec{N,Float64}, b::Vec{N,Float64}, n::Vec{N,Int},
                        values, bc::Type{T}=Natural) =
    CubicSplines(a, b, n, values, T())

function CubicSplines{N,T<:Flag}(approx::Approximation{N}, values,
                                 bc::Union{T, Type{T}}=Natural)
    CubicSplines(approx.a, approx.b, approx.n, values, bc)
end

# set_values! TODO
# evaluate(::SplineInterpoland, xs::Real...) TODO


# provided API methods on `AbstractDoloInterpoland`
# TODO: implement the other ones documented above
grid{T<:AbstractDoloInterpoland}(dr::T) = grid(T, dr.a, dr.b, dr.n)
mgrid{T<:AbstractDoloInterpoland}(dr::T) = mgrid(T, dr.a, dr.b, dr.n)
