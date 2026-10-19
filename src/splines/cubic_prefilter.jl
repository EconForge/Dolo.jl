### First Version of Prefilter

@inline function solve_coefficients!(
            bands::AbstractMatrix{Float64},
            bb::AbstractVector{T},
            coefs::AbstractVector{T}
    ) where T

    M = size(coefs,1) - 2

    # Solve interpolating equations
    # First and last rows are different

    bands[1, 2] /= bands[1, 1]
    bands[1, 3] /= bands[1, 1]
    bb[1] /= bands[1, 1]
    bands[1, 1] = 1.0
    bands[2, 2] -= bands[2, 1] * bands[1, 2]
    bands[2, 3] -= bands[2, 1] * bands[1, 3]
    bb[2] -= bands[2, 1] * bb[1]
    bands[1, 1] = 0.0
    bands[2, 3] /= bands[2, 2]
    bb[2] /= bands[2, 2]
    bands[2, 2] = 1.0

    # Now do rows 2 through M+1
    for row = 3:(M+1)
        bands[row, 2] -= bands[row, 1] * bands[row - 1, 3]
        bb[row] -= bands[row, 1] * bb[row - 1]
        bands[row, 3] /= bands[row, 2]
        bb[row] /= bands[row, 2]
        bands[row, 1] = 0.0
        bands[row, 2] = 1.0
    end
    # Do last row
    bands[M + 2, 2] -= bands[M + 2, 1] * bands[M, 3]
    bb[M + 2] -= bands[M + 2, 1] * bb[M]
    bands[M + 2, 3] -= bands[M + 2, 2] * bands[M + 1, 3]
    bb[M + 2] -= bands[M + 2, 2] * bb[M + 1]
    bb[M + 2] /= bands[M + 2, 3]
    bands[M + 2, 3] = 1.0

    coefs[M + 2] = bb[M + 2]
    # Now back substitute up
    # for row in (M+1):2
    for ii in 1:M
        row = M+2 - ii
        coefs[row] = bb[row] - bands[row, 3] * coefs[row + 1]
    end
    # Finish with first row
    coefs[1] = bb[1] - bands[1, 2] * coefs[2] - bands[1, 3] * coefs[3]

    return nothing
end

@inline function fill_bands!(M::Int, bands::AbstractMatrix, bb::AbstractVector{T}, data::AbstractVector{T}) where T

    # boundary conditions
    @inbounds bands[1, 1] = 1.0
    @inbounds bands[1, 2] = -2.0
    @inbounds bands[1, 3] = 1.0
    @inbounds bands[M+2, 1] = 1.0
    @inbounds bands[M+2, 2] = -2.0
    @inbounds bands[M+2, 3] = 1.0

    for k=2:M+1
        @inbounds bands[k,1] = 1.0/6.0
        @inbounds bands[k,2] = 2.0/3.0
        @inbounds bands[k,3] = 1.0/6.0
    end

    @inbounds bb .= data
    @inbounds bb[1] = zero(T)
    @inbounds bb[M+2] = zero(T)

    return nothing

end


# TODO: there are certainly some cases where we don't want MVectors there

function prefilter_1!(data::AbstractVector{T}) where T

    N = length(data)

    bands = zero(MMatrix{N, 3, Float64, N*3})
    bb = zero(MVector{N, T})

    fill_bands!(N-2, bands, bb, data)
    solve_coefficients!(bands, bb, data)

end

###
### Second Version of Prefilter
###
### in-place without memory allocation but


const zq_0 = sqrt(3.0)-2.0
const λ_0 = 6.0

getprecision(g::AbstractArray{T}) where T<:Number = T
getprecision(g::AbstractArray{T}) where T<:SVector = eltype(g).types[1].types[1]


function prefilter_2!(c::AbstractVector{T}) where T
# function prefilter!(c::AbstractVector{T}) where T


    prec = getprecision(c)

    λ = convert(prec, λ_0)
    zq = convert(prec, zq_0)
    
    N = length(c)
    Horizon = (min(20, N-2))

    f1 = c[2]
    ff = c[N-1]

    t1 = zero(eltype(c)) #0.0
    zn = convert(prec,1.0)
    for n in 0:Horizon
        t1 += λ*(2*f1-c[2+n])*zn
        zn *= zq 
    end

    c[2] = t1
    c[1] = (c[2]-λ*f1)/zq

    @inbounds for k in 3:N
        c[k] = λ*c[k] + zq*c[k-1]
    end

    c[N-1] = ff
    @inbounds for kk in 1:(N-2)
        k = N-1-kk
        c[k] = zq*(c[k+1] - c[k])
    end
    c[N] = 2*ff-c[N-2]

    return nothing

end


###
### Multidimensional arrays
###
### θ holds the values to interpolate, with one point of padding on each side in each
### dimension. The 1D filter is applied along each dimension, on every line of θ.

using KernelAbstractions: @kernel, @index, get_backend, synchronize, CPU

"""
    prefilter!(θ)

Replace the values in `θ` by the coefficients of the cubic spline interpolating them.
The sequential algorithm (`prefilter_2!`) is used on CPU, the chunked one elsewhere (e.g. GPU).
"""
prefilter!(θ::AbstractArray) =
    get_backend(θ) isa CPU ? prefilter_sequential!(θ) : prefilter_chunked!(θ, similar(θ))

# the line of A along dimension `dim`, at index I of the other dimensions
@inline line(A, ::Val{dim}, I::CartesianIndex) where dim =
    view(A, ntuple(i -> i < dim ? I[i] : i == dim ? Colon() : I[i-1], Val(ndims(A)))...)

# size of A without dimension `dim`
other_size(A, dim) = (size(A)[1:dim-1]..., size(A)[dim+1:end]...)

# sequential algorithm along each line; independent lines are processed by different threads
function prefilter_sequential!(θ::AbstractArray)
    for dim in 1:ndims(θ)
        lines = CartesianIndices(other_size(θ, dim))
        if length(lines) == 1
            prefilter_2!(line(θ, Val(dim), first(lines)))
        else
            Threads.@threads for I in lines
                prefilter_2!(line(θ, Val(dim), I))
            end
        end
    end
    return θ
end


## Chunked version
#
# Each sweep of prefilter_2! is a linear recurrence with pole zq ≈ -0.268: the influence of
# its initial condition decays like |zq|^k, below machine precision after WARMUP=32 steps.
# A line can therefore be cut in chunks of CHUNK points, each filtered independently (by its own
# thread) after a warm-up on the WARMUP preceding points. Same result as prefilter_2!, to machine
# precision, but with one thread per chunk instead of one per line.
#
# It works out of place: the causal sweep writes in a buffer, the anticausal one back in θ.

const CHUNK = 128
const WARMUP = 32

# Causal sweep on chunk b of a line: y[k] = λ x[k] + zq y[k-1] for k ≥ 3, from the
# initial condition of prefilter_2!.
@inline function causal_chunk!(y, x, b)
    N = length(x)
    Tp = eltype(eltype(x))
    λ = convert(Tp, λ_0)
    zq = convert(Tp, zq_0)
    s = 3 + (b-1)*CHUNK                 # chunk: s:e
    e = min(N, s + CHUNK - 1)
    p = max(3, s - WARMUP)              # start of the warm-up
    acc = zero(eltype(x))
    if p == 3                           # the warm-up reaches the start of the line: exact initial value
        f1 = x[2]
        zn = one(Tp)
        for n in 0:min(20, N-2)
            acc += λ*(2*f1 - x[2+n])*zn
            zn *= zq
        end
        if b == 1
            y[2] = acc
            y[1] = (acc - λ*f1)/zq
        end
    end
    for k in p:e
        acc = λ*x[k] + zq*acc
        k >= s && (y[k] = acc)
    end
    e >= N-1 && (y[N-1] = x[N-1])       # last value, used by the anticausal sweep
    return nothing
end

# Anticausal sweep on chunk b of a line: c[k] = zq (c[k+1] - y[k]) for k = N-2 … 1,
# from c[N-1] = y[N-1]; then c[N] = 2c[N-1] - c[N-2].
@inline function anticausal_chunk!(c, y, b)
    N = length(y)
    zq = convert(eltype(eltype(y)), zq_0)
    e = N - 2 - (b-1)*CHUNK             # chunk: e down to s
    s = max(1, e - CHUNK + 1)
    p = min(N - 2, e + WARMUP)          # start of the warm-up
    ff = y[N-1]
    acc = p == N - 2 ? ff : zero(eltype(y))
    for k in p:-1:s
        acc = zq*(acc - y[k])
        k <= e && (c[k] = acc)
    end
    if b == 1
        c[N-1] = ff
        c[N] = 2*ff - c[N-2]
    end
    return nothing
end

# one thread per (chunk, line)
@kernel function causal_kernel!(Y, X, dim)
    c = @index(Global, Cartesian)
    I = CartesianIndex(Base.tail(Tuple(c)))
    causal_chunk!(line(Y, dim, I), line(X, dim, I), c[1])
end

@kernel function anticausal_kernel!(C, Y, dim)
    c = @index(Global, Cartesian)
    I = CartesianIndex(Base.tail(Tuple(c)))
    anticausal_chunk!(line(C, dim, I), line(Y, dim, I), c[1])
end

"Chunked prefilter, using `buffer` (same size as θ)."
function prefilter_chunked!(θ::AbstractArray, buffer::AbstractArray)
    backend = get_backend(θ)
    for dim in 1:ndims(θ)
        ndrange = (cld(size(θ, dim) - 2, CHUNK), other_size(θ, dim)...)
        causal_kernel!(backend)(buffer, θ, Val(dim); ndrange)
        anticausal_kernel!(backend)(θ, buffer, Val(dim); ndrange)
    end
    synchronize(backend)
    return θ
end
