using StaticArrays
import LoopVectorization: VectorizationBase
import Base: getindex

getindex(A::Vector{Float64}, i::VectorizationBase.Vec{4,Int64}) = VectorizationBase.Vec{4, Float64}(A[i(1)], A[i(2)], A[i(3)], A[i(4)])


# ## TODO : rewrite the following
# it is currently slower than the custom versions below

@inline function reduce_tensors(λ::SVector{d1,U}, v::SArray{T,V,d2,W}) where d1 where d2 where T where U where V where W
    vv = tuple( (SVector(1-e, e) for e in λ)...) 
    xx = SVector( (prod(e) for e in Iterators.product(vv...))...)
    return sum( xx .* v[:])
end

matextract(v::AbstractArray{T,3}, i,j,k) where T = SArray{Tuple{2, 2, 2}, T, 3, 8}(
    v[i,  j  ,k],
    v[i+1,j  ,k],
    v[i  ,j+1,k],
    v[i+1,j+1,k],
    v[i  ,j  ,k+1],
    v[i+1,j  ,k+1],
    v[i  ,j+1,k+1],
    v[i+1,j+1,k+1]
)

matextract(v::AbstractArray{T,2}, i, j) where T = SArray{Tuple{2,2}, T, 2, 4}(
    v[i,j],
    v[i+1,j],
    v[i,j+1],
    v[i+1,j+1]
)
matextract(v::AbstractArray{T,1}, i) where T = SArray{Tuple{2}, T, 1, 2}(
    v[i],
    v[i+1]
)

function interp(ranges::NTuple{d, Tuple{Float64, Float64, Int64}}, values::AbstractArray{T,d}, x::SVector{d, U}) where d where T where U
    
    a = SVector( (e[1] for e in ranges)... )
    b = SVector( (e[2] for e in ranges)... )
    n = SVector( (e[3] for e in ranges)... )

    δ = (b.-a)./(n.-1)

    i = div.( (x.-a), δ)

    i = max.(min.(i, n.-2), 0)

    λ = (x.-(a .+ δ.*i))./δ

    i_ = floor.(Int, i) .+ 1

    # inds = tuple( (SVector(e,e+1) for e in i_)...)
    # return reduce_tensors( λ, values[inds...] )

    return reduce_tensors( λ, matextract(values, i_...) )

end

function interp(ranges::NTuple{d, Tuple{Float64, Float64, Int64}}, values::AbstractArray{T,d}, x::Vararg{U}) where d where T where U
    xx = SVector(x...)
    interp(ranges, values, xx)
end


function vecinterp_1(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    for n=1:N
        out[n] = interp(ranges, values, x[n])
    end
    return out
end


function vecinterp_2(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    out .= (u->interp(ranges, values, u)).(x)
    return out
end

using LoopVectorization

function vecinterp_3(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    @simd for n=1:N
        out[n] = interp(ranges, values, x[n])
    end
    return out
end


using LoopVectorization

function vecinterp_4(ranges, values::Vector{T}, X::Vector{<:Number}) where T
    
    K = length(X)
    out = zeros(T,K)

    a = (ranges[1][1])
    b = (ranges[1][2])
    N = (ranges[1][3])

    δ = (b-a)/(N-1)

    @turbo for n=1:K

        x = X[n]

        i = div( (x-a), δ)

        i = max(min(i, N-2), 0)

        λ = (x-(a + δ*i))/δ

        i_ = floor(Int, i) + 1

        v0 = values[i_]
        v1 = values[i_+1]
        # v0 = getindex(values, i_)
        # v1 = getindex(values, i_+1)

        out[n] = (1-λ)*v0 + λ*v1
    end
    return out
end

function vecinterp_5(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    @turbo for n=1:N
        ret = interp(ranges, values, x[n])
        out[n] = ret
    end
    return out
end


