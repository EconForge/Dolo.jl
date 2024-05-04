using StaticArrays
import LoopVectorization: VectorizationBase
import Base: getindex

getindex(A::Vector{Float64}, i::VectorizationBase.Vec{4,Int64}) = VectorizationBase.Vec{4, Float64}(A[i(1)], A[i(2)], A[i(3)], A[i(4)])


# ## TODO : rewrite the following
## this version allocates for d>=2
@inline function reduce_tensors(λ::SVector{d1,U}, v::SArray{T,V,d2,W}) where d1 where d2 where T where U where V where W
    vv = tuple( (SVector(1-e, e) for e in λ)...) 
    xx = SVector( (prod(e) for e in Iterators.product(vv...))...)
    return sum( xx .* v[:])
end

@inline function reduce_tensors(λ::SVector{2,U}, M::SArray{T,V,2,W}) where d1 where d2 where T where U where V where W
    v1 = SVector(1-λ[1],λ[1])
    v2 = SVector(1-λ[2],λ[2])
    return M[1,1]*v1[1]*v2[1] + M[1,2]*v1[1]*v2[2] + M[2,1]*v1[2]*v2[1] + M[2,2]*v1[2]*v2[2]
    # vv = tuple( (SVector(1-e, e) for e in λ)...) 
    # xx = SVector( (prod(e) for e in Iterators.product(vv...))...)
    # return sum( xx .* v[:])
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

    M = matextract(values, i_...)
    
    return reduce_tensors(λ, M)
    # return reduce_tensors( λ, matextract(values, i_...) )

end

function interp(ranges::NTuple{d, Tuple{Float64, Float64, Int64}}, values::AbstractArray{T,d}, x::Vararg{U}) where d where T where U
    xx = SVector(x...)
    interp(ranges, values, xx)
end
