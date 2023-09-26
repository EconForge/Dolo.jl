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

# @inline function reduce_tensors(λ::SVector{3,U}, v::SArray{T,Float64,3,V}) where d1 where d2 where T where U where V
#     (1.0 - λ[3]) * (             
#         (1.0 - λ[2])*(
#             (1.0 - λ[1])*v[1,1,1] + λ[1]*v[2,1,1]
#         ) + λ[2] * (
#             (1.0 - λ[1])*v[1,2,1] + λ[1]*v[2,2,1]
#         )
#     )
#     + λ[3] * (
#         (1.0 - λ[2])*(
#             (1.0 - λ[1])*v[1,1,2] + λ[1]*v[2,1,2]
#         ) + λ[2] * (
#             (1.0 - λ[1])*v[1,2,2] + λ[1]*v[2,2,2]
#         )
#     )
# end




# @inline function interp(ranges::Tuple{Tuple{Float64, Float64, Int64}}, values::AbstractArray{T}, x::U) where T where U
# # function interp(ranges::Tuple{Tuple{Float64, Float64, Int64}}, values, x) where T where U 

#     a = (ranges[1][1])
#     b = (ranges[1][2])
#     n = (ranges[1][3])

#     δ = (b-a)/(n-1)

#     i = div( (x-a), δ)

#     i = max(min(i, n-2), 0)

#     λ = (x-(a + δ*i))/δ

#     i_ = floor(Int, i) + 1

#     v0 = values[i_]
#     v1 = values[i_+1]

#     (1.0 - λ)*v0 + λ*v1

# end

# @inline function interp(ranges::Tuple{Tuple{Float64, Float64, Int64},Tuple{Float64, Float64, Int64}}, values::AbstractArray{T}, x_1::U, x_2::U) where T where U
#     # function interp(ranges::Tuple{Tuple{Float64, Float64, Int64}}, values, x) where T where U 
    
#         a_1 = (ranges[1][1])
#         b_1 = (ranges[1][2])
#         n_1 = (ranges[1][3])
    
#         a_2 = (ranges[2][1])
#         b_2 = (ranges[2][2])
#         n_2 = (ranges[2][3])

#         δ_1 = (b_1-a_1)/(n_1-1)
#         δ_2 = (b_2-a_2)/(n_2-1)
    
#         i_1 = div( (x_1-a_1), δ_1)
#         i_2 = div( (x_2-a_2), δ_2)
    
#         i_1 = max(min(i_1, n_1-2), 0)
#         i_2 = max(min(i_2, n_2-2), 0)
    
#         λ_1 = (x_1-(a_1 + δ_1*i_1))/δ_1
#         λ_2 = (x_2-(a_2 + δ_2*i_2))/δ_2
    
#         i_1_ = floor(Int, i_1) + 1
#         i_2_ = floor(Int, i_2) + 1

#         v00 = values[i_1_, i_2_]
#         v01 = values[i_1_, i_2_+1]
#         v10 = values[i_1_+1, i_2_]
#         v11 = values[i_1_+1, i_2_+1]
    
#         res = (1.0 - λ_2)*(
#             (1.0 - λ_1)*v00 + λ_1*v10
#         ) + λ_2 * (
#             (1.0 - λ_1)*v01 + λ_1*v11
#         )
    
#         return res

# end


# @inline function interp(ranges::Tuple{
#             Tuple{Float64, Float64, Int64},
#             Tuple{Float64, Float64, Int64},
#             Tuple{Float64, Float64, Int64}
#         }, values::AbstractArray{T}, x_1::U, x_2::U, x_3::U) where T where U
    
#         a_1 = (ranges[1][1])
#         b_1 = (ranges[1][2])
#         n_1 = (ranges[1][3])
    
#         a_2 = (ranges[2][1])
#         b_2 = (ranges[2][2])
#         n_2 = (ranges[2][3])

#         a_3 = (ranges[3][1])
#         b_3 = (ranges[3][2])
#         n_3 = (ranges[3][3])

#         δ_1 = (b_1-a_1)/(n_1-1)
#         δ_2 = (b_2-a_2)/(n_2-1)
#         δ_3 = (b_3-a_3)/(n_3-1)
    
#         i_1 = div( (x_1-a_1), δ_1)
#         i_2 = div( (x_2-a_2), δ_2)
#         i_3 = div( (x_3-a_3), δ_3)


#         i_1 = max(min(i_1, n_1-2), 0)
#         i_2 = max(min(i_2, n_2-2), 0)
#         i_3 = max(min(i_3, n_3-2), 0)
    
#         λ_1 = (x_1-(a_1 + δ_1*i_1))/δ_1
#         λ_2 = (x_2-(a_2 + δ_2*i_2))/δ_2
#         λ_3 = (x_3-(a_3 + δ_3*i_3))/δ_3
    
#         i_1_ = floor(Int, i_1) + 1
#         i_2_ = floor(Int, i_2) + 1
#         i_3_ = floor(Int, i_3) + 1

#         v000 = values[i_1_,   i_2_,   i_3_]
#         v001 = values[i_1_,   i_2_,   i_3_+1]
#         v010 = values[i_1_,   i_2_+1, i_3_]
#         v011 = values[i_1_,   i_2_+1, i_3_+1]
#         v100 = values[i_1_+1, i_2_,   i_3_]
#         v101 = values[i_1_+1, i_2_,   i_3_+1]
#         v110 = values[i_1_+1, i_2_+1, i_3_]
#         v111 = values[i_1_+1, i_2_+1, i_3_+1]


#         res = (1.0 - λ_3) * (             
#                 (1.0 - λ_2)*(
#                     (1.0 - λ_1)*v000 + λ_1*v100
#                 ) + λ_2 * (
#                     (1.0 - λ_1)*v010 + λ_1*v110
#                 )
#             )
#             + λ_3 * (
#                 (1.0 - λ_2)*(
#                     (1.0 - λ_1)*v001 + λ_1*v101
#                 ) + λ_2 * (
#                     (1.0 - λ_1)*v011 + λ_1*v111
#                 )
#             )

#         return res

# end


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





# @code_warntype interp(ranges, values, 2.0)

# xvec = [range(-0.1, 2.1; length=1000)...]
# tvec = xvec .^ 2

# yvec = [(interp(ranges, values, x)) for x in xvec]


# using ForwardDiff

# ForwardDiff.jacobian(u->[interp(ranges, values, u[1])], [0.3])

# using Plots

# plot(xvec, tvec)
# plot!(xvec, yvec)

# using LoopVectorization: VectorizationBase


# a,b,N = ranges[1]
# x = VectorizationBase.Vec{4, Float64}(4, 3, 2, 5)
# δ = (b-a)/(N-1)

# i = div( (x-a), δ)

# i = max(min(i, N-2), 0)

# λ = (x-(a + δ*i))/δ

# i_ = floor(Int, i) + 1

# v0 = values[i_]
# v1 = values[i_+1]

# out[n] = (1-λ)*v0 + λ*v1



# import Base: getindex
# getindex(A::Vector{T}, i::VectorizationBase.Vec{4,Int64}) where T where d = VectorizationBase.Vec{4, T}(A[i(1)], A[i(2)], A[i(3)], A[i(4)])

# getindex(values, i_)


# values[i_]

# i_(1)