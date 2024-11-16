include("macros.jl")

using StaticArrays


# old style call
function eval_UC_spline(a, b, orders, C::AbstractArray{Tf}, S::Matrix{Tf}) where Tf

    d,N = size(S)

    n_x = size(C,1)
    dims = tuple(size(C)[2:end]...)
    out = zeros(n_x,N)
    CC = reshape(reinterpret(SVector{n_x,Tf},vec(C)),dims)
    SS = reshape(reinterpret(SVector{d,Tf},vec(S)),(N,))
    V  = reshape(reinterpret(SVector{n_x,Tf},vec(out)),(N,))
    eval_UC_spline!(a, b, orders, CC, SS, V)
    return out
end


function eval_UC_spline(a, b, orders, C::AbstractArray{T, d}, S::Vector{SVector{d,Tf}}) where T where d where Tf
    N = length(S)
    vals = zeros(T,N)
    eval_UC_spline!(a, b, orders, C, S, vals)
    return vals
end


@generated function eval_UC_spline!(a, b, orders, C::AbstractArray{T, d}, S::AbstractVector{SVector{d}}, V::AbstractVector{T}) where d where T
    # d = C.parameters[2]-1 # first dimension of C indexes the splines
    # the speed penalty of extrapolating points when iterating over point
    # seems very small so this is the default
    Tf = eltype(a)

    fun = (create_function(d,"natural"; Tf=Tf))
    return fun.args[2].args[2]
end


@generated function eval_UC_spline_(a, b, orders, C::AbstractArray{T, d}, S::SVector{d,U}) where d where T where U
    Tf = eltype(a)
    fun = create_function(d,"natural"; Tf=Tf, vectorize=false)
    return fun.args[2].args[2]
end

function eval_UC_spline(a, b, orders, C::AbstractArray{T, d}, S::SVector{d,U}) where d where U where T 
    eval_UC_spline_(a, b, orders, C, S)
end


function cPhi(λ)

    T = typeof(λ)

    λ0 = convert(T,1)
    λ1 = λ
    λ2 = λ*λ1
    λ3 = λ*λ2

    if λ1<0
        Phi_4 = convert(T,-0.5) * λ1 + convert(T, 0.16666666666666666)
        Phi_3 = convert(T, 0.0) * λ1 + convert(T, 0.6666666666666666)
        Phi_2 = convert(T, 0.5) * λ1 + convert(T, 0.16666666666666666)
        Phi_1 = convert(T, 0.0) * λ1 + convert(T, 0.0)
    elseif λ1>1
        # Phi_4 = (3 * -0.16666666666666666 + 2 * 0.5 + -0.5) * (1 - 1) + (-0.16666666666666666 + 0.5 + -0.5 + 0.16666666666666666)
        # Phi_3 = (3 * 0.5 + 2 * -1.0 + 0.0) * (λ1 - 1) + (0.5 + -1.0 + 0.0 + 0.6666666666666666)
        # Phi_2 = (3 * -0.5 + 2 * 0.5 + 0.5) * (λ1 - 1) + (-0.5 + 0.5 + 0.5 + 0.16666666666666666)
        # Phi_1 = (3 * 0.16666666666666666 + 2 * 0.0 + 0.0) * (λ1 - 1) + (0.16666666666666666 + 0.0 + 0.0 + 0.0)
        Phi_4 = convert(T,0.0)
        Phi_3 = (-convert(T,0.5)) * (λ1 - convert(T,1)) + (-convert(T,0.5) + convert(T,0.6666666666666666))
        Phi_2 = convert(T, (0.5 + 0.16666666666666666))
        Phi_1 = 3 * convert(T,0.16666666666666666)  * (λ1 - 1) + convert(T,0.16666666666666666)

    else
        Phi_1 = convert(T,-0.16666666666666666) * λ3 + convert(T,0.5) * λ2 + convert(T,-0.5) * λ1 + convert(T,0.16666666666666666) * λ0
        Phi_2 = convert(T, 0.5) * λ3 + convert(T,-1.0) * λ2 + convert(T,0.6666666666666666) * λ0
        Phi_3 = convert(T,-0.5) * λ3 + convert(T,0.5) * λ2 + convert(T,0.5) * λ1 + convert(T,0.16666666666666666) * λ0
        Phi_4 = convert(T,0.16666666666666666) * λ3
    end

    
    return SVector(Phi_4, Phi_3, Phi_2, Phi_1)
end


function mextract(c::AbstractArray{T,d}, inds::NTuple{d,<:Int}, ::Val{n}) where T where d where n
    ranges = tuple((i:i+n-1 for i in inds)...)
    return @inbounds SArray{NTuple{d,n},T,d,n^d}(view(c, ranges...))
end

@inline function reduce_tensors( v::SArray, Φ::NTuple{d,SVector{n,U}}) where d where n where U
    xx = SVector( (prod(e) for e in Iterators.product(Φ...))...)
    return sum( xx .* v[:])
end

function reduce_tensors(v::SArray, Φ::Tuple{SVector{4,U}}) where U
    @inbounds v[1]*Φ[1][1] + v[2]*Φ[1][2] + v[3]*Φ[1][3] + v[4]*Φ[1][4]
end


function reduce_tensors(v::SArray, Φ::Tuple{SVector{4,U},SVector{4,U}}) where U
    @inbounds v[1,1]*Φ[1][1]*Φ[2][1] +
              v[1,2]*Φ[1][1]*Φ[2][2] +
              v[1,3]*Φ[1][1]*Φ[2][3] +
              v[1,4]*Φ[1][1]*Φ[2][4] +
              v[2,1]*Φ[1][2]*Φ[2][1] +
              v[2,2]*Φ[1][2]*Φ[2][2] +
              v[2,3]*Φ[1][2]*Φ[2][3] +
              v[2,4]*Φ[1][2]*Φ[2][4] +
              v[3,1]*Φ[1][3]*Φ[2][1] +
              v[3,2]*Φ[1][3]*Φ[2][2] +
              v[3,3]*Φ[1][3]*Φ[2][3] +
              v[3,4]*Φ[1][3]*Φ[2][4] +
              v[4,1]*Φ[1][4]*Φ[2][1] +
              v[4,2]*Φ[1][4]*Φ[2][2] +
              v[4,3]*Φ[1][4]*Φ[2][3] +
              v[4,4]*Φ[1][4]*Φ[2][4]    
end

function reduce_tensors(v::SArray, Φ::Tuple{SVector{4,U},SVector{4,U},SVector{4,U}}) where U
    a = @inbounds  v[1,1,1]*Φ[1][1]*Φ[2][1]*Φ[3][1] +
                   v[1,1,2]*Φ[1][1]*Φ[2][1]*Φ[3][2] +
                   v[1,1,3]*Φ[1][1]*Φ[2][1]*Φ[3][3] +
                   v[1,1,4]*Φ[1][1]*Φ[2][1]*Φ[3][4] +
                   v[1,2,1]*Φ[1][1]*Φ[2][2]*Φ[3][1] +
                   v[1,2,2]*Φ[1][1]*Φ[2][2]*Φ[3][2] +
                   v[1,2,3]*Φ[1][1]*Φ[2][2]*Φ[3][3] +
                   v[1,2,4]*Φ[1][1]*Φ[2][2]*Φ[3][4] +
                   v[1,3,1]*Φ[1][1]*Φ[2][3]*Φ[3][1] +
                   v[1,3,2]*Φ[1][1]*Φ[2][3]*Φ[3][2] +
                   v[1,3,3]*Φ[1][1]*Φ[2][3]*Φ[3][3] +
                   v[1,3,4]*Φ[1][1]*Φ[2][3]*Φ[3][4] +
                   v[1,4,1]*Φ[1][1]*Φ[2][4]*Φ[3][1] +
                   v[1,4,2]*Φ[1][1]*Φ[2][4]*Φ[3][2] +
                   v[1,4,3]*Φ[1][1]*Φ[2][4]*Φ[3][3] +
                   v[1,4,4]*Φ[1][1]*Φ[2][4]*Φ[3][4]
    b =  @inbounds v[2,1,1]*Φ[1][2]*Φ[2][1]*Φ[3][1] +
                   v[2,1,2]*Φ[1][2]*Φ[2][1]*Φ[3][2] +
                   v[2,1,3]*Φ[1][2]*Φ[2][1]*Φ[3][3] +
                   v[2,1,4]*Φ[1][2]*Φ[2][1]*Φ[3][4] +
                   v[2,2,1]*Φ[1][2]*Φ[2][2]*Φ[3][1] +
                   v[2,2,2]*Φ[1][2]*Φ[2][2]*Φ[3][2] +
                   v[2,2,3]*Φ[1][2]*Φ[2][2]*Φ[3][3] +
                   v[2,2,4]*Φ[1][2]*Φ[2][2]*Φ[3][4] +
                   v[2,3,1]*Φ[1][2]*Φ[2][3]*Φ[3][1] +
                   v[2,3,2]*Φ[1][2]*Φ[2][3]*Φ[3][2] +
                   v[2,3,3]*Φ[1][2]*Φ[2][3]*Φ[3][3] +
                   v[2,3,4]*Φ[1][2]*Φ[2][3]*Φ[3][4] +
                   v[2,4,1]*Φ[1][2]*Φ[2][4]*Φ[3][1] +
                   v[2,4,2]*Φ[1][2]*Φ[2][4]*Φ[3][2] +
                   v[2,4,3]*Φ[1][2]*Φ[2][4]*Φ[3][3] +
                   v[2,4,4]*Φ[1][2]*Φ[2][4]*Φ[3][4]
    c = @inbounds  v[3,1,1]*Φ[1][3]*Φ[2][1]*Φ[3][1] +
                   v[3,1,2]*Φ[1][3]*Φ[2][1]*Φ[3][2] +
                   v[3,1,3]*Φ[1][3]*Φ[2][1]*Φ[3][3] +
                   v[3,1,4]*Φ[1][3]*Φ[2][1]*Φ[3][4] +
                   v[3,2,1]*Φ[1][3]*Φ[2][2]*Φ[3][1] +
                   v[3,2,2]*Φ[1][3]*Φ[2][2]*Φ[3][2] +
                   v[3,2,3]*Φ[1][3]*Φ[2][2]*Φ[3][3] +
                   v[3,2,4]*Φ[1][3]*Φ[2][2]*Φ[3][4] +
                   v[3,3,1]*Φ[1][3]*Φ[2][3]*Φ[3][1] +
                   v[3,3,2]*Φ[1][3]*Φ[2][3]*Φ[3][2] +
                   v[3,3,3]*Φ[1][3]*Φ[2][3]*Φ[3][3] +
                   v[3,3,4]*Φ[1][3]*Φ[2][3]*Φ[3][4] +
                   v[3,4,1]*Φ[1][3]*Φ[2][4]*Φ[3][1] +
                   v[3,4,2]*Φ[1][3]*Φ[2][4]*Φ[3][2] +
                   v[3,4,3]*Φ[1][3]*Φ[2][4]*Φ[3][3] +
                   v[3,4,4]*Φ[1][3]*Φ[2][4]*Φ[3][4]
    d = @inbounds  v[4,1,1]*Φ[1][4]*Φ[2][1]*Φ[3][1] +
                   v[4,1,2]*Φ[1][4]*Φ[2][1]*Φ[3][2] +
                   v[4,1,3]*Φ[1][4]*Φ[2][1]*Φ[3][3] +
                   v[4,1,4]*Φ[1][4]*Φ[2][1]*Φ[3][4] +
                   v[4,2,1]*Φ[1][4]*Φ[2][2]*Φ[3][1] +
                   v[4,2,2]*Φ[1][4]*Φ[2][2]*Φ[3][2] +
                   v[4,2,3]*Φ[1][4]*Φ[2][2]*Φ[3][3] +
                   v[4,2,4]*Φ[1][4]*Φ[2][2]*Φ[3][4] +
                   v[4,3,1]*Φ[1][4]*Φ[2][3]*Φ[3][1] +
                   v[4,3,2]*Φ[1][4]*Φ[2][3]*Φ[3][2] +
                   v[4,3,3]*Φ[1][4]*Φ[2][3]*Φ[3][3] +
                   v[4,3,4]*Φ[1][4]*Φ[2][3]*Φ[3][4] +
                   v[4,4,1]*Φ[1][4]*Φ[2][4]*Φ[3][1] +
                   v[4,4,2]*Φ[1][4]*Φ[2][4]*Φ[3][2] +
                   v[4,4,3]*Φ[1][4]*Φ[2][4]*Φ[3][3] +
                   v[4,4,4]*Φ[1][4]*Φ[2][4]*Φ[3][4]     
    a+b+c+d           
end

function eval_spline( ranges, C::AbstractArray{T,d}, x::SVector{d,U}) where d where U where T 

    
    a = SVector( (e[1] for e in ranges)... )
    b = SVector( (e[2] for e in ranges)... )
    n = SVector( (e[3] for e in ranges)... )
    
    Tf = eltype(a)

    δ = (b.-a)./(n.-convert(Tf,1))
    i = div.( (x.-a), δ)
    i = max.(min.(i, n.-2), 0)
    λ = (x.-(a .+ δ.*i))./δ

    # i_ = floor.(Int, i) .+ 1
    i_ = unsafe_trunc.(Int, i) .+ 1
    
    Φ = tuple( (cPhi(e) for e in λ)... )

    M = mextract(C, tuple( (ii for ii in i_)...), Val(4))

    # return sum(M) + sum(Φ[1])+ sum(Φ[2]) + sum(Φ[3]) 
    # reduce_tensors(M,Φ)
    @inline reduce_tensors(M,Φ)

end