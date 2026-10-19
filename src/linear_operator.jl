import Base: *
using LinearMaps
using BlockDiagonals
using LinearAlgebra: UniformScaling


struct LL{G,D_,F}
    grid::G
    D::D_
    φ::F
end

struct I_L{_L}
    L::_L
end

-(::UniformScaling, L::LL) = I_L(L)
*(A::I_L, b) = b-A.L*b

convert(::Type{LinearMap}, IL::I_L) = I - convert(LinearMap,IL.L)


const LF = LL

import Base: /,*
import LinearAlgebra: mul!, ldiv!


function mul!(L::LL, x::Number)

    map!(
        e->tuple((
            (;F_x=x*d.F_x, S=d.S) 
            for d in e
        )...),
        L.D,
        L.D
    )

end

function mul!(L::LL, x::Number, ::Nothing)

    mul!(L, x)

end

function mul!(dr, L2::LF, x)
    engine=get_backend(L2.D)
    mul!(dr, L2::LF, x, engine)
end


function mul!(dr, L2::LL, x::GVector, ::Nothing)

    D = L2.D
    dφ = L2.φ

    fit!(dφ, x)
    for n=1:length(D)
        t0 = zero(eltype(dr))
        for k=1:length(D[n])
            (;F_x, S) = D[n][k]      # profile: slow ?
            t0 += F_x*dφ(S)          # TODO : add back complementarity
        end
        dr[n] = t0
    end
end


function *(L::LL, x0)

    # BUG: with deepcopy, doesn' t work
    # TODO: why do we even copy here?
    # println("did you call me?")
    r = duplicate(x0) ###### BUG: this doesn't wokr on the GPU
    Tf = getprecision(x0)
    r.data .*= convert(Tf, 0.0)
    mul!(r, L, x0)
    r
    
end

function *(L::LL, v::AbstractArray{T}) where T<:Union{Float32, Float64}

    # println("did you call me?")
    r = copy(v) ###### BUG: this doesn't wokr on the GPU
    mul!(r, L, v)
    return r
    
end


function mul!(r_::AbstractArray{T}, L::LL, v_::AbstractArray{T}) where T<:Union{Float32, Float64}

    mul!(r_, L, v_, nothing)

end

function mul!(r_::AbstractArray{T}, L::LL, v_::AbstractArray{T}, t_engine) where T<:Union{Float32, Float64}

    elt = typeof(L.D[1][1].F_x)
    Tf = eltype(elt)
    p,q = size(elt)
    N = length(L.grid)
    typ = SVector{q, Tf}

    r = GArray(L.grid, reinterpret(typ, r_))
    x = GArray(L.grid, reinterpret(typ, v_))

    mul!(r, L, x, t_engine)

    r_
    
end


function ldiv!(J, L::Dolo.LL)
    for n in 1:length(L.D)
        E = L.D[n]
        L.D[n]=  tuple((
            (;F_x=J.data[n]\d.F_x, S=d.S) 
            for d in E
        )...)
    end
end

ldiv!(J, L::Dolo.LL, ::Nothing) = ldiv!(J, L)

# this takes 0.2 s !
function \(J, L::Dolo.LL) 
    R = deepcopy(L)
    ldiv!(J, R)
    return R
end

   
function dF_2(dmodel, xx, φ)

    res = [
            let
                r_F = ForwardDiff.jacobian(
                    r->complementarities(dmodel.model, s,x,r),
                    sum( w*arbitrage(dmodel,s,x,S,φ(S)) for (w,S) in τ(dmodel, s, x) ),
                )
                tuple(
                    (
                            (;
                                F_x=w*r_F*ForwardDiff.jacobian(u->Dolo.arbitrage(dmodel,s,x,S,u), φ(S)),
                                S=S
                            )
                    for (w,S) in Dolo.τ(dmodel, s, x)
                    )
                ...)
            end
           for (s,x) in zip(Dolo.enum(dmodel.grid), xx) 
        ]
    
    LL(
        dmodel.grid,
        res,
        deepcopy(φ)
    )

end

function dF_2!(L, dmodel, xx::GArray, φ::DFun)
    dF_2!(L, dmodel, xx, φ, CPUS())
end

function dF_2!(L, dmodel, xx::GArray, φ::DFun, ::Nothing)


    for n=1:length(xx)

        s = dmodel.grid[n+im]
        x = xx[n]

        r_F = ForwardDiff.jacobian(
            r->complementarities(dmodel.model, s,x,r),
            sum( w*arbitrage(dmodel,s,x,S,φ(S)) for (w,S) in τ(dmodel, s, x) ),
        )

        L.D[n] = tuple(
                (
                    (;
                        F_x=w*r_F*ForwardDiff.jacobian(u->Dolo.arbitrage(dmodel,s,x,S,u), φ(S)),
                        S=S
                    )
                for (w,S) in Dolo.τ(dmodel, s, x)
                )
        ...)
    end
    nothing


end


function neumann(L2::LF, r0; K=1000, τ_η=1e-10, t_engine=nothing)
    
    # TODO
    # dx = deepcopy(r0)
    # du = deepcopy(r0)
    # dv = deepcopy(r0)

    dx = duplicate(r0)
    du = duplicate(r0)
    dv = duplicate(r0)

    mem = (;du, dv)

    neumann!(dx, L2, r0, t_engine;
        mem = mem, K=K, τ_η=τ_η)

    return dx

end

# computes dx = (I-L)^(-1) r0 = Σ_k L^k r0, returns the number of terms
function neumann!(dx, L2::LF, r0, t_engine=nothing; mem=(;du=duplicate(r0), dv=duplicate(r0)), K=1000, τ_η=1e-10, tmr=NO_TIMER)
    
    (; du, dv) = mem
    
    dx.data .= r0.data
    du.data .= r0.data
    dv.data .= r0.data

    k = 0
    while k < K
        k += 1
        @timeit tmr "eval L" mul!(du, L2, dv, t_engine)
        dx.data .+= du.data
        η = norm(du)
        if η<τ_η
            break
        end
        du,dv=dv,du
    end
    
    return k

end

convert(::Type{BlockDiagonal}, ga::GArray{G, Vector{T}}) where T where G = BlockDiagonal(ga.data)


convert(::Type{Matrix}, L::LF) = convert(Matrix, convert(LinearMap,L))


convert(::Type{LinearMap}, L::LF, engine=nothing) = let
    elt = typeof(L.D[1][1].F_x)
    Tf = eltype(elt)
    p,q = size(elt)
    N = length(L.grid)
    typ = SVector{q, Tf}
    LinearMap(
        (v,u)-> mul!(v,L,u, engine), # AWESOME
        p*N,
        q*N;
        ismutating=true
    )
end

