import Base: *

struct LL{G,D_,F}
    grid::G
    D::D_
    φ::F
end


const LF = LL

import Base: /,*
import LinearAlgebra: mul!

function mul!(L::LL, x::Number)

    for n in 1:length(L.D)
        L.D[n]  =tuple((
            (;F_x=x*d.F_x, S=d.S) 
            for d in L.D[n]
        )...)
    end
end


using LoopVectorization

function mul!(dr, L2::LL, x)
    
    D = L2.D
    # N,K = size(M_ij)
    dφ = L2.φ

    fit!(dφ, x)
    for n=1:length(D)
    # Threads.@threads for n=1:length(D)
        t0 = zero(eltype(dr))
        for k=1:length(D[n])
            (;F_x, S) = D[n][k]
            t0 += F_x*dφ(S) # TODO : add back complementarity
        end
        dr[n] = t0
    end
end

function *(L::LL, x0)
    r = deepcopy(x0)
    r.data .= r.data .* 0.0
    mul!(r, L, x0)
    return r
end

function \(J, L::Dolo.LL) 
    Dolo.LL(
        L.grid,
        [ tuple((
            (;F_x=J.data[n]\d.F_x, S=d.S) 
            for d in L.D[n]
        )...)
        for n in 1:length(L.D) ],
        deepcopy(L.φ)
    )
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


function dF_2!(L, model, xx::GArray, φ::DFun)

    for (n,(s,x)) in enumerate(zip(Dolo.enum(dmodel.grid), xx))
        L.D[n] = tuple(
                (
                    (;
                        F_x=w*ForwardDiff.jacobian(u->Dolo.arbitrage(dmodel,s,x,S,u), φ(S)),
                        S=S
                    )
                for (w,S) in Dolo.τ(dmodel, s, x)
                )
        ...)
    end
    
    LL(
        dmodel.grid,
        res,
        deepcopy(φ)
    )

end

# # this is for SGrid \times CGrid
# function dF_2(model, x::GArray, φ::DFun )

#     res = []
#     for (s,x) in zip(enum(model.grid), x)
#         l = []
#         for (w, S) in τ(model, s, x)
#             el = w*ForwardDiff.jacobian(u->arbitrage(model,s,x,S,u), φ(S))
#             # push!(l, (el, S.loc[2]))
#             push!(l, (el, S.loc))
#         end
#         push!(res, l)
#     end
#     N = length(res)
#     J = length(res[1])
#     tt = res[1][1]
#     M_ij = Array{typeof(tt[1])}(undef,N,J)
#     S_ij = Array{typeof(tt[2])}(undef,N,J)
#     for n=1:N
#         for j=1:J
#             M_ij[n,j] = res[n][j][1]
#             S_ij[n,j] = res[n][j][2]
#         end
#     end

#     L = LF(model.grid, M_ij, S_ij, deepcopy(φ))
#     return L
# end

# function dF_2!(L, model, xx::GArray, φ::DFun )
#     for (n,(s,x)) in enumerate(zip(enum(model.grid), xx))
#         for (j,(w, S)) in enumerate(τ(model, s, x))
#             el = w*ForwardDiff.jacobian(u->arbitrage(model,s,x,S,u), φ(S))
#             L.M_ij[n,j] = el
#             L.S_ij[n,j] = S.loc
#         end
#     end
# end




# function mul!(dr, L2::LF, x)
#     (;M_ij, S_ij) = L2
#     N,K = size(M_ij)
#     dφ = L2.φ
#     fit!(L2.φ, x)
#     for n=1:N
#         t0 = dr[n]*0.0
#         for k=1:K
#             F_x = M_ij[n,k]
#             S = S_ij[n,k]
#             t0 += F_x*dφ(S)
#         end
#         dr[n] = t0
#     end
# end


function neumann(L2::LF, r0; K=1000, τ_η=1e-10)
    
    dx = deepcopy(r0)
    du = deepcopy(r0)
    dv = deepcopy(r0)

    mem = (;du, dv)

    neumann!(dx, L2, r0, mem; K=K, τ_η=τ_η)

    return dx

end

function neumann!(dx, L2::LF, r0, mem=(;du=deepcopy(r0), dv=deepcopy(r0)); K=1000, τ_η=1e-10)
    
    (; du, dv) = mem
    
    dx.data .= r0.data
    du.data .= r0.data
    dv.data .= r0.data

    for k=1:K
        mul!(du, L2, dv)
        dx.data .+= du.data
        η = norm(du)
        if η<τ_η
            break
        end
        du,dv=dv,du
    end
    # return dx
end


# function *(L::LF, x0)
#     r = deepcopy(x0)
#     r.data .= r.data .* 0.0
#     mul!(r, L, x0)
#     return r
# end

# \(J, L::LF) = LF(L.grid, J.data .\ L.M_ij, L.S_ij, L.φ)
# \(L::LF, r_) = solve(L, r_)

using LinearMaps
using BlockDiagonals
convert(::Type{BlockDiagonal}, ga::GArray{G, Vector{T}}) where T where G = BlockDiagonal(ga.data)


convert(::Type{Matrix}, L::LF) = convert(Matrix, convert(LinearMap,L))

convert(::Type{LinearMap}, L::LF) = let
    # elt = eltype(L.D[1][1].F_x)
    elt = typeof(L.D[1][1].F_x)
    p,q = size(elt)
    N = length(L.grid)
    typ = SVector{q, Float64}
    fun = u->(ravel(L*GArray(L.grid, reinterpret(typ, u))))
    LinearMap(
        fun,
        p*N,
        q*N
    )
end



