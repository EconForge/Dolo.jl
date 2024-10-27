import Base: *

struct LL{G,D_,F}
    grid::G
    D::D_
    φ::F
end


const LF = LL

import Base: /,*
import LinearAlgebra: mul!

# function mul!(L::LL, x::Number)

#     for n in 1:length(L.D)
#         L.D[n]  =tuple((
#             (;F_x=x*d.F_x, S=d.S) 
#             for d in L.D[n]
#         )...)
#     end
# end

function mul!(L::LL, x::Number)

    map!(
        e->tuple((
            (;F_x=x*d.F_x, S=d.S) 
            for d in e
        )...),
        L.D,
        L.D
    )
    # for n in 1:length(L.D)
    #     L.D[n]  =tuple((
    #         (;F_x=x*d.F_x, S=d.S) 
    #         for d in L.D[n]
    #     )...)
    # end
end

function mul!(dr, L2::LF, x)
    engine=get_backend(L2.D)
    mul!(dr, L2::LF, x, engine)
end


function mul!(dr, L2::LL, x, ::Nothing)
 
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
    r = duplicate(x0) ###### BUG: this doesn't wokr on the GPU
    Tf = getprecision(x0)
    r.data .*= convert(Tf, 0.0)
    mul!(r, L, x0)
    r
    
end

# this takes 0.2 s !
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

function dF_2!(L, dmodel, xx::GArray, φ::DFun)
    dF_2!(L, dmodel, xx, φ, CPU())
end

function dF_2!(L, dmodel, xx::GArray, φ::DFun, ::Nothing)

    for (n,(s,x)) in enumerate(zip(Dolo.enum(dmodel.grid), xx))
        (i,j) = Dolo.from_linear(dmodel.grid, n)
        s_ = dmodel.grid[i,j]
        s = QP((i,j), s_)

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





function neumann(L2::LF, r0; K=1000, τ_η=1e-10)
    
    # TODO
    # dx = deepcopy(r0)
    # du = deepcopy(r0)
    # dv = deepcopy(r0)

    dx = duplicate(r0)
    du = duplicate(r0)
    dv = duplicate(r0)

    mem = (;du, dv)

    neumann!(dx, L2, r0, mem; K=K, τ_η=τ_η)

    return dx

end

function neumann!(dx, L2::LF, r0, mem=(;du=duplicate(r0), dv=duplicate(r0)); K=1000, τ_η=1e-10)
    
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
    
    nothing

end



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



