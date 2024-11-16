# This starts making a difference when grids have more than 200000 points...

using KernelAbstractions
import KernelAbstractions.Extras: @unroll

import KernelAbstractions: get_backend

get_backend(g::GArray) = get_backend(g.data)

function F!(r, model, x, φ, engine)

    @kernel function FF_(r, @Const(dm), @Const(x), @Const(φ))

        n = @index(Global, Linear)
    
        # TODO: special function here
        s = dm.grid[n+im]
    
        xx = x[n]
        
        rr = Dolo.F(dm, s, xx, φ)
    
        r[n] = rr
    
    end

    
    sz = size(model.grid)
    
    fun_cpu = FF_(engine, 1000)

    fun_cpu(r, model, x, φ; ndrange=sz)
    synchronize(engine)


end



function dF_1!(out, model, controls::GArray, φ::Union{GArray, DFun}, engine)

    @kernel function FF_(r,@Const(dm), @Const(x),@Const(φ) )


        n = @index(Global, Linear)
    
        # TODO: special function here
        s = dm.grid[n+im]
    
        xx = x[n]
    
   
        rr = ForwardDiff.jacobian(u->F(model, s, u, φ), xx)
    
        r[n] = rr

    end

    fun_cpu = FF_(engine)

    sz = size(model.grid)

    res = fun_cpu(out, model, controls, φ; ndrange=sz)
    # wait(res)

end    #### no alloc
    
function dF_2!(out, dmodel, controls::GArray, φ::DFun, engine)

    @kernel function FF_(L, @Const(dm), @Const(x),@Const(φ) )


        n = @index(Global, Linear)
    
        s = dm.grid[n+im]
    
        xx = x[n]

        # r0 = sum( w*arbitrage(dmodel,s,xx,S,φ(S)) for (w,S) in τ(dmodel, s, xx) )
        r0 = F(dmodel, s, xx, φ)
        r_F = ForwardDiff.jacobian(
            r->complementarities(dmodel.model, s,xx,r),
            r0
        )

        tt = tuple(
            (
                (;
                    F_x=w*r_F*ForwardDiff.jacobian(u->Dolo.arbitrage(dmodel,s,xx,S,u), φ(S)),
                    S=S
                )
            for (w,S) in Dolo.τ(dmodel, s, xx)
            )
        ...)

        L.D[n] = tt

    end

    fun = FF_(engine)
    sz = size(dmodel.grid)

    fun(out, dmodel, controls, φ; ndrange=sz)

    nothing



end


# function dF_2!(L, dmodel, xx::GArray, φ::DFun, ::CPU)

#     # for (n,(s,x)) in enumerate(zip(Dolo.enum(dmodel.grid), xx))
#     Threads.@threads for n=1:length(dmodel.grid)

#         i,j = Dolo.from_linear(dmodel.grid, n)
            
#         s_ = dmodel.grid[i,j]
#         s = QP((i,j), s_)
#         x = xx[n]

#         L.D[n] = tuple(
#                 (
#                     (;
#                         F_x=w*ForwardDiff.jacobian(u->Dolo.arbitrage(dmodel,s,x,S,u), φ(S)),
#                         S=S
#                     )
#                 for (w,S) in Dolo.τ(dmodel, s, x)
#                 )
#         ...)

#     end

#     nothing


# end

function mul!(dr, L2::Dolo.LL, x, engine)
 
    D = L2.D
    dφ = L2.φ

    fit!(dφ, x)

    @kernel function tm_(dr, @Const(D), @Const(dφ))
        n = @index(Global, Linear)
        t0 = zero(eltype(dr))
        @unroll for k=1:length(D[n])    # @unroll provides  a 20% gain
            (;F_x, S) = D[n][k]      # profile: slow ?
            t0 += F_x*dφ(S) # TODO : add back complementarity
        end
        dr[n] = t0
    end

    fun = tm_(engine)

    K = length(D)
    res = fun(dr, D, dφ; ndrange=K)

    synchronize(engine)
    # wait(res)
    nothing

end


function \(J, L::Dolo.LL) 


    @kernel function tm_(@Const(J), L)
        n = @index(Global)
        
        elt = L.D[n]
        M = J.data[n]
        tt = tuple((
            (;F_x=M\d.F_x, S=d.S) 
            for d in elt
        )...)
        L.D[n] = tt
    
    end

    engine = get_backend(J.data)
    fun = tm_(engine)

    K = length(J.data)

    nL = deepcopy(L)
    res = fun(J, nL; ndrange=K)

    return nL

end
