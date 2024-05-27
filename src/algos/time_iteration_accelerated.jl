# This starts making a difference when grids have more than 200000 points...

using KernelAbstractions



# # This is for a CGrid model only 
# function F!(r, model, x, φ, ::CPU)

#     @kernel function FF_(r, @Const(model), @Const(x), @Const(φ))

#         i = @index(Global, Linear)

#         s_ = model.grid[i]
#         s = QP((i,), s_)
#         xx = x[i]
        
#         rr = sum(
#             w*Dolo.arbitrage(model,s,xx,S,φ(S)) 
#             for (w,S) in Dolo.τ(model, s, xx)
#         )      
        
#         r[i] = rr

#     end

#     fun_cpu = FF_(CPU())


#     p = length(model.grid)
    
#     # p,q = size(x)

#     res = fun_cpu(r, model, x, φ; ndrange=(p,))
#     wait(res)

# end

# This is for a PGrid model only 
function F!(r, model, x, φ, engine)

    @kernel function FF_(r, @Const(dm), @Const(x), @Const(φ))

        n = @index(Global, Linear)
        ind = @index(Global, Cartesian)
        (i,j) = ind.I
    
        # k = length(dm.grid.grids[1])
        # i_ = (n-1)÷k
        # j_ = (n-1) - (i)*κ
        # i = i_+1
        # j = j_+1
        
        # TODO: special function here
        s_ = dm.grid[n]
        s = Dolo.QP((i,j), s_)
    
        xx = x[n]
    
        # (i,j) = @inline Dolo.from_linear(model.grid, n)    
    
        rr = Dolo.F(dm, s, xx, φ)
    
        r[n] = rr
    
    end

    fun_cpu = FF_(engine)


    # p = length(model.grid.grids[1])
    # q = length(model.grid.grids[2])
    # p,q = size(x)
    if typeof(model.grid)<:CGrid
        p = model.grid.ranges[1][3]
        q = model.grid.ranges[2][3]
    else
        p = length(model.grid.grids[1])
        q = length(model.grid.grids[2])
    end
    
    res = fun_cpu(r, model, x, φ; ndrange=(p,q))
    # wait(res)

end



function dF_1!(out, model, controls::GArray, φ::Union{GArray, DFun}, engine)

    @kernel function FF_(r,@Const(dm), @Const(x),@Const(φ) )


        n = @index(Global, Linear)
        ind = @index(Global, Cartesian)
        (i,j) = ind.I
    
        # k = length(dm.grid.grids[1])
        # i_ = (n-1)÷k
        # j_ = (n-1) - (i)*κ
        # i = i_+1
        # j = j_+1
        
        # TODO: special function here
        s_ = dm.grid[n]
        s = Dolo.QP((i,j), s_)
    
        xx = x[n]
    
        # (i,j) = @inline Dolo.from_linear(model.grid, n)    
    
        rr = ForwardDiff.jacobian(u->F(model, s, u, φ), xx)
    
        r[n] = rr

    end

    fun_cpu = FF_(engine)

    if typeof(model.grid)<:CGrid
        p = model.grid.ranges[1][3]
        q = model.grid.ranges[2][3]
    else
        p = length(model.grid.grids[1])
        q = length(model.grid.grids[2])
    end
    # p,q = size(out)

    res = fun_cpu(out, model, controls, φ; ndrange=(p,q))
    # wait(res)

end    #### no alloc
    

function dF_2!(L, dmodel, xx::GArray, φ::DFun, ::CPU)

    # for (n,(s,x)) in enumerate(zip(Dolo.enum(dmodel.grid), xx))
    Threads.@threads for n=1:length(dmodel.grid)

        i,j = Dolo.from_linear(dmodel.grid, n)
            
        s_ = dmodel.grid[i,j]
        s = QP((i,j), s_)
        x = xx[n]

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

    nothing


end

using KernelAbstractions
import KernelAbstractions.Extras: @unroll

function mul!(dr, L2::Dolo.LL, x, ::CPU)
 
    D = L2.D
    dφ = L2.φ

    fit!(dφ, x)


    @kernel function tm_(dr, @Const(L2), @Const(x))
        n = @index(Global)
        t0 = zero(eltype(dr))
        @unroll for k=1:length(D[n])    # @unroll provides  a 20% gain
            (;F_x, S) = D[n][k]      # profile: slow ?
            t0 += F_x*dφ(S) # TODO : add back complementarity
        end
        dr[n] = t0
    end

    fun = tm_(CPU())

    K = length(D)
    res = fun(dr, L2, x; ndrange=K)
    # wait(res)
    res

end
