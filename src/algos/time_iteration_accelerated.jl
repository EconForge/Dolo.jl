# This starts making a difference when grids have more than 200000 points...

using KernelAbstractions



# This is for a CGrid model only 
function F!(r, model, x, φ, ::CPU)

    @kernel function FF_(r, @Const(model), @Const(x), @Const(φ))

        i = @index(Global, Linear)

        s_ = model.grid[i]
        s = QP((i,), s_)
        xx = x[i]
        
        rr = sum(
            w*Dolo.arbitrage(model,s,xx,S,φ(S)) 
            for (w,S) in Dolo.τ(model, s, xx)
        )      
        
        r[i] = rr

    end

    fun_cpu = FF_(CPU())


    p = length(model.grid)
    
    # p,q = size(x)

    res = fun_cpu(r, model, x, φ; ndrange=(p,))
    wait(res)

end

# # This is for a PGrid model only 
# function F!(r, model, x, φ, ::CPU)

#     @kernel function FF_(r, @Const(model), @Const(x), @Const(φ))

#         c = @index(Global, Cartesian)

#         i,j = c.I

#         s_ = model.grid[i,j]
#         s = QP((i,j), s_)
#         xx = x[i,j]
        
#         rr = sum(
#             w*Dolo.arbitrage(model,s,xx,S,φ(S)) 
#             for (w,S) in Dolo.τ(model, s, xx)
#         )      
        
#         r[i,j] = rr

#     end

#     fun_cpu = FF_(CPU())


#     p = length(model.grid.g1)
#     q = length(model.grid.g2)
#     # p,q = size(x)

#     res = fun_cpu(r, model, x, φ; ndrange=(p,q))
#     wait(res)

# end



function dF_1!(out, model, controls::GArray, φ::Union{GArray, DFun}, ::CPU)

    @kernel function FF_(r,@Const(model), @Const(x),@Const(φ) )

        c = @index(Global, Cartesian)

        i,j = c.I

        s_ = model.grid[i,j]
        s = QP((i,j), s_)
        xx = x[i,j]
        
        rr = ForwardDiff.jacobian(u->F(model, s, u, φ), xx)
        
        r[i,j] = rr

    end

    fun_cpu = FF_(CPU())

    p = length(model.grid.g1)
    q = length(model.grid.g2)
    # p,q = size(out)

    res = fun_cpu(out, model, controls, φ; ndrange=(p,q))
    wait(res)

end    #### no alloc
    
