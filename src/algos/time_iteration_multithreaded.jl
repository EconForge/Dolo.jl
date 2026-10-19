## In theory this is brilliant !
## In practice with Julia 1.12, Threads.@threds works equally well

# This starts making a difference when grids have more than 200000 points...

# using KernelAbstractions
# import KernelAbstractions.Extras: @unroll
# import KernelAbstractions: get_backend
# get_backend(g::GArray) = get_backend(g.data)

import Base: Threads

struct CPUS end # multithreaded engine

function F!(r, dm, x, φ, engine::CPUS)

    Threads.@threads for n=1:length(dm.grid)
    
        # TODO: special function here
        s = dm.grid[n+im]
        xx = x[n]
        rr = Dolo.F(dm, s, xx, φ)
        r[n] = rr
    
    end

end



function dF_1!(out, dmodel, controls::GArray, φ::Union{GArray, DFun}, engine::CPUS)

    Threads.@threads for n=1:length(dmodel.grid)

        s = dmodel.grid[n+im]
        xx = controls[n]
        rr = dF_1(dmodel, s, xx, φ)
        # rr = ForwardDiff.jacobian(u->F(dmodel, s, u, φ), xx)
        out[n] = rr

    end


end    #### no alloc
    
function dF_2!(L, dmodel, controls::GArray, φ::DFun, engine::CPUS)

    x = controls
    Threads.@threads for n=1:length(dmodel.grid)

        s = dmodel.grid[n+im]
        xx = x[n]
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

end


function mul!(dr::GVector, L2::LL, x::GVector, engine::CPUS)
 

    D = L2.D
    dφ = L2.φ

    fit!(dφ, x)

    Threads.@threads for n=1:length(D)

        t0 = zero(eltype(dr))
        # TODO add unroll
        for k=1:length(D[n])
            (;F_x, S) = D[n][k]
            t0 += F_x*dφ(S) # TODO : add back complementarity
        end
        dr[n] = t0
    end

end

# function \(J, L::Dolo.LL) 

#     Threads.@threads for n=1:length(D)

#         elt = L.D[n]
#         M = J.data[n]
#         tt = tuple((
#             (;F_x=M\d.F_x, S=d.S) 
#             for d in elt
#         )...)
#         L.D[n] = tt
        
#     end

# end

function ldiv!(J, L::Dolo.LL, engine::CPUS)
    Threads.@threads for n in 1:length(L.D)
        E = L.D[n]
        L.D[n]=  tuple((
            (;F_x=J.data[n]\d.F_x, S=d.S) 
            for d in E
        )...)
    end
end


function mul!(L::LL, x::Number, engine::CPUS)

    Threads.@threads for n=1:length(L.D)
        e = L.D[n]
        L.D[n] = tuple((
            (;F_x=x*d.F_x, S=d.S) 
            for d in e
        )...)
    end

end
