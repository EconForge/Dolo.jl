# function F(model, s, x::SVector, φ::Policy)
#     tot = 0.0
#     for (w,S) in τ(model, s, x)
#         X = φ(S)
#        tot +=  w*arbitrage(model,s,x,S,X) 
#     end
#     tot
# end

# F(model, s, x::SVector, φ::Policy) = let
#     r = sum(
#          w*arbitrage(model,s,x,S,φ(S)) 
#          for (w,S) in τ(model, s, x)
#     )
#     # r
#     r = complementarities(model.model, s,x,r)
#     r
# end



function F(dmodel::ADModel, s::QP, x::SVector{d,T}, φ::Union{Policy, GArray, DFun}) where d where T

    r = zero(SVector{d,T})

    for (w,S) in τ(dmodel, s, x)
        # return w,S
        r += w*arbitrage(dmodel,s,x,S,φ(S)) 
    end

    # TODO: why does the following allocate ?
    # strange: if reloaded it doesn't allocate anymore
    # r += sum(
    #      w*arbitrage(model,s,x,S,φ(S)) 
    #      for (w,S) in τ(model, s, x)
    # )

    r::SVector{d,T}
    r = complementarities(dmodel.model, s,x,r)
    r
end


F(model::ADModel, controls::GArray, φ::Union{GArray, DFun}) =
    GArray(
        model.grid,
        [
            F(model,s,x,φ)
            for (s,x) in zip(enum(model.grid), controls)
        ],
    )

function F!(out, model::ADModel, controls::XT, φ::Union{GArray, DFun}) where XT<:GArray{PG, Vector{SVector{d,T}}} where PG where d where T
    
    for s in enum(model.grid)
        x = controls[s.loc...]
        out[s.loc...] = F(model,s,x,φ)
    end

    return nothing
end


function F!(out, model::ADModel, controls::XT, φ::Union{GArray, DFun}, ::Nothing) where XT<:GArray{PG, Vector{SVector{d,T}}} where PG where d where T
    
    for s in enum(model.grid)
        x = controls[s.loc...]
        out[s.loc...] = F(model,s,x,φ)
    end

    return nothing
end


## no alloc
dF_1(model::ADModel, s, x::SVector, φ) = ForwardDiff.jacobian(u->F(model, s, u, φ), x)

dF_1(model::ADModel, controls::GArray, φ::Union{GArray, DFun}) =
    GArray(    # this shouldn't be needed
        model.grid,
        [
            dF_1(model, s, x, φ)
            for (s,x) in zip(enum(model.grid), controls) 
        ]
    )


## no alloc

function dF_1!(out, model::ADModel, controls::GArray, φ::Union{GArray, DFun}, ::Nothing)
    dF_1!(out, model, controls::GArray, φ::Union{GArray, DFun})
end

function dF_1!(out, model::ADModel, controls::GArray, φ::Union{GArray, DFun})

    i = 0
    for s in enum(model.grid)
        i+=1
        x = controls[i]
        # x = controls[s.loc...]
        # out[s.loc...] = dF_1(model, s, x, φ)
        out.data[i] = dF_1(model, s, x, φ)
    end
    return nothing
end    #### no alloc
    


# warning: the versions of dF_2 don't have complementarities

dF_2(model::ADModel, s, x::SVector, φ::GArray, dφ::GArray) = 
    sum(
            w*ForwardDiff.jacobian(u->arbitrage(model,s,x,S,u), φ(S))* dφ(S)
            for (w, S) in τ(model, s, x)
    )   ### no alloc


dF_2(model::ADModel, controls::GArray, φ::GArray, dφ::GArray) =
    GArray(
        model.grid,
        [(dF_2(model,s,x,φ,dφ)) for (s,x) in zip(enum(model.grid), controls) ],
    )

function dF_2!(out::GArray, model::ADModel, controls::GArray, φ::GArray, dφ::GArray)
    for (n, (s,x)) in enumerate(zip(enum(model.grid), controls))
        out[n] = dF_2(model, s, x, φ, dφ)
    end
end   


using LinearMaps
using Adapt

function time_iteration_workspace(dmodel; interp_mode=:linear, improve=false, dest=Array)

    T = eltype(dmodel)

    x0 = (Dolo.initial_guess(dmodel))
    x1 = deepcopy(x0)
    x2 = deepcopy(x0)
    r0 = deepcopy(x0)
    dx = deepcopy(x0)
    N = length(dx)
    n = length(dx.data[1])
    J = GArray(
        dmodel.grid,
        zeros(SMatrix{n,n,T,n*n}, N)
    )
    vars = variables(dmodel.model.controls)
    φ = DFun(dmodel.model.states, x0, vars; interp_mode=interp_mode)

    if improve
        L = Dolo.dF_2(dmodel, x1, φ)
        tt = (;x0, x1, x2, r0, dx, J, L, φ)
    else
        tt = (;x0, x1, x2, r0, dx, J, φ)
    end

    return adapt(dest, tt)

end

function newton_workspace(dmodel::ADModel; interp_mode=:linear)

    
    res = time_iteration_workspace(dmodel; interp_mode=interp_mode)
    T =  Dolo.dF_2(dmodel, res.x0, res.φ)
    res = merge(res, (;T=T,memn=(;du=deepcopy(res.x0), dv=deepcopy(res.x0))))
    return res
end

function time_iteration(model::YModel; kwargs...)
    discr_options = get(kwargs, :discretization, Dict())
    interp_mode = get(kwargs, :interpolation, :cubic)
    improve = get(kwargs, :improve, false)
    dmodel = discretize(model, discr_options)
    wksp = time_iteration_workspace(dmodel; interp_mode=interp_mode, improve=improve)
    kwargs2 = pairs(NamedTuple( k=>v for (k,v) in kwargs if !(k in (:discretization, :interpolation))))
    time_iteration(dmodel, wksp; kwargs2...)
end

function time_iteration(model::DYModel,
    workspace=time_iteration_workspace(model);
    T=500,
    improve_K=1000,
    improve_wait=10,
    tol_ε=1e-8,
    tol_η=1e-6,
    newton_steps=10,
    mbsteps = 5,
    verbose=true,
    trace=false,
    improve=false,
    engine=:none
)


    log = IterationLog(
            it = ("n", Int),
            err =  ("ϵₙ=|F(xₙ,xₙ)|", Float64),
            sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
            lam = ("λₙ=ηₙ/ηₙ₋₁",Float64),
            elapsed = ("Time", Float64),
            k_newton = ("Newton", Int64),
            k = ("Backsteps", Int64)
        )    
    initialize(log, verbose=verbose; message="Time Iteration")

    # mem = typeof(workspace) <: Nothing ? time_iteration_workspace(model) : workspace

    (;x0, x1, x2, r0, dx, J, φ) = workspace

    
    local η_0 = NaN
    local η
    convergence = false
    iterations = T
    
    Tf = eltype(eltype(x0))

    lam = convert(Tf, 0.5)

    if engine in (:cpu, :gpu)
        t_engine = get_backend(workspace.x0)
    else
        t_engine = nothing
    end

    ti_trace = trace ? IterationTrace(typeof(φ)[]) : nothing

    if improve
        J_2 = workspace.L
    end
    local k_newton_ = 0
    local k_ = 0
    local ε

    for t=1:T
        
        t1 = time_ns()

        Dolo.fit!(φ, x0)

        trace && push!(ti_trace.data, deepcopy(φ))

        F!(r0, model, x0, φ, t_engine)

        # r0 = F(model, x0, φ)

        ε = norm(r0)

        if ε<tol_ε
            convergence = true
            iterations = t
            break
        end

        # solve u->F(u,x0) 
        # result in x1
        # J and r0 are modified

        x1.data .= x0.data

 
        for k_newton=1:newton_steps
            
            #newton steaps

            F!(r0, model, x1,  φ, t_engine)

            ε_n = norm(r0)
            if ε_n<tol_ε
                iterations = t
                break
            end

            dF_1!(J, model, x1,  φ, t_engine)
            
            # TODO: fix regression or revert to Julia 1.10
            @inline dx.data .= J.data .\ r0.data
            # for k =1:length(dx.data)
            #     dx.data[k] = J.data[k][1,1] \ r0.data[k]
            # end


            for k=0:mbsteps

                # fails on GPU
                x2.data .= x1.data .- dx.data .* convert(Tf,lam^k)
                F!(r0, model, x2,  φ, t_engine)
                ε_b = norm(r0)
                if ε_b<ε_n
                    iterations = t
                    k_newton_ = k_newton
                    k_ = k
                    break
                end
                
            end

            x1.data .= x2.data

        end


        η = distance(x0, x1)
        gain = η/η_0

        # verbose ? println("$t: $ε : $η: ") : nothing

        elapsed = time_ns() - t1
        elapsed /= 1000000000

        verbose ? append!(log; verbose=verbose, it=t-1, err=ε, sa=η_0, lam=gain, elapsed=elapsed, k_newton=k_newton_, k=k_) : nothing

        if η < tol_η
            iterations = t
            break
        end
        η_0 = η


        if improve && t>improve_wait
            # x = T(x)
            # x1 = T(x0)
            # x - x1 = -T'(x) (x - x0)
            # x = x1 - T' (x - x0)
            # x = (I-T')\(x1 - T' x0)

            J_1 = J

            Dolo.dF_1!(J_1, model, x1, φ, t_engine)

            Dolo.dF_2!(J_2, model, x1, φ, t_engine)
            # J_2 = Dolo.dF_2(model, x1, φ)
            
            mul!(J_2, convert(Tf,-1.0))
            
            # J_2.M_ij[:] *= -1.0
            Tp = J_1 \ J_2
            
            d = (x1-x0)
            
            # return (J_1, J_2, d)
            rhs =  neumann(Tp, d; K=improve_K, t_engine=t_engine)

            x0.data .+= rhs.data

        else
            x0.data .= x1.data
        end

    end

    verbose ?  finalize(log, verbose=verbose) : nothing

    TimeIterationResult(
        φ,
        iterations,
        tol_η,
        η,
        log,
        ti_trace
    )
end


function newton(model::ADModel, workspace=newton_workspace(model);
    K=20, tol_ε=1e-8, tol_η=1e-6, verbose=false, improve=false, interp_mode=:cubic,
    engine=:none
    )

    # mem = typeof(workspace) <: Nothing ? time_iteration_workspace(model) : workspace

    (;x0, x1, x2, dx, r0, J, φ, T, memn) = workspace

    if engine in (:cpu, :gpu)
        t_engine = get_backend(workspace.x0)
    else
        t_engine = nothing
    end


    for t=1:K
        
        Dolo.fit!(φ, x0)

        F!(r0, model, x0, φ, t_engine)

        ε = norm(r0)
        verbose ? println("$t: $ε") : nothing

        if ε<tol_ε
            return (;message="Solution found", solution=x0, n_iterations=t-1, dr=φ)
        end

        x1.data .= x0.data


        dF_1!(J, model, x0, φ, t_engine)
        dF_2!(T, model, x0, φ, t_engine)

        mul!(T, -1.0)

        ldiv!(J, T)
        # T = J \ T

        # T.M_ij .*= -1.0
        # T.M_ij .= J.data .\ T.M_ij 
        
        r0.data .= J.data .\ r0.data
        # return (dx, T, r0)
        neumann!(dx, T, r0, memn; K=K)

        x0.data .= x1.data .- dx.data
        x1.data .= x0.data

        # return x0

        # end


    end

    return (;solution=x0, message="No Convergence") # The only allocation when workspace is preallocated

end


function newton_clinic(model::ADModel, workspace, tmr;
    K=20, tol_ε=1e-8, tol_η=1e-6, verbose=false, improve=false, interp_mode=:cubic,
    engine=:none
    )

    # mem = typeof(workspace) <: Nothing ? time_iteration_workspace(model) : workspace

    (;x0, x1, x2, dx, r0, J, φ, T, memn) = workspace

    if engine in (:cpu, :gpu)
        t_engine = get_backend(workspace.x0)
    else
        t_engine = nothing
    end


    for t=1:K
        
        Dolo.fit!(φ, x0)

        F!(r0, model, x0, φ, t_engine)

        ε = norm(r0)
        verbose ? println("$t: $ε") : nothing

        if ε<tol_ε
            return (;message="Solution found", solution=x0, n_iterations=t-1, dr=φ)
        end

        x1.data .= x0.data


        dF_1!(J, model, x0, φ, t_engine)
        dF_2!(T, model, x0, φ, t_engine)

        mul!(T, -1.0)

        ldiv!(J, T)
        # T = J \ T

        # T.M_ij .*= -1.0
        # T.M_ij .= J.data .\ T.M_ij 
        
        r0.data .= J.data .\ r0.data
        # return (dx, T, r0)
        neumann!(dx, T, r0, memn; K=K)

        x0.data .= x1.data .- dx.data
        x1.data .= x0.data

        # return x0

        # end


    end

    return (;solution=x0, message="No Convergence") # The only allocation when workspace is preallocated

end


# function time_iteration_1(model;
#     T=500,
#     K=10,
#     tol_ε=1e-8,
#     tol_η=1e-6,
#     verbose=false,
# )

#     N = length(model.grid)
#     x0 = GArray(model.grid, [SVector(model.calibration.x) for n=1:N])
#     x1 = deepcopy(x0)

#     local x0
#     local x1
#     local t


#     for t=1:T

#         r0 = F(model, x0, x0)

#         ε = norm(r0)

#         if ε<tol_ε
#             return (;message="Solution found", solution=x0, n_iterations=t)
#         end
#         # println("Iteration $t: $ε")
#         if verbose
#             println("ϵ=$(ε)")
#         end

#         x1.data .= x0.data
        
#         for k=1:K

#             r1 = F(model, x1, x0)
#             J = dF_1(model, x1, x0)
#             # dF!(J, model, x1, x0)

#             dx = GArray(model.grid, J.data .\ r1.data)

#             η = norm(dx)
            
#             x1.data .-= dx.data

#             verbose ? println(" - $k: $η") : nothing

#             if η<tol_η
#                 break
#             end

#         end

#         x0 = x1

#     end

#     return (;solution=x0, message="No Convergence", n_iterations=t)

# end

# using NLsolve

# function time_iteration_2(model;
#     T=500,
#     K=10,
#     tol_ε=1e-8,
#     tol_η=1e-6,
#     verbose=false
# )

#     N = length(model.grid)
#     x0 = GArray(model.grid, [SVector(model.calibration.x) for n=1:N])
#     x1 = deepcopy(x0)

#     local x0
#     local x1
#     local t


#     for t=1:T

#         function fun(u::AbstractVector{Float64})
#             x = unravel(x0, u)
#             r = F(model, x, x0)
#             return ravel(r)
#         end

#         function dfun(u::AbstractVector{Float64})
#             x = unravel(x0, u)
#             dr = dF_1(model, x, x0)
#             J = convert(Matrix,dr)
#             return J
#         end

#         u0 = ravel(x0)
#         sol = nlsolve(fun, dfun, u0)
#         u1 = sol.zero

#         η = maximum(u->abs(u[1]-u[2]), zip(u0,u1))
        
#         verbose ? println("Iteration $t: $η") : nothing
        
#         x0 = unravel(x0, u1)

#         if η<tol_η
#             return (;solution=x0, message="Convergence", n_iterations=t)
#         end

#     end

#     return (;message="No Convergence", n_iterations=T)

# end

# using LinearAlgebra

# function time_iteration_3(model;
#     T=500,
#     K=10,
#     tol_ε=1e-8,
#     tol_η=1e-8,
#     verbose=false,
#     improve=true,
#     x0=nothing,
#     interp_mode=:linear
# )

#     N = length(model.grid)
#     if x0===nothing
#         x0 = initial_guess(model)
#     else
#         x0 = deepcopy(x0)
#     end
#     # x0 = GArray(model.grid, [SVector(model.calibration.x) for n=1:N])

#     x1 = deepcopy(x0)

#     # local x0
#     local x1
#     local t


#     for t=1:T

#         φ = DFun(model, x0; interp_mode=interp_mode)

#         r0 = F(model, x0, φ)
#         ε = norm(r0)

#         if ε<tol_ε
#             return (;message="Solution found", solution=x0, n_iterations=t)
#         end

#         # φ = x0
#         x1.data .= x0.data
        
#         for k=1:K

#             r1 = F(model, x1, φ)
#             J = dF_1(model, x1, φ)
#             # dF!(J, model, x1, x0)

#             dx = GArray(model.grid, J.data .\ r1.data)

#             η = norm(dx)
            
#             x1.data .-= dx.data

#             verbose ? println(" - $k: $η") : nothing

#             if η<tol_η
#                 break
#             end

#         end

#         # ε = norm(F(model, x1, φ))
#         η = distance(x1,x0)

#         if !(improve)
#             x0.data[:] = x1.data[:]
#         else
#             # x = T(x)
#             # xnn = T(xn)
#             # x - xnn = -T'(x) (x - xn)
#             # x = xnn - T' (x - xn)
#             # x = (I-T')\(xnn - T' xn)

#             # # this version assumes same number of shocks

#             J_1 = Dolo.dF_1(model, x1, φ)
#             J_2 =  Dolo.dF_2(model, x1, x0)
#             J_2.M_ij[:] *= -1.0
#             Tp = J_1 \ J_2
#             r = x1 - Tp * x0
#             x0 = invert(r, Tp; K=500)

#         end
        
#         ε = maximum(abs, ravel(F(model, x0, x0)))

#         verbose ? println("Iteration $t : $η : $ε") : nothing

#         if η<tol_η
#             return (;solution=x0, dr=φ, message="Convergence", n_iterations=t)
#         end

#     end

#     return (;solution=x0, message="No Convergence", n_iterations=T)

# end



# # function time_iteration(model;
# #     T=500,
# #     K=10,
# #     tol_ε=1e-8,
# #     tol_η=1e-6,
# #     verbose=false
# # )


# #     N = length(model.grid)

# #     x0 = GArray(model.grid, [SVector(model.x) for n=1:N])
# #     x1 = deepcopy(x0)
# #     dx = deepcopy(x0)

# #     r0 = x0*0
    
# #     J = dF0(model, x0, x0)[2]

# #     local x0
# #     local x1

# #     for t=1:T
# #         # r0 = F(model, x0, x0)
# #         F!(r0, model, x0, x0)
# #         ε = norm(r0)
# #         if ε<tol_ε
# #             break
# #         end
# #         if verbose
# #             println("ϵ=$(ε)")
# #         end
# #         x1.data .= x0.data
# #         for k=1:K
# #             # r = F(model, x1, x0)
# #             F!(r0, model, x1, x0)
# #             # J = dF(model, x1, x0)
# #             dF!(J, model, x1, x0)
# #             # dx = J\r0
# #             for n=1:length(r0)
# #                 dx.data[n] = J.data[n]\r0.data[n]
# #             end
# #             e = norm(dx)
# #             # println("e=$(e)")
# #             x1.data .-= dx.data
# #             if e<tol_η
# #                 break
# #             end
# #         end
# #         x0 = x1

# #     end
# #     return x0
# # end



# # # Alternative implementations

# # function F0(model, s, x::SVector, xfut::GArray)
# #     tot = SVector((x*0)...)
# #     for (w, S) in τ(model, s, x)
# #         ind = (S[1], S[3])
# #         X = xfut(ind...)
# #         tot += w*arbitrage(model,s,x,S,X)
# #     end
# #     return tot
# # end


# # F(model, controls::GArray, φ::GArray) =
# #     GArray(
# #         model.grid,
# #         [F(model,s,x,φ) for (s,x) in zip(iti(model.grid), controls) ],
# #     )

# # function F!(out, model, controls, φ) 
# #     # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
# #     n=0
# #     for s in iti(model.grid)
# #         n += 1
# #         x = controls.data[n]
# #         out.data[n] = F(model,s,x,φ)
# #     end
# #     # end
# # end

# # dF(model, controls::GArray, φ::GArray) =
# #     GArray(    # this shouldn't be needed
# #         model.grid,
# #         [
# #             ForwardDiff.jacobian(u->F(model, s, u, φ), x)
# #             for (s,x) in zip(iti(model.grid), controls) 
# #         ]
# #     )

# # function dF!(out, model, controls, φ) 
# #     # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
# #     n=0
# #     for s in iti(model.grid)
# #         n += 1
# #         x = controls.data[n]
# #         out.data[n] = ForwardDiff.jacobian(u->F(model, s, u, φ), x)
# #     end
# #     # end
# # end

# # # function dF2!(out, model, controls, φ) 
# # #     # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
# # #     n=0
# # #     for s in iti(model.grid)
# # #         n += 1
# # #         x = controls.data[n]
# # #         out.data[n] = ForwardDiff.jacobian(u->F(model, s, x, φ), φ)
# # #     end
# # #     # end
# # # end


# # FdF(model, controls::GArray, φ::GArray) =
# #     GArray(
# #         model.grid,
# #         [
# #             (F(model,s,x,φ), ForwardDiff.jacobian(u->F(model, s, u, φ), x))
# #             for (s,x) in zip(iti(model.grid), controls) 
# #         ]
# #     )


# # function F0(model, controls::GArray, xfut::GArray)

# #     N = length(controls)
# #     res = GArray(
# #         model.grid,
# #         zeros(typeof(controls[1]), N)
# #     )
# #     for (i,(s,x)) in enumerate(zip(iti(model.grid), controls))
# #         res[i] = F(model,s,x,xfut)
# #     end
# #     return res
# # end



# # function dF0(model, controls::GArray, xfut::GArray)

# #     N = length(controls)
# #     res = deepcopy(controls)
# #     dres = GArray(
# #         model.grid,
# #         zeros(typeof(res[1]*res[1]'), N)
# #     )
# #     for (i,(s,x)) in enumerate(zip(iti(model.grid), controls))
# #         res[i] = F(model,s,x,xfut)
# #         dres[i] = ForwardDiff.jacobian(u->F(model, s, u, xfut), x)
# #     end
# #     return res, dres
# # end





