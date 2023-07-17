function τ(model::Dolo.YModel, s0::QP, a::SVector)
    (transition(model, s0, a))
end


### transition function

function τ(dmodel::Dolo.DYModel{M}, ss::T, a::SVector) where M<:Union{Dolo.YModel{<:Dolo.VAR1},Dolo.YModel{<:Dolo.MarkovChain}}  where T<:QP


    (i,_) = ss.loc
    s_ = ss.val
    

    # TODO: replace following block by one nonallocating function
    Q = dmodel.dproc.Q
    P = dmodel.dproc.P

    it = (
        (
            P[i,j],
            let 
                S_exo = Q[j]
                S_endo = SVector(transition(dmodel.model, s_, a, Q[j])...)
                S = SVector(S_exo..., S_endo...)
                QP((j,S_endo),S)
                # (loc=(j,S_endo),val=S)
            end
        )
        for j in 1:size(P, 2)
    )

    it

end


function τ(dmodel::Dolo.DYModel{M}, ss::T, a::SVector)  where M<:Dolo.YModel{<:Dolo.MvNormal} where T<:QP


    # ind = ss.loc
    s_ = ss.val
    

    # TODO: replace following block by one nonallocating function

    # m,s = Dolo.split_states(dmodel.model, s_)

    x = dmodel.dproc.x
    w = dmodel.dproc.w

    it = (
        (
            w[j],
            let 
                # S_exo = Q[j]
                S = transition(dmodel.model, s_, a, x[j])
                QP(S,S)
                # (loc=S,val=S)
            end
        )
        for j in 1:length(w)
    )

    it

end


function trembling__hand(g::CGrid{1}, xv)
    x = xv[1]
    r = g.ranges[1]
    u = (x-r[1])/(r[2]-r[1])
    n = r[3]

    i_ = floor(Int, u*(n-1))
    i_ = max(min(i_,n-2),0)
    λ = u*(n-1)-i_

    λ = min(max(λ, 0.0), 1.0)
    
    (
        (1-λ, i_+1),
        (λ, i_+2)
    )
end

using ResumableFunctions

@resumable function τ_fit(model, ss::Tuple, a::SVector; linear_index=false)

    p = model.calibration.p
    P = model.transition
    Q = model.grid.g1.points

    i = ss[1][1]

    n_m = length(model.calibration.m)

    (i,_),(s_) = ss

    # TODO: replace following block by one nonallocating function
    k  = length(model.calibration.m)
    l = length(model.calibration.s)

    m = SVector((s_[__i] for __i=1:k)...)
    s = SVector((s_[__i] for __i=k+1:(k+l))...)

    for j in 1:size(P, 2)

        M = Q[j]
        S = transition(model, m, s, a, M, p)

        for (w, i_S) in trembling__hand(model.grid.g2, S)

            res = (
                P[i,j]*w,

                (
                    (linear_index ? to__linear_index(model.grid, (j,i_S)) : (j,i_S)),

                    SVector(M..., model.grid.g2[i_S]...)
                )
            )
            if linear_index
                res
                # res::Tuple{Float64, Tuple{Int64, Tuple{SVector{1},SVector{1}}}}
            else
                res
                # res::Tuple{Float64, Tuple{Tuple{Int64, Int64}, Tuple{SVector{1},SVector{1}}}}
            end
            @yield res
        end
    end
end

function G(model, μ::GDist{T}, x) where T
    μ1 = GArray(μ.grid, zeros(Float64, length(μ)))
    for ss in iti(model.grid)
        a = x[ss[1]...]
        for (w, (ind, _)) in τ_fit(model, ss, a)
            μ1[ind...] += w*μ[ind...]
        end
    end

    μ1

end

# TODO: write interpolation version of G


function transition_matrix(model, x)

    N = length(x)
    P = zeros(N,N)

    for (ss,a) in zip(enum(model.grid),x)
        ind_i = ss[1]
        i = to__linear_index(model.grid, ind_i)
        for (w, (ind_j, _)) in τ_fit(model, ss, a)
            j = to__linear_index(model.grid, ind_j)
            P[i,j] = w
        end
    end

    P

end


using LinearAlgebra: I

function ergodic_distribution(model, x::GArray)

    P = transition_matrix(model, x)

    PP = [ (P-I)[:,1:end-1] ;;  ones(size(P,1))]
    R = zeros(size(PP,1))
    R[end] = 1
    μ = PP'\R
    
    ergo = GDist(model.grid, μ)

    return ergo

end


## TODO: some ways to plot the ergo dist...

τ(model, ss::Tuple, φ::DFun) = τ(model, ss, φ(ss))




## TODO: some simulation

# TODO
# function simulate(model, ss::Tuple, φ; T=20)
#     Typ = typeof(ss)
#     sim = Typ[ss]
#     for t = 1:T
#     end
# end