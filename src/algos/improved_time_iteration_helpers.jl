

# function euler_residuals_noalloc(model, s::ListOfPoints, x::MSM{Point{n_x}}, dr, dprocess, parms::SVector; keep_J_S=false, out=nothing) where n_x #, jres=nothing, S_ij=nothing)
#     #
#     # if set_dr ==true
#     #   set_values!(dr,x)
#     # end
#     if out === nothing
#         res = zeros_like(x)::MSM{Point{n_x}}
#     else
#         res = out
#     end

#     res.data .*= 0.0

#     N_s = length(s) # Number of gris points for endo_var
#     n_s = length(s[1]) # Number of states

#     n_ms = length(x.views)  # number of exo states today
#     n_mst = n_inodes(dprocess,1)  # number of exo states tomorrow
#     d = length(s[1])

#     # TODO: allocate properly...
    
#     # if keep_J_S
#     #     jres = zeros((n_ms,n_mst,N_s,n_x,n_x))
#     #     fut_S = zeros((n_ms,n_mst,N_s,n_s))
#     #     J_ij = Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}[to_LOJ(jres[i,j,:,:,:]) for i=1:size(jres,1), j=1:size(jres,2)]
#     #     S_ij = Vector{Point{d}}[to_LOP(fut_S[i,j,:,:]) for i=1:size(fut_S,1), j=1:size(fut_S,2)]
#     # end


#     # xx = x.views[1]
#     for i_ms in 1:n_ms
#         m = node(Point, dprocess,i_ms)
#         for I_ms in 1:n_mst
#             M = inode(Point, dprocess, i_ms, I_ms)
#             w = iweight(dprocess, i_ms, I_ms)
#             for n=1:N_s
#                 s_n = s[n]
#                 x_n = dr(i_ms, s_n)
#                 S_n = transition(model, m, s_n, x_n, M, parms)
#                 X_n = dr(i_ms, I_ms, S_n)
#                 r_n = arbitrage(model, m, s_n, x_n, M,S_n,X_n, parms)
#                 res.views[i_ms][n] += w*r_n
#             end
#         end
#     end


#     if keep_J_S
#         return res,J_ij,S_ij
#     else
#         return res
#     end
# end

####################################
# one in place filtering step: M.r #
####################################


function d_filt_dx!(Π_i::MSM{Point{n_x}}, π_i::MSM{Point{n_x}}, M_ij, S_ij::Matrix{ListOfPoints{d}}, dumdr::Dolo.CachedDecisionRule) where n_x where d
    n_m,n_mt = size(M_ij)
    Dolo.set_values!(dumdr,π_i)
    for i in 1:n_m
        Π_i.views[i][:] .*= 0.0
        for j in 1:n_mt
            A = M_ij[i,j]
            B = dumdr(i,j,S_ij[i,j])
            addmul!(Π_i.views[i],A,B)
        end
    end
end

function d_filt_dx!(Π_i::Vector{ListOfPoints{n_x}}, π_i::AbstractVector{<:AbstractVector{Point{n_x}}}, M_ij, S_ij::Matrix{ListOfPoints{d}}, dumdr::Dolo.CachedDecisionRule) where n_x where d
    n_m,n_mt = size(M_ij)
    Dolo.set_values!(dumdr,π_i)
    for i in 1:n_m
        Π_i[i][:] *= 0
        for j in 1:n_mt
            A = M_ij[i,j]
            B = dumdr(i,j,S_ij[i,j])
            addmul!(Π_i[i],A,B)
        end
    end
end



function preinvert!(π_i, D_i, J_ij, S_ij)

    for i=1:length(D_i)
        invert!(D_i[i])
    end
    for i=1:size(J_ij,1)
        for j=1:size(J_ij,2)
            premult!(D_i[i],J_ij[i,j])
        end
    end
    for i=1:length(D_i)
        premult!(D_i[i], π_i[i])
    end

    return π_i, J_ij, S_ij

end

function invert_jac(R_i, M_ij, S_ij, dumdr; tol::Float64=1e-10,
            maxit::Int=1000, verbose::Bool=false)

    dprocess = dumdr.process

    lam = -1.0
    lam_max = -1.0

    if verbose==true
    print("Starting inversion")
    end

    verbose && println(repeat("-", 35))
    verbose && @printf "%-6s%-12s%-5s\n" "err" "gain" "gain_max"
    verbose && println(repeat("-", 35))
    precomputed=false

    π_i = deepcopy(R_i)
    Π_i = deepcopy(π_i)
    tot = deepcopy(π_i) # should premult by D_i

    err_0 = maxabs(π_i)
    err = 10000

    errors = [err_0]

    it = 0
    while it<maxit && err>tol
        it +=1
        d_filt_dx!(Π_i, π_i, M_ij, S_ij, dumdr)
        π_i = Π_i
        err = maxabs(Π_i)
        push!(errors, err)
        lam = err/err_0
        lam_max = max(lam_max, lam)
        err_0 = err
        verbose && @printf "%-6f%-12.2e%-5.2e\n" err lam lam_max
        for i=1:length(tot)
            tot[i] += Π_i[i]
        end
    end

    return tot, it, lam, errors

end

####
# The LinearThing represents
# - (I-M) with (`L*...``)
# - M with (`mul(L,...)``)
####

mutable struct LinearThing
   M_ij
   S_ij
   I
   counter::Int
end

import LinearAlgebra
import LinearAlgebra: *, mul!
using LinearAlgebra

LinearThing(M,S,I) = LinearThing(M,S,I,0)


eltype(L::LinearThing) = Float64
function shape(L::LinearThing)
   n_m = size(L.M_ij,1)
   N = size(L.M_ij[1,1],1)
   n_x = size(L.M_ij[1,1][1],1)
   return (n_x, N, n_m)
end
size(L::LinearThing,d) = prod(shape(L))


function *(L::LinearThing,x::MSM{Point{n_x}}) where n_x
    xx = x*1.0
    Dolo.d_filt_dx!(xx, x, L.M_ij, L.S_ij, L.I)
    L.counter += 1
    return xx
 end
 


# function *(L::LinearThing,x::AbstractVector{<:AbstractVector{Point{n_x}}}) where n_x
#    xx = [copy(e) for e in x]
#    Dolo.d_filt_dx!(xx, x, L.M_ij, L.S_ij, L.I)
#    L.counter += 1
#    return xx
# end

# function *(L::LinearThing,m::AbstractArray{Float64, 3})
#    n_x,N,n_m = size(m)
#    # TODO remove copy there
#    x = [copy(reshape(reinterpret(SVector{n_x,Float64},vec(m[:,:,i])),(N,))) for i=1:n_m]
#    y = deepcopy(x)
#    xx = L*y
#    rr = [reshape(reinterpret(Float64, vec(xx[i])), (n_x,N)) for i=1:length(xx)]
#    rrr = cat(rr...; dims=3)
#    return reshape(rrr, n_x,N,n_m)
# end

# function *(L::LinearThing,v::AbstractVector{Float64})
#    m = copy(v)
#    sh = shape(L)
#    n_x = sh[1]
#    vv = reshape(reinterpret(Point{n_x},vec(m)),(sh[2],sh[3]))
#    # x = [view(vv,:,i) for i=1:sh[3]]
#    x = [vv[:,i] for i=1:sh[3]]
#    y = x-L*x
#    mm = reshape(m, sh...)
#    mmm = mm-L*mm
#    return mmm[:]
# end

# function mul(L::LinearThing,v::AbstractVector{Float64})
#    m = copy(v)
#    sh = shape(L)
#    mm = reshape(m, sh...)
#    mmm = L*mm
#    return mmm[:]
# end


# function A_mul_B!(w::AbstractVector{Float64},L::LinearThing,v::AbstractVector{Float64})
#    w[:] = L*v
#    return w
# end


# function mul!(w,L::LinearThing,v)
#    w[:] = L*v
#    return w
# end


####

mutable struct ImprovedTimeIterationResult
  dr::AbstractDecisionRule
  N::Int
  f_x::Float64
  d_x::Float64
  x_converged::Bool
  complementarities::Bool
  dprocess::AbstractDiscretizedProcess
  # Time_search::
  tol::Float64
  Lambda::Float64
  N_invert::Float64
  N_search::Float64
  radius::Float64
  trace_data
end


converged(r::ImprovedTimeIterationResult) = r.x_converged

function Base.show(io::IO, r::ImprovedTimeIterationResult)
  @printf io "Results of Improved Time Iteration Algorithm\n"
  @printf io " * Number of iterations: %s\n" string(r.N)
  @printf io " * Complementarities: %s\n" string(r.complementarities)
  @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
  @printf io " * Convergence: %s\n" converged(r)
  @printf io " * Contractivity: %s\n" string(r.Lambda)
  @printf io "   * |x - x'| < %.1e: %s\n" r.tol r.x_converged
end
