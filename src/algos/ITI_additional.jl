function euler_residuals(model, s::ListOfPoints, x::Vector{<:ListOfPoints}, dr, dprocess, parms::SVector; keep_J_S=false, set_dr=true) #, jres=nothing, S_ij=nothing)
    #
    if set_dr ==true
      set_values!(dr,x)
    end

    N_s = length(s) # Number of gris points for endo_var
    n_s = length(s[1]) # Number of states
    n_x = length(x[1][1]) # Number of controls

    n_ms = length(x)  # number of exo states today
    n_mst = n_inodes(dprocess,1)  # number of exo states tomorrow
    d = length(s[1])

    # TODO: allocate properly...
    res = deepcopy(x)
    for i_m=1:length(res)
        res[i_m][:] *= 0.0
    end

    if keep_J_S
        jres = zeros((n_ms,n_mst,N_s,n_x,n_x))
        fut_S = zeros((n_ms,n_mst,N_s,n_s))
        J_ij = Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}[to_LOJ(jres[i,j,:,:,:]) for i=1:size(jres,1), j=1:size(jres,2)]
        S_ij = Vector{Point{d}}[to_LOP(fut_S[i,j,:,:]) for i=1:size(fut_S,1), j=1:size(fut_S,2)]
    end

    for i_ms in 1:n_ms
        m = SVector(node(dprocess,i_ms)...)
        for I_ms in 1:n_mst
            M = SVector(inode(dprocess, i_ms, I_ms)...)
            w = iweight(dprocess, i_ms, I_ms)
            S = transition(model, m, s, x[i_ms], M, parms)
            X = dr(i_ms, I_ms, S)
            if keep_J_S==true
                rr, rr_XM = arbitrage(model,(Val(0),Val(6)),m,s,x[i_ms],M,S,X,parms)
                J_ij[i_ms,I_ms][:] = w*rr_XM
                S_ij[i_ms,I_ms][:] = S
            else
                rr = arbitrage(model, m,s,x[i_ms],M,S,X,parms)
            end
            res[i_ms][:] += w*rr
        end
    end
    if keep_J_S
        return res,J_ij,S_ij
    else
        return res
    end
end

####################################
# one in place filtering step: M.r #
####################################

function d_filt_dx!(Π_i::Vector{ListOfPoints{n_x}}, π_i::Vector{ListOfPoints{n_x}}, M_ij, S_ij::Matrix{ListOfPoints{d}}, dumdr::Dolo.CachedDecisionRule) where n_x where d
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

LinearThing(M,S,I) = LinearThing(M,S,I,0)

eltype(L::LinearThing) = Float64
function shape(L::LinearThing)
   n_m = size(L.M_ij,1)
   N = size(L.M_ij[1,1],1)
   n_x = size(L.M_ij[1,1][1],1)
   return (n_x, N, n_m)
end
size(L::LinearThing,d) = prod(shape(L))


function *(L::LinearThing,x::Vector{ListOfPoints{n_x}}) where n_x
   xx = deepcopy(x)
   Dolo.d_filt_dx!(xx, x, L.M_ij, L.S_ij, L.I)
   L.counter += 1
   return xx
end

function *(L::LinearThing,m::Array{Float64, 3})
   n_x,N,n_m = size(m)
   x = [reinterpret(SVector{n_x,Float64},m[:,:,i],(N,)) for i=1:n_m]
   y = deepcopy(x)
   xx = L*y
   rr = [reinterpret(Float64, xx[i], (n_x,N)) for i=1:length(xx)]
   rrr = cat(3,rr...)
   return reshape(rrr, n_x,N,n_m)
end

function *(L::LinearThing,v::AbstractVector{Float64})
   m = copy(v)
   sh = shape(L)
   n_x = sh[1]
   vv = reinterpret(Point{n_x},m,(sh[2],sh[3]))
   # x = [view(vv,:,i) for i=1:sh[3]]
   x = [vv[:,i] for i=1:sh[3]]
   y = x-L*x
   mm = reshape(m, sh...)
   mmm = mm-L*mm
   return mmm[:]
end

function mul(L::LinearThing,v::AbstractVector{Float64})
   m = copy(v)
   sh = shape(L)
   mm = reshape(m, sh...)
   mmm = L*mm
   return mmm[:]
end


function A_mul_B!(w::AbstractVector{Float64},L::LinearThing,v::AbstractVector{Float64})
   w[:] = L*v
   return w
end



####

type ImprovedTimeIterationResult
  dr::AbstractDecisionRule
  N::Int
  f_x::Float64
  d_x::Float64
  x_converged::Bool
  complementarities::Bool
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
