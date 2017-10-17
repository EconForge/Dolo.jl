function absmax(s::ListOfPoints)
    t = 0.0
    for p in s
        t = max(t, maximum(p))
    end
    t
end
absmax(x::Vector{<:ListOfPoints}) = maximum(absmax.(x))

function euler_residuals(model, s::ListOfPoints, x::Vector{<:ListOfPoints}, dr, dprocess, parms::SVector; with_jres=false, set_dr=true) #, jres=nothing, S_ij=nothing)

    if set_dr ==true
      set_values!(dr,x)
    end

    N_s = length(s) # Number of gris points for endo_var
    n_s = length(s[1]) # Number of states
    n_x = length(x[1][1]) # Number of controls

    n_ms = length(x)  # number of exo states today
    n_mst = n_inodes(dprocess,1)  # number of exo states tomorrow
    d = length(s[1])

    # res = zeros(n_ms, N_s, n_x)
    res = deepcopy(x)

    if with_jres == true
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
            if with_jres==true
                rr, rr_XM = arbitrage(model,(Val(0),Val(6)),m,s,x[i_ms],M,S,X,parms)
                J_ij[i_ms,I_ms][:] = w*rr_XM
                S_ij[i_ms,I_ms][:] = S
            else
                rr = arbitrage(model, m,s,x[i_ms],M,S,X,parms)
            end
            res[i_ms][:] += w*rr
        end
    end
    if with_jres
        return res,J_ij,S_ij
    else
        return res
    end
end

#
#
# function euler_residuals_2(model, s::ListOfPoints, x::Vector{<:ListOfPoints}, dr, dprocess, parms::SVector) #; with_jres=false, set_dr=true) #, jres=nothing, S_ij=nothing)
#     #
#     # if set_dr ==true
#     #   set_values!(dr,x)
#     # end
#
#     N_s = length(s) # Number of gris points for endo_var
#     n_s = length(s[1]) # Number of states
#     n_x = length(x[1][1]) # Number of controls
#
#     n_ms = length(x)  # number of exo states today
#     n_mst = n_inodes(dprocess,1)  # number of exo states tomorrow
#     d = length(s[1])
#
#     # res = zeros(n_ms, N_s, n_x)
#     R_i = deepcopy(x)
#     D_i = [Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}(N_s) for i=1:n_ms]
#
#     for i_ms in 1:n_ms
#         m = SVector(node(dprocess,i_ms)...)
#         for I_ms in 1:n_mst
#             M = SVector(inode(dprocess, i_ms, I_ms)...)
#             w = iweight(dprocess, i_ms, I_ms)
#             S, S_x = transition(model, (Val(0), Val(3)), m, s, x[i_ms], M, parms)
#             X = dr(i_ms, I_ms, S)
#             rr = arbitrage(model, m,s,x[i_ms],M,S,X,parms)
#             R_i[i_ms][:] += w*rr
#         end
#     end
#
#     return R_i, D_i
#
# end

#############################################################################
# I am still not sure we need it
function euler_residuals(model, x::Array{Float64,2}, dr,
                         dprocess, parms::AbstractArray; with_jres=false, set_dr=true,
                         jres=nothing, S_ij=nothing)
   N_m = n_nodes(dprocess.grid)
   x_reshaped = destack0(x,N_m)
   return euler_residuals(model,s,x_reshaped,dr,dprocess,parms; with_jres=false, set_dr=true,
                          jres=nothing, S_ij=nothing)
end



######################
function SerialDifferentiableFunction(f, epsilon=1e-6)

    function df(x::AbstractVector)

      v0 = f(x)

      n_m = size(v0,1)
      N_s = size(v0,2)
      n_v = size(v0,3)
      assert(size(x[1],1) == N_s)
      n_x = size(x[1],2)

      dv = zeros(n_m*N_s,  n_v, n_x)
      for i in 1:n_x
        xi = deepcopy(cat(1,x...))
        xi[:,i] += epsilon
        # You could also use f(xi)
        vi = f(destack0(xi,n_m))
        dd=(vi+(-1*v0))./epsilon
        # dd=permutedims(dd, [2,1,3])
        dv[:,:, i] = reshape(dd,N_s*n_m,n_x) # (1dim) corresponds to equations, in raws you first stuck derivatives wrt 1rst exo state, 2nd, etc
      end
      # dv_AA = AxisArray(dv, Axis{:N}(1:n_m*N_s), Axis{:n_v}(1:n_v), Axis{:n_x}(1:n_x))
      return [v0, dv]
    end


    function df(x::Array{Float64,2})

      v0 = f(x)

      # n_m = size(v0,1)
      N_s = size(v0,1)
      n_v = size(v0,2)
      assert(size(x,1) == N_s)
      n_x = size(x,2)

      dv = zeros(N_s,  n_v, n_x)
      for i in 1:n_x
         xi = deepcopy(x)
         xi[:,i] += epsilon
         vi = f(xi)
         dd=(vi+(-1*v0))./epsilon
         dv[:,:, i] = dd
      end
      # dv_AA = AxisArray(dv, Axis{:N}(1:N_s), Axis{:n_v}(1:n_v), Axis{:n_x}(1:n_x))
      return [v0, dv]

    end

    function df(x::Array{Float64,1})

      v0 = f(x)

      n_v = size(v0,1)
      # assert(size(x,1) == N_s)
      n_x = size(x,1)

      dv = zeros(n_v,n_x)
      for i in 1:n_x
        xi = deepcopy(x)
        xi[i,:]+= epsilon
        vi = f(xi)
        dd=(vi+(-1*v0))./epsilon
        dv[:,i] = dd
      end
      # dv_AA = AxisArray(dv, Axis{:N}(1:N_s), Axis{:n_v}(1:n_v), Axis{:n_x}(1:n_x))
      return [v0, dv]

    end
    df
end


function ssmul(A::Array{Float64,3},B::Array{Float64,2})
  # simple serial_mult (matrix times vector)
  N,a,b = size(A)
  NN,b = size(B)
  O = zeros(N,a)
  for n in 1:N
    for k in 1:a
      for l in 1:b
        O[n,k] += A[n,k,l]*B[n,l]
      end
    end
  end
  return O

end

function destack0(x::Array{Float64,3},n_m::Int)
  xx=collect(x)
  return [xx[i, :, :] for i=1:n_m]
end

####


function reorder_data(dprocess, res, dres, jres, fut_S)

    R_i = [to_LOP(res[i,:,:]) for i=1:size(res,1)]
    D_i =  [to_LOJ(dres[i,:,:,:]) for i=1:size(res,1)]
    J_ij = typeof(D_i[1])[to_LOJ(jres[i,j,:,:,:]) for i=1:size(jres,1), j=1:size(jres,2)]

    d = size(fut_S,4)
    S_ij = Vector{Point{d}}[to_LOP(fut_S[i,j,:,:]) for i=1:size(fut_S,1), j=1:size(fut_S,2)]

    Π_i = deepcopy(R_i)
    π_i = deepcopy(R_i)

    return R_i, D_i, J_ij, S_ij

end

function invert!(A)
    # A[i] <- (A[i])^(-1)
    N = length(A)
    for n=1:N
        A[n] = inv(A[n])
    end
end

function premult!(A,B)
    # B[i] <- A[i]*B[i]
    N = length(A)
    for n=1:N
        B[n] = A[n]*B[n]
    end
end

function addmul!(O,A,B)
    # O[i] <- A[i]*B[i]
    N = length(A)
    for n=1:N
        O[n] += A[n]*B[n]
    end
end


import Base
function Base.maxabs(v::ListOfPoints{d}) where d
    m = 0.0
    for n = 1:length(v)
        mm = maximum(abs, v[n])
        if mm>m
            m = mm
        end
    end
    m
end

function Base.maxabs(v::Vector{ListOfPoints{d}}) where d
    maximum(maxabs.(v))
end


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


function invert_jac(res::Array{Float64,3},dres::Array{Float64,4},jres::Array{Float64,5},
            fut_S::Array{Float64,4}, dumdr; tol::Float64=1e-10,
            maxit::Int=1000, verbose::Bool=false)

    dprocess = dumdr.process
    R_i, D_i, J_ij, S_ij = reorder_data(dprocess, res, dres, jres, fut_S)
    sol, it, lam, errors =  invert_jac(R_i, D_i, J_ij, S_ij, dumdr, tol=tol, maxit=maxit, verbose=verbose)
    n_ms = length(sol)
    N = length(sol[1])
    n_x = length(sol[1][1])
    rsol = cat(1, [reshape(from_LOP(tt),1,N,n_x) for tt in sol]...)
    return rsol, it, lam, errors
end

function invert_jac(R_i, D_i, J_ij, S_ij, dumdr; tol::Float64=1e-10,
            maxit::Int=1000, verbose::Bool=false)

    dprocess = dumdr.process

    Dinv = deepcopy(D_i)
    for i=1:length(Dinv)
        invert!(Dinv[i])
    end

    M_ij = deepcopy(J_ij)
    for i=1:size(M_ij,1)
        for j=1:size(M_ij,2)
            premult!(Dinv[i],M_ij[i,j])
            # M_ij[i,j][:] *= Dolo.iweight(dprocess,i,j)
        end
    end

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
    for i=1:length(Dinv)
        premult!(Dinv[i], π_i[i])
    end
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


###### radius_jac
function radius_jac(res::AbstractArray,dres::AbstractArray,jres::AbstractArray,
  S_ij::AbstractArray,
  dumdr; tol=1e-08, maxit=1000, verbose=false, precomputed = false)

  n_m, N_s, n_x = size(res)
  err0 = 0.0
  ddx = rand(size(res))*10^9

  lam = 0.0
  lam_max = 0.0

  lambdas = zeros(maxit)
  if verbose==true
    print("Starting inversion. Radius_Jac")
  end

  for nn in 1:maxit
    ddx /= maximum(abs, ddx)
    d_filt_dx(ddx,jres,S_ij,dumdr; precomputed=precomputed)
    lam = maximum(abs, ddx)
    lam_max = max(lam_max, lam)
    lambdas[nn] = lam
  end

  return (lam, lam_max, lambdas)

end


function inplace(Phi::Array{Float64,2}, J::Array{Float64,4})
  a,b,c,d = size(J)

  for i_a in 1:a
    for i_b in 1:b
      for i_c in 1:c
        [J[i_a,i_b,i_c,i_d] *= Phi[i_b,i_c] for i_d in 1:d]
      end
    end
  end
  return J
end

function smooth_right(res::Array{Float64,2}, dres::Array{Float64,3}, jres::Array{Float64,4},
                          dx::Array{Float64,2}; pos::Float64=1.0)

   BIG = 1e20
   dxinf = -dx .<= -BIG   #  isinf(a) |

   n_x = size(res,2)
   sq = sqrt(res.^2 + dx.^2 )
   H = res + dx - sq
   H[dxinf] = res[dxinf]
   Phi_a = 1 - res./sq
   Phi_b = 1 - dx./sq

   Phi_a[dxinf] = 1.0
   Phi_b[dxinf] = 0.0

   H_x = reshape(Phi_a,(size(Phi_a,1),size(Phi_a,2),1)).*dres

   for i_x in 1:n_x
       H_x[:,i_x,i_x] += Phi_b[:,i_x].*pos
   end

   jres = inplace(Phi_a, jres)

   return H, H_x, jres

end


function smooth_right(res::Array{Float64,2},  dx::Array{Float64,2})
    BIG = 1e20
    dxinf = -dx .<= -BIG
    n_x = size(res,2)
    sq = sqrt(res.^2 + dx.^2 )
    H = res + dx - sq
    H[dxinf] = res[dxinf]
    return H
end
