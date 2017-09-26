


function euler_residuals(model, s::AbstractArray, x::Array{Array{Float64,2},1}, dr,
                         dprocess, parms::AbstractArray; with_jres=false, set_dr=true,
                         jres=nothing, S_ij=nothing)

    if set_dr ==true
      set_values!(dr,x)
    end

    # s = model.grid.nodes
    N_s = size(s,1) # Number of gris points for endo_var
    n_s = size(s,2) # Number of states
    n_x = size(x[1],2) # Number of controls

    # P = dprocess.values
    # Q = dprocess.transitions

    # n_ms = Dolo.n_nodes(dprocess)  # number of exo states today
    n_ms = size(x,1)  # number of exo states today
    n_mst = n_inodes(dprocess,1)  # number of exo states tomorrow
    # n_mv = size(Dolo.node(dprocess, 1),1)  # number of exo variable

    res = zeros(n_ms, N_s, n_x)

    if with_jres == true
      if jres== nothing
        jres = zeros((n_ms,n_mst,N_s,n_x,n_x))
      end
      if S_ij== nothing
        S_ij = zeros((n_ms,n_mst,N_s,n_s))
      end
    end


    for i_ms in 1:n_ms
       m = node(dprocess,i_ms)
       for (w, M, I_ms) in get_integration_nodes(dprocess,i_ms)
           # Update the states
           S = transition(model, m, s, x[i_ms], M, parms)
           X = dr(i_ms, I_ms, S)
           if with_jres==true
               ff = SerialDifferentiableFunction(u->arbitrage(model, m,s,x[i_ms],M,S,u,parms))
               rr, rr_XM = ff(X)
               jres[i_ms,I_ms,:,:,:] = w*rr_XM
               S_ij[i_ms,I_ms,:,:] = S
           else
               rr = arbitrage(model, m,s,x[i_ms],M,S,X,parms)
           end
           res[i_ms,:,:] += w*rr
        end
    end


    # res_AA = AxisArray(res, Axis{:n_m}(1:n_ms), Axis{:N_s}(1:N_s), Axis{:n_x}(1:n_x))
    if with_jres==true
        return res, jres, S_ij
    else
        return res
    end

end
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
function SerialDifferentiableFunction(f, epsilon=1e-8)

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

function d_filt_dx(res::Array{Float64,3},jres::Array{Float64,5},S_ij::Array{Float64,4},
  dumdr; precomputed::Bool=false)
  n_m=size(jres,1)
  n_mt=size(jres,2)
  set_values!(dumdr,destack0(res, n_m))
  for i in 1:n_m
    res[i,:,:] = 0
    for j in 1:n_mt
      A = jres[i,j,:,:,:]
      if precomputed== false
        # B = dumdr(j,S_ij[i,j,:,:])
        B = dumdr(i,j,S_ij[i,j,:,:])
      else
        B = dumdr(i,j) #,fut_S[i,j,:,:])
      end
      res[i,:,:] += ssmul(A,B)
    end
  end
  return res
end


function invert_jac(res::Array{Float64,3},dres::Array{Float64,4},jres::Array{Float64,5},
  fut_S::Array{Float64,4}, dumdr; tol::Float64=1e-10,
  maxit::Int=1000, verbose::Bool=false)
  n_m, N_s, n_x = size(res)
  ddx = zeros(n_m,N_s,n_x)
  # A=deepcopy(dres)
  # B = deepcopy(res)
  # for i_m in 1:n_m
  #     for n in 1:N_s
  #        ddx[i_m,n,:]= invert(A[i_m,n,:,:],collect(B[i_m,n,:]))[2]
  #     end
  # end
  A=copy(dres)
  B = copy(res)
  for i_m in 1:n_m
    for n in 1:N_s
      ddx[i_m,n,:]= A[i_m,n,:,:]\B[i_m,n,:]
    end
  end


  # if filt == nothing
  #   error("No filter supplied.")
  # else
  #   dumdr = filt
  # end

  lam = -1.0
  lam_max = -1.0
  err_0 = maximum(abs, ddx)
  tot = deepcopy(ddx)

  if verbose==true
    print("Starting inversion")
  end
  err=err_0


  verbose && println(repeat("-", 35))
  verbose && @printf "%-6s%-12s%-5s\n" "err" "gain" "gain_max"
  verbose && println(repeat("-", 35))
  precomputed=false

  it = 0
  while it<maxit && err>tol
    it +=1

    ddx = d_filt_dx(ddx,jres,fut_S,dumdr; precomputed=precomputed)
    # might also work
    # d_filt_dx(ddx,jres,fut_S,n_m,N,n_x,dumdr; precomputed=precomputed)

    err = maximum(abs, ddx)
    lam = err/err_0
    lam_max = max(lam_max, lam)
    tot += ddx
    err_0 = err
    verbose && @printf "%-6f%-12.2e%-5.2e\n" err lam lam_max
  end
  if err<tol
    ddx = d_filt_dx(ddx,jres,fut_S,dumdr; precomputed=precomputed)
    err = maximum(abs, ddx)

    lam = err/err_0
    lam_max = max(lam_max, lam)
  end

  tot += ddx*lam/(1-lam)
  return tot, it, lam
end

function invert_jac(res::AbstractArray,dres::AbstractArray,jres::Array{Float64,5},
  fut_S::Array{Float64,4}; tol::Float64=1e-10,
  maxit::Int=1000, verbose::Bool=false)
  return error("No filter supplied.")
end


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
