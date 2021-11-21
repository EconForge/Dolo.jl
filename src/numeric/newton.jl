function add_epsilon!(x::ListOfPoints{d}, i, epsilon) where d
  ei = SVector{d,Float64}([(j==i ? epsilon : 0.0) for j=1:d])
  for i=1:length(x)
    x[i] += ei
  end
end

function DiffFun(fun, x0::Vector{ListOfPoints{n_x}}, epsilon=1e-6) where n_x

    xi = deepcopy(x0)
    N = length(x0[1])
    n_m = length(x0)
    r0 = fun(x0)::Vector{ListOfPoints{n_x}}
    JMat = [zeros(n_x, n_x, N) for i=1:n_m]
    for i_x=1:n_x
      for i_m=1:n_m
        add_epsilon!(xi[i_m], i_x, epsilon)
      end
      fi = fun(xi)::Vector{ListOfPoints{n_x}}
      di = (fi-r0)/epsilon
      for i_m=1:n_m
        JMat[i_m][:,i_x,:] = reshape(reinterpret(Float64, vec(di[i_m])), (n_x, N))
        add_epsilon!(xi[i_m], i_x, -epsilon)
      end
    end
    J = [reshape(reinterpret(SMatrix{n_x,n_x,Float64,n_x^2},vec(JMat[i])),(N,)) for i=1:n_m]

    return (r0,J) #::Tuple{Vector{ListOfPoints{n_x}},Vector{Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}}}

end


function DiffFun(fun, x0::MSM{Point{n_x}}, epsilon=1e-6) where n_x

  xi = deepcopy(x0)

  N = length(x0.data)

  r0 = fun(x0)::MSM{Point{n_x}}

  zz = zeros( SMatrix{n_x, n_x, Float64, n_x*n_x}, N)

  JMat = MSM(zz, x0.sizes)::MSM{SMatrix{n_x,n_x,Float64,n_x*n_x}}
  
  Jarray = reshape(reinterpret(Float64, JMat.data), (n_x, n_x, N))

  for i_x=1:n_x

    add_epsilon!(xi.data, i_x, epsilon)

    fi = fun(xi)::MSM{Point{n_x}}
    di = (fi-r0)/epsilon

    Jarray[:,i_x,:] = reshape(reinterpret(Float64, di.data), (n_x, N))

    add_epsilon!(xi.data, i_x, -epsilon)

  end

  return r0, JMat

end

function DiffFun!(fun!, x0::MSM{Point{n_x}}, epsilon=1e-6, out=nothing) where n_x


  f0,fi,xi,JMat = out
  # xi = deepcopy(x0)

  xi.data .= x0.data

  N = length(x0.data)

  fun!(f0, x0)
  
  Jarray = reshape(reinterpret(Float64, JMat.data), (n_x, n_x, N))

  for i_x=1:n_x

    add_epsilon!(xi.data, i_x, epsilon)

    fun!(fi, xi)
    di = (fi-f0)/epsilon

    Jarray[:,i_x,:] = reshape(reinterpret(Float64, di.data), (n_x, N))

    add_epsilon!(xi.data, i_x, -epsilon)

  end

  return f0, JMat

end
#     for i_m=1:n_m
#       JMat[:,i_x,:] = reshape(reinterpret(Float64, vec(di[i_m])), (n_x, N))
#       add_epsilon!(xi[i_m], i_x, -epsilon)
#     end
#   end




struct NewtonResult
    solution
    iterations::Int64
    maxit::Int64
    err_ϵ::Float64
    tol_ϵ::Float64
    err_η::Float64
    tol_η::Float64
    errors::Vector{Float64}
end

function newton(fun::Function, x0::Vector{ListOfPoints{n_x}}, a::Union{Vector{ListOfPoints{n_x}},Nothing}=nothing, b::Union{Vector{ListOfPoints{n_x}},Nothing}=nothing; maxit=10, verbose=false, n_bsteps=5, lam_bsteps=0.5) where n_x

    steps = (lam_bsteps).^collect(0:n_bsteps)

    n_m = length(x0)
    N = length(x0[1])
    x = x0::Vector{ListOfPoints{n_x}}
    err_e_0 = 10.0
    tol_e = 1e-8
    err_η = 10.0
    tol_η = -1.0
    it=0
    errors = Float64[]

    while (it<maxit) && (err_e_0>tol_e)
        it += 1
        R_i, D_i = DiffFun(fun, x) # ::Tuple{Vector{ListOfPoints{n_x}},Vector{Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}}}
        if !(a isa Nothing)
            PhiPhi!(R_i,x,a,b,D_i)
        end
        err_e_0 = maxabs(R_i)
        push!(errors, err_e_0)
        err_e_1 = err_e_0
        dx = [[D_i[i][n]\R_i[i][n] for n=1:N] for i=1:n_m]::Vector{ListOfPoints{n_x}}
        err_η = maxabs(dx)
        i_bckstps = 0
        new_x = x
        while err_e_1>=err_e_0 && i_bckstps<length(steps)
            i_bckstps += 1
            new_x = x-dx*steps[i_bckstps]
            new_res = fun(new_x)::Vector{ListOfPoints{n_x}} # no diff
            if !(a isa Nothing)
                new_res = [PhiPhi0.(new_res[i],new_x[i],a[i],b[i]) for i=1:n_m]
            end
            err_e_1 = maxabs(new_res)
        end
        x = new_x::Vector{ListOfPoints{n_x}}
        if verbose
            println(i_bckstps, " : ", err_x, " : ", err_e_1, )
        end
    end

    res = NewtonResult(
        x,
        it,
        maxit,
        err_e_0,
        tol_e,
        err_η,
        tol_η,
        errors
    )
    return res
end

# TODO: cleanup derivative calculations

function newton(fun::Function, x0::MSM{Point{n_x}}; diff=true, maxit=10, verbose=false, n_bsteps=5, lam_bsteps=0.5) where n_x

  steps = (lam_bsteps).^collect(0:n_bsteps)

  n_m = length(x0)

  # N = length(x0[1])
  x = x0
  err_e_0 = 10.0
  tol_e = 1e-8
  err_η = 10.0
  tol_η = -1.0
  it=0
  errors = Float64[]

  while (it<maxit) && (err_e_0>tol_e)

      it += 1

      if diff
        R_i, D_i = DiffFun(fun, x) # ::Tuple{Vector{ListOfPoints{n_x}},Vector{Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}}}
      else
        R_i, D_i = fun(x)
      end

      err_e_0 = maxabs(R_i)
      push!(errors, err_e_0)
      err_e_1 = err_e_0
      dx = D_i \ R_i

      err_η = maxabs(dx)
      i_bckstps = 0
      new_x = x
      while err_e_1>=err_e_0 && i_bckstps<length(steps)
          i_bckstps += 1
          new_x = x-dx*steps[i_bckstps]
          if diff
          new_res = fun(new_x)::MSM{Point{n_x}} # no diff
          else
            new_res = fun(new_x)[1]
          end
          err_e_1 = maxabs(new_res)
      end
      x = new_x::MSM{Point{n_x}}
      if verbose
          println(i_bckstps, " : ", err_η, " : ", err_e_1, )
      end
  end

  res = NewtonResult(
      x,
      it,
      maxit,
      err_e_0,
      tol_e,
      err_η,
      tol_η,
      errors
  )
  return res
end
