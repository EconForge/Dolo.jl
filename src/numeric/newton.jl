function add_epsilon!(x::ListOfPoints{d}, i, epsilon) where d
  ei = SVector{d,Float64}([(j==i?epsilon:0.0) for j=1:d])
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
        JMat[i_m][:,i_x,:] = reinterpret(Float64, di[i_m], (n_x, N))
        add_epsilon!(xi[i_m], i_x, -epsilon)
      end
    end
    J = [reinterpret(SMatrix{n_x,n_x,Float64,n_x^2},JMat[i],(N,)) for i=1:n_m]
    return (r0,J)::Tuple{Vector{ListOfPoints{n_x}},Vector{Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}}}
end


function newton(fun::Function, x0::Vector{ListOfPoints{n_x}}, a::Vector{ListOfPoints{n_x}}, b::Vector{ListOfPoints{n_x}}; maxit=10, verbose=false, n_bsteps=5, lam_bsteps=0.5) where n_x

    steps = (lam_bsteps).^collect(0:n_bsteps)

    n_m = length(x0)
    N = length(x0[1])
    x = x0::Vector{ListOfPoints{n_x}}
    err_e_0 = 10
    tol_e = 1e-8
    it=0
    while (it<maxit) && (err_e_0>tol_e)
        it += 1
        R_i, D_i = DiffFun(fun, x)::Tuple{Vector{ListOfPoints{n_x}},Vector{Vector{SMatrix{n_x,n_x,Float64,n_x*n_x}}}}
        PhiPhi!(R_i,x,a,b,D_i)
        err_e_0 = maxabs(R_i)
        err_e_1 = err_e_0
        dx = [[D_i[i][n]\R_i[i][n] for n=1:N] for i=1:n_m]::Vector{ListOfPoints{n_x}}
        err_x = maxabs(dx)
        i_bckstps = 0
        new_x = x
        while err_e_1>=err_e_0 && i_bckstps<length(steps)
            i_bckstps += 1
            new_x = x-dx*steps[i_bckstps]
            new_res = fun(new_x)::Vector{ListOfPoints{n_x}} # no diff
            new_res = [PhiPhi0.(new_res[i],new_x[i],a[i],b[i]) for i=1:n_m]
            err_e_1 = maxabs(new_res)
        end
        x = new_x::Vector{ListOfPoints{n_x}}
        if verbose
            println(i_bckstps, " : ", err_x, " : ", err_e_1, )
        end
    end

    return (x, it)::Tuple{Vector{ListOfPoints{n_x}}, Int}

end
