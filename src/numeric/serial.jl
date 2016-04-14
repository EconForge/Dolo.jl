function serial_multiply(A::Array{Float64,3},B::Array{Float64,3})
    N,I,J = size(A)
    N,JJ,K = size(B)
    assert(J==JJ)

    C = zeros(N,I,K)
    for n = 1:N
        C[n,:,:] = slice(A,n,:,:)*slice(B,n,:,:)
    end
    return C
end

function serial_solver(fun::Function, x0::Array{Float64,2}, maxit::Int64)

    N = size(x0,1)
    n_x = size(x0,2)

    tol = 1e-8
    eps = 1e-8

    err = 1;
    it = 0;

    while (err>tol) && (it<maxit)
        res = fun(x0)
       N = size(res,1)

         # compute numerical gradient
        jac = zeros(N, n_x, n_x)
        for i = 1:n_x
            xx = x0
            xx[:,i] = xx[:,i] + eps
            dres = (fun(xx) - res)/eps
            jac[:,:,i] = dres
        end
        jac = permutedims(jac,[2,3,1])

        dx = zeros( size( x0 ) )
        for i = 1:N
            mat = jac[:,:,i]
            dx[i,:] = mat \ res[i,:]'
        end
        x = x0 - dx

        err = maximum(abs(res))
        it = it + 1
        x0 = x
    end

    return (x0,it)

end
