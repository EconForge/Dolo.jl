# function smooth(x::Array{Float64,1},a::Array{Float64,1},b::Array{Float64,1},f::Array{Float64,1},J::Array{Float64,2})
#
#     BIG = 1e20
#
#     da = a - x
#     db = b - x
#
#
#     if a<=-BIG
#         pval = fx
#         dpdy = 1.0 + da*0
#         dpdz = 0.0 + da*0
#     else
#         sq1 = sqrt.( fx.^2 + da.^2)
#         pval = fx + sq1 + da
#         dpdy = 1.0 + fx./sq1
#         dpdz = 1.0 + da./sq1
#     end
#     if b>=BIG
#         fxnew = pval
#         dmdy = 1.0
#         dmdz = 0
#     else
#         sq2 = sqrt(pval.^2 + db.^2)
#         fxnew = pval - sq2 + db
#         dmdy = 1.0 - pval./sq2
#         dmdz = 1.0 - db./sq2
#     end
#     ff = dmdy.*dpdy
#     xx = dmdy.*dpdz + dmdz
#
# end


function mcp_smooth(x::Vector, f::Vector, J::Matrix, lower::Vector, upper::Vector)

    fx = copy(f)

    for i = 1:length(x)
        if isfinite(upper[i])
            fx[i] += (x[i]-upper[i]) + sqrt(fx[i]^2+(x[i]-upper[i])^2)
        end
        if isfinite(lower[i])
            fx[i] += (x[i]-lower[i]) - sqrt(fx[i]^2+(x[i]-lower[i])^2)
        end
    end

    return fx

end

function mcp_diff(x::Vector, f::Vector, J::Matrix, lower::Vector, upper::Vector)

    fx = copy(f)
    gx = copy(J)

    for i = 1:length(x)
        if isfinite(upper[i])
            fx[i] += (x[i]-upper[i]) + sqrt(fx[i]^2+(x[i]-upper[i])^2)
        end
        if isfinite(lower[i])
            fx[i] += (x[i]-lower[i]) - sqrt(fx[i]^2+(x[i]-lower[i])^2)
        end
    end

    # Derivatives of phiplus
    sqplus = sqrt(fx.^2+(x-upper).^2)

    dplus_du = 1 + fx./sqplus

    dplus_dv = similar(x)
    for i = 1:length(x)
        if isfinite(upper[i])
            dplus_dv[i] = 1 + (x[i]-upper[i])/sqplus[i]
        else
            dplus_dv[i] = 0
        end
    end

    # Derivatives of phiminus
    phiplus = copy(fx)
    for i = 1:length(x)
        if isfinite(upper[i])
            phiplus[i] += (x[i]-upper[i]) + sqplus[i]
        end
    end

    sqminus = sqrt(phiplus.^2+(x-lower).^2)

    dminus_du = 1-phiplus./sqminus

    dminus_dv = similar(x)
    for i = 1:length(x)
        if isfinite(lower[i])
            dminus_dv[i] = 1 - (x[i]-lower[i])/sqminus[i]
        else
            dminus_dv[i] = 0
        end
    end

    # Final computations
    for i = 1:length(x)
        for j = 1:length(x)
            gx[i,j] *= dminus_du[i]*dplus_du[i]
        end
        gx[i,i] += dminus_dv[i] + dminus_du[i]*dplus_dv[i]
    end

    return fx,gx

end

function serial_solver(fun::Function, x0::Array{Float64,2}; maxit=10, verbose=true, a=zeros(0,0), b=zeros(0,0))


    N = size(x0,1)
    n_x = size(x0,2)

    if size(a) != (N,n_x)
        a = -ones(N,n_x)*Inf
    end
    if size(b) != (N,n_x)
        b = ones(N,n_x)*Inf
    end

    tol = 1e-8
    eps = 1e-8

    err = 1;
    it = 0;

    n_bsteps = 5
    backsteps = 0.5.^(0:(n_bsteps-1))

    x = x0
    res = fun(x0)
    N = size(res,1)
    err = maximum(abs(res))
    err_0 = err
    if verbose
        println("Initial error: ", err_0)
    end


    while (err>tol) && (it<maxit)

         # compute numerical gradient
        jac = zeros(N, n_x, n_x)
        for i = 1:n_x
            xx = copy(x0)
            xx[:,i] +=  eps
            jac[:,:,i] = (fun(xx) - res)/eps
        end

        for n=1:N
            r1, r2 = mcp_diff(x0[n,:], res[n,:], jac[n,:,:], a[n,:], b[n,:])
            res[n,:] = r1
            jac[n,:,:] = r2
        end

        dx = zeros( size( x0 ) )
        for n = 1:size(x0,1)
            mat = jac[n,:,:]
            dx[n,:] = mat \ res[n,:]
        end

        i = 0
        for i=1:n_bsteps
            lam = backsteps[i]
            x = x0 - lam*dx
            try
                res = fun(x)
                for n=1:N
                    r1 = mcp_smooth(x0[n,:], res[n,:], jac[n,:,:], a[n,:], b[n,:])
                    res[n,:] = r1
                end
                err = maximum(abs(res))
            catch
                err = Inf
            end
            if err<err_0
                break
            end
        end
        it = it + 1

        if verbose
            println("It: ", it, " ; Err: ", err)
        end

        err_0 = err
        x0 = x

    end

    return (x0,it)

end
