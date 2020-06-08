function numdiff(f, x0; ϵ=1e-6)
    f0 = f(x0)
    p = size(f0,1)
    q = size(x0,1)
    J = zeros(p,q)
    for j=1:q
        ei = zeros(q)
        ei[j] = 1.0
        J[:,j] = (f(x0+ϵ*ei) - f0)/ϵ
    end
    return J
end

function get_ss_derivatives(model)
    m, s, x, p = model.calibration[:exogenous, :states, :controls, :parameters]
    g_diff = transition(model,Val{(1,2,3,4)}, m, s, x, m, p)
    f_diff = arbitrage(model, Val{(1,2,3,4,5,6)}, m, s, x, m, s, x, p)
    # g_diff = [numdiff(f,u) for (f,u) in [
    #     (u->transition(model, u, s, x, m, p),m),
    #     (u->transition(model, m, u, x, m, p),s),
    #     (u->transition(model, m, s, u, m, p),x),
    #     (u->transition(model, u, s, x, u, p),m)
    # ]]

    # f_diff = [numdiff(f,u) for (f,u) in [
    #     (u->arbitrage(model, u, s, x, m, s, x, p), m),
    #     (u->arbitrage(model, s, u, x, m, s, x, p), s),
    #     (u->arbitrage(model, s, s, u, m, s, x, p), x),
    #     (u->arbitrage(model, s, s, x, u, s, x, p), m),
    #     (u->arbitrage(model, s, s, x, m, u, x, p), s),
    #     (u->arbitrage(model, s, s, x, m, s, u, p), x),
    # ]]
    g_diff, f_diff
end

perturb(p::IIDExogenous) = (zeros(0), zeros(0, 0))
perturb(p::VAR1) = (p.mu, p.R)

mutable struct PerturbationResult
    dr::BiTaylorExpansion
    generalized_eigenvalues::Union{Vector, Nothing}
    stable::Union{Bool, Nothing}     # biggest e.v. lam of solution is < 1
    determined::Union{Bool, Nothing} # next eigenvalue is > lam + epsilon (MOD solution well defined)
    unique::Union{Bool, Nothing}     # next eigenvalue is > 1
end

function Base.show(io::IO, pbr::PerturbationResult)
    @printf io "Perturbation Results\n"
    @printf io " * Decision Rule type: %s\n" string(typeof(pbr.dr))
    # @printf io " * Blanchard-Kahn: %s\n" blanchard_kahn(pbr)
    if !(typeof(pbr.stable)<:Nothing)
        @printf io " * stable < %s\n" pbr.stable
    end
    if !(typeof(pbr.determined)<:Nothing)
        @printf io " * determined < %s\n" pbr.determined
    end
end

blanchard_kahn(fos::PerturbationResult) = fos.stable && fos.unique

function perturb_first_order(g_s, g_x, f_s, f_x, f_S, f_X)

    eigtol = 1.0+1e-6

    ns = size(g_s, 2)
    nx = size(g_x, 2)

    nv = ns + nx

    A = [Matrix(1.0I,ns,ns) zeros(ns, nx);
         -f_S     -f_X]

    B = [g_s g_x;
         f_s f_x]

    # do orderd QZ decomposition
    gs = schur(A, B)

    genvals = (abs.(gs.α) ./ abs.(gs.β))
    sort!(genvals, rev=true)
    n_keep = ns # number of eigenvalues to keep
    diff = genvals[n_keep+1] - genvals[n_keep]
    eigtol = genvals[n_keep] + diff/2

    select = (abs.(gs.α) .> eigtol*abs.(gs.β))

    ordschur!(gs, select)
    S, T, Q, Z = gs.S, gs.T, gs.Q, gs.Z
    diag_S = diag(S)
    diag_T = diag(T)
    eigval = abs.(diag_S./diag_T)

    # now obtain first order solution
    Z11 = Z[1:ns, 1:ns]
    Z12 = Z[1:ns, ns+1:end]
    Z21 = Z[ns+1:end, 1:ns]
    Z22 = Z[ns+1:end, ns+1:end]
    S11 = S[1:ns, 1:ns]
    T11 = T[1:ns, 1:ns]

    # first order solution
    C = (Z11'\Z21')'
    return C, genvals
end

function get_gf_derivatives(model::AbstractModel)

    g_diff, f_diff = get_ss_derivatives(model)
    _f_m, _f_s, _f_x, _f_M, _f_S, _f_X = f_diff
    _g_m, _g_s, _g_x, _g_M = g_diff

    (M, R) = perturb(model.exogenous)

    f_x = _f_x
    f_X = _f_X
    _m, _s, x, p = model.calibration[:exogenous, :states, :controls, :parameters]

    if size(R, 1)>0
        f_s = [_f_m _f_s]
        f_S = [_f_M _f_S]
        g_s = [R zeros(size(R, 1), size(_g_s, 2)); _g_m _g_s]
        g_x = [zeros(size(_g_m, 1), size(_g_x, 2)); _g_x]
        s = cat(_m, _s; dims=1)
    else
        f_s = _f_s
        f_S = _f_S
        g_s = _g_s
        g_x = _g_x
        s = _s
    end
    return g_s, g_x, f_s, f_x, f_S, f_X
end

"""
TBD
"""
function perturb(model::Model; method=:qz)

    g_s, g_x, f_s, f_x, f_S, f_X = get_gf_derivatives(model)
    nx = size(g_x, 2)

    (M, R) = perturb(model.exogenous)
    _m, _s, x, p = model.calibration[:exogenous, :states, :controls, :parameters]

    if size(R, 1)>0
        s = cat(_m, _s; dims=1)

    else
        s = _s
    end

    n_s = size(g_s, 1)
    if method==:qz
        tol = 1e-6 # minimum distance betweel lam_n and lam_{n+1}
        C, genvals = perturb_first_order(g_s, g_x, f_s, f_x, f_S, f_X)
        sort!(genvals)
        stable = (genvals[n_s]<1)
        deter = (genvals[n_s+1]-genvals[n_s]>tol)
        uni = (genvals[n_s+1]>1)
    elseif method==:linear_time_iteration
        C, nit = linear_time_iteration(g_s, g_x, f_s, f_x, f_S, f_X)
        genvals = nothing
        stable = nothing
        deter = nit>0
        uni = nothing
    end
        

    if size(R, 1)>0
        C_exo = C[:, 1:length(_m)]
        C_endo = C[:, (length(_m)+1):end]
        dr = BiTaylorExpansion{nx}(_m, _s, x, C_exo, C_endo)
    else
        dr = BiTaylorExpansion{nx}(_m, _s, x, zeros(nx, length(_m)), C)
    end


    PerturbationResult(
        dr,
        genvals,
        stable,
        deter,
        uni
    )

end




function linear_time_iteration(model; maxit=1000, tol_ϵ=1e-8, improve=20, verbose=false)
    g_s, g_x, f_s, f_x, f_S, f_X = Dolo.get_gf_derivatives(model)
    return linear_time_iteration(g_s, g_x, f_s, f_x, f_S, f_X; maxit=maxit)
end


function L(P, Q, X, Y, u)
    return P*u*Q
end

function linear_time_iteration(g_s, g_x, f_s, f_x, f_S, f_X; maxit=1000, tol_ϵ=1e-8, improve=20, verbose=false)

    n_x, n_s = size(f_s)

    K = f_s + f_S*g_s

    X_0 = 1 .+rand(n_x, n_s)
    local X
    for it=1:maxit
        Y = X_0
        B = -K - f_X*Y*g_s
        A = f_x + f_S*g_x + f_X*Y*g_x
        X = A\B
        ϵ = maximum(abs.(X-X_0))
        if ϵ<tol_ϵ
            return (X,it)
        end
        if verbose
            println(it, "\t: ", ϵ)
        end
        if (improve>0)
            S = (X-X_0)
            δ = (X-X_0)
            P = -(f_x + f_S*g_x + f_X*Y*g_x)\f_X
            Q = g_s + g_x*X
            for k=1:improve
                S = δ + L(P, Q, X, Y, S)
            end
            X_0 = X_0 + S
        else
            X_0 = X
        end
    end

    return (X,-1)

end

function ttest(g_s, g_x, f_s, f_x, f_S, f_X, X)
    res = f_s + f_x*X + f_S*(g_s+g_x*X) + f_X*X*(g_s+g_x*X)
    return res
end
