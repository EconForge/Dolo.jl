function get_ss_derivatives(model)
    m,s,x,p = model.calibration[:exogenous,:states,:controls,:parameters]

    g_diff = eval_jacobians(model, transition!, length(s), (m, s, x, m), p)
    f_diff = eval_jacobians(model, arbitrage!, length(x), (m, s, x, m, s, x), p)

    return g_diff,f_diff
end

perturbate(p::IIDExogenous) = (zeros(0), zeros(0,0))
perturbate(p::VAR1) = (p.mu, p.R)

type PerturbationResult
    solution::BiTaylorExpansion
    generalized_eigenvalues::Vector
    stable::Bool     # biggest e.v. lam of solution is < 1
    determined::Bool # next eigenvalue is > lam + epsilon (MOD solution well defined)
    unique::Bool     # next eigenvalue is > 1
end

function Base.show(io::IO, pbr::PerturbationResult)
    @printf io "Perturbation Results\n"
    @printf io " * Decision Rule type: %s\n" string(typeof(PerturbationResult))
    @printf io "   * %s\n" show(pbr.solution)
    @printf io " * Blanchard-Kahn: %s\n" blanchard_kahn(pbr)
    @printf io "   * stable < %s\n" pbr.stable
    @printf io "   * determined < %s\n" pbr.determined
    @printf io "   * unique < %s\n" pbr.unique
end

blanchard_kahn(fos::PerturbationResult) = fos.stable && fos.unique

function perturbate_first_order(g_s, g_x, f_s, f_x, f_S, f_X)

    eigtol = 1.0+1e-6

    ns = size(g_s,1)
    nx = size(g_x,1)
    nv = ns + nx

    A = [eye(ns) zeros(ns, nx);
         -f_S     -f_X]

    B = [g_s g_x;
         f_s f_x]

    # do orderd QZ decomposition
    gs = schurfact(A, B)

    genvals = (abs(gs[:alpha]) ./ abs(gs[:beta]))
    sort!(genvals, rev=true)
    n_keep = ns # number of eigenvalues to keep
    diff = genvals[n_keep+1] - genvals[n_keep]
    eigtol = genvals[n_keep] + diff/2

    select = (abs(gs[:alpha]) .> eigtol*abs(gs[:beta]))

    ordschur!(gs, select)
    S, T, Q, Z = gs[:S], gs[:T], gs[:Q], gs[:Z]
    diag_S = diag(S)
    diag_T = diag(T)
    eigval = abs(diag_S./diag_T)

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

function get_gf_derivatives(model::AbstractNumericModel)

    g_diff,f_diff = get_ss_derivatives(model)
    _f_m,_f_s,_f_x,_f_M,_f_S,_f_X = f_diff
    _g_m,_g_s,_g_x,_g_M = g_diff

    (M, R) = perturbate(model.exogenous)

    f_x = _f_x
    f_X = _f_X
    _m,_s,x,p = model.calibration[:exogenous,:states,:controls,:parameters]

    if size(R,1)>0
        f_s = [_f_m _f_s]
        f_S = [_f_M _f_S]
        g_s = [R zeros(size(R,1),size(_g_s,2)); _g_m _g_s]
        g_x = [zeros(size(_g_m,1),size(_g_x,2)); _g_x]
        s = cat(1,_m, _s)
    else
        f_s = _f_s
        f_S = _f_S
        g_s = _g_s
        g_x = _g_x
        s = _s
    end
    return g_s, g_x, f_s, f_x, f_S, f_X
end


function perturbate(model::AbstractNumericModel; infos=false)

    g_s, g_x, f_s, f_x, f_S, f_X = get_gf_derivatives(model)
    nx = size(g_x, 2)

    (M, R) = perturbate(model.exogenous)
    _m,_s,x,p = model.calibration[:exogenous,:states,:controls,:parameters]

    if size(R,1)>0
        s = cat(1,_m, _s)
    else
        s = _s
    end

    C, genvals = perturbate_first_order(g_s, g_x, f_s, f_x, f_S, f_X)
    sort!(genvals)

    if size(R,1)>0
        C_exo = C[:,1:length(_m)]
        C_endo = C[:,(length(_m)+1):end]
        dr = BiTaylorExpansion(_m, _s, x, C_exo, C_endo)
    else
        dr = BiTaylorExpansion(_m, _s, x, zeros(nx, length(_m)), C)
    end

    tol = 1e-6 # minimum distance betweel lam_n and lam_{n+1}

    if !infos
        return dr
    else
        n_s = size(g_s,1)
        PerturbationResult(
            dr,
            genvals,
            genvals[n_s]<1,
            genvals[n_s+1]-genvals[n_s]>tol,
            genvals[n_s+1]>1
        )
    end



end
