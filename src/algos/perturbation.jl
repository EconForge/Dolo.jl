
# define custom exception types that hold some data for us to catch/explore
immutable GeneralizedEigenvaluesError <: Exception
    msg::AbstractString
    diag_S::Vector{Float64}
    diag_T::Vector{Float64}
end

immutable BlanchardKahnError <: Exception
    msg::AbstractString
    n_found::Int
    n_expected::Int
end

function _check_gev(diag_S, diag_T, tol)
    ok = sum((abs(diag_S) .< tol) .* (abs(diag_T) .< tol)) == 0
    if !ok
        msg = "Eigenvalues are not uniquely defined."
        throw(GeneralizedEigenvaluesError(msg, diag_S, diag_T))
    end
    true
end

function _check_bk_conditions(eigval, tol, n_expected)
    n_big = sum(eigval .> tol)
    if n_expected != n_big
        msg = """
        There are $(n_big) eigenvalues greater than one.
        There should be exactly $(n_expected) to meet Blanchard-Kahn conditions.
        """
        throw(BlanchardKahnError(msg, n_big, n_expected))
    end
    true
end

function get_ss_derivatives(model)
    m,s,x,p = model.calibration[:exogenous,:states,:controls,:parameters]

    g_diff = eval_jacobians(model, transition!, length(s), (m, s, x, m), p)
    f_diff = eval_jacobians(model, arbitrage!, length(x), (m, s, x, m, s, x), p)

    return g_diff,f_diff
end

perturbate(p::IIDExogenous) = (zeros(0), zeros(0,0))

perturbate(p::VAR1) = (p.M, p.R)

function perturbate(model::AbstractNumericModel, eigtol::Float64=1.0+1e-6)


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

    ns = length(s)
    nx = length(x)
    nv = ns + nx

    A = [eye(ns) zeros(ns, nx);
         -f_S     -f_X]
    B = [g_s g_x;
         f_s f_x]

    # do orderd QZ decomposition
    gs = schurfact(A, B)
    select = (abs(gs[:alpha]) .> eigtol*abs(gs[:beta]))
    ordschur!(gs, select)
    S, T, Q, Z = gs[:S], gs[:T], gs[:Q], gs[:Z]
    diag_S = diag(S)
    diag_T = diag(T)
    eigval = abs(diag_S./diag_T)

    # check eigen values and BK conditions
    _check_gev(diag_S, diag_T, 1e-10)
    _check_bk_conditions(eigval, eigtol, nx)

    # now obtain first order solution
    Z11 = Z[1:ns, 1:ns]
    Z12 = Z[1:ns, ns+1:end]
    Z21 = Z[ns+1:end, 1:ns]
    Z22 = Z[ns+1:end, ns+1:end]
    S11 = S[1:ns, 1:ns]
    T11 = T[1:ns, 1:ns]

    # first order solution
    C = (Z11'\Z21')'
    # P = (S11'\Z11')'*(Z11'\T11')'
    # Q = g_e
    # A = g_s + g_x*C
    # B = g_e
    C_exo = C[:,1:length(_m)]
    C_endo = C[:,(length(_m)+1):end]

    if size(R,1)>0
        BiTaylorExpansion(_m, _s, x, C_exo, C_endo)
    else
        BiTaylorExpansion(zeros(0), _s, x, C_exo, C_endo)
    end
end
