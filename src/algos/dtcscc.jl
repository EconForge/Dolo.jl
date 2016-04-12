function solve_steady_state(m::DTCSCCModel, mc::ModelCalibration=m.calibration)
    p, s0, x0 = mc[:parameters, :states, :controls]
    e_ = zeros(mc[:shocks])
    ns = length(s0)
    nx = length(x0)

    function obj!(sx, out)
        s = sub(sx, 1:ns)
        x = sub(sx, ns+1:ns+nx)
        s_out = sub(out, 1:ns)
        x_out = sub(out, ns+1:ns+nx)

        # update state part of residual
        evaluate!(m.functions.transition, s, x, e_, p, s_out)
        broadcast!(-, s_out, s_out, s)

        # now update control part
        evaluate!(m.functions.arbitrage, s, x, e_, s, x, p, x_out)
        out
    end

    sol = nlsolve(obj!, vcat(s0, x0))
    !converged(sol) &&  error("Nonlinear solver failed to find steady state")

    # otherwise set controls
    out = deepcopy(mc)
    out[:states] = sol.zero[1:ns]
    out[:controls] = sol.zero[ns+1:end]
    out
end

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

function linear_solve(m::DTCSCCModel, calib::ModelCalibration=m.calibration;
                      verbose::Bool=false, eigtol::Float64=1.0+1e-6)
    f = m.functions.arbitrage
    g = m.functions.transition

    p, s, x, e = calib[:parameters, :states, :controls, :shocks]
    ns = length(s)
    nx = length(x)
    ne = length(e)
    nv = ns + nx

    # evaluate Jacobians and stack into A,B matrix for qz
    g_s, g_x, g_e = eval_jacobians(g, (s, x, e), p)
    f_s, f_x, f_e, f_S, f_X = eval_jacobians(f, (s, x, e, s, x), p)

    A = [eye(ns) zeros(ns, nx)
         -f_S     -f_X]
    B = [g_s g_x
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

    dr = TaylorExpansion(s, x, C)
    # dr.A = A
    # dr.B = B
    # dr.sigma = get(m.distribution, :Normal, zeros(ne, ne))
end
