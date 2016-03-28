function solve_steady_state(m::DTCSCCModel, mc::ModelCalibration=m.calibration)
    p, s0, x0 = mc["parameters", "states", "controls"]
    e = zeros(mc["shocks"])
    ns = length(s)
    nx = length(x)

    function obj!(sx, out)
        s = sub(sx, 1:ns)
        X = sub(sx, ns+1:ns+nx)
        x_out = sub(out, ns+1:ns+nx)

        # update state part of residual
        S = evaluate(m.functions.transition, s, x, e, p)
        out[1:ns] = S-s

        # now update control part
        evaluate!(m.functions.arbitrage, s, x, e, S, X, p, x_out)
        out
    end

    sol = nlsolve(obj!, vcat(s0, x0))
    !converged(sol) &&  error("Nonlinear solver failed to find steady state")

    # otherwise set controls
    out = deepcopy(mc)
    out[m.symbolic.symbols[:states]...] = sol.zero[1:ns]
    out[m.symbolic.symbols[:controls]...] = sol.zero[ns+1:end]
    out
end
