function solve_steady_state(m::DTCSCCModel, mc::ModelCalibration=m.calibration)
    p, s0, x0 = mc["parameters", "states", "controls"]
    e_ = zeros(mc["shocks"])
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
    out["states"] = sol.zero[1:ns]
    out["controls"] = sol.zero[ns+1:end]
    out
end
