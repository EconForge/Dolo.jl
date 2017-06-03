function steady_state_residuals(model::AModel, calibration::ModelCalibration)
    m = calibration[:exogenous]
    s = calibration[:states]
    x = calibration[:controls]
    p = calibration[:parameters]
    res = Dolo.arbitrage(model, m, s, x, m, s, x, p)
    S = Dolo.transition(model, m, s, x, m, p)
    return Dict(:arbitrage=>res, :transition=>S-s)
end

function steady_state_residuals(model::AModel)
    return steady_state_residuals(model, model.calibration)
end

function find_deterministic_equilibrium(model::AModel, calibration::ModelCalibration)
    m, p, s0, x0 = calibration[:exogenous, :parameters, :states, :controls]
    ns = length(s0)
    nx = length(x0)

    function obj!(sx, out)
        s = view(sx, 1:ns)
        x = view(sx, ns+1:ns+nx)
        s_out = view(out, 1:ns)
        x_out = view(out, ns+1:ns+nx)

        # update state part of residual
        transition!(model, s_out, m, s, x, m, p)
        broadcast!(-, s_out, s_out, s)

        # now update control part
        arbitrage!(model, x_out, m, s, x, m, s, x, p)
        out
    end

    sol = nlsolve(obj!, vcat(s0, x0))
    NLsolve.converged(sol) || error("Nonlinear solver failed to find steady state")

    # otherwise set controls
    out = deepcopy(calibration)
    out[:states] = sol.zero[1:ns]
    out[:controls] = sol.zero[ns+1:end]
    out
end

function find_deterministic_equilibrium(model::AModel)
    return find_deterministic_equilibrium(model, model.calibration)
end

# ---------- #
# Docstrings #
# ---------- #

"""
    find_deterministic_equilibrium(model::AModel, [calib::ModelCalibration])

Solve for the steady state equilibrium of `model`, data in `cailb` to fill
in parameter values and provide an initial guess for the states and controls.
When no calibration is passed `model.calibration` is used

The `exogenous` variables at time t-1 (`m`) and t (`M`) are set to
`calib[:exogenous]`.

The deterministic equilibrium is found by solving for vectors `s` and `x`, such
that

1. `s = transition(m, s, x, m, p)`
2. `0 = arbitrage(m, s, x, m, s, x, p)`

"""
find_deterministic_equilibrium


"""
    steady_state_residuals(model::AModel, [calib::ModelCalibration])::Dict

Compute the steady state residuals for the aribtrage and transition equations
of `model`, when these functions are evaluated at the data in `calib`. If no
`calib` is provided, `model.calibration` will be used.

See the docstring for `find_deterministic_equilibrium` for more information
"""
steady_state_residuals
