"""
Computes the residuals of the arbitrage and transition equations at the steady state of the model.

If the calibration for the model is not explicitly provided, the calibration is that associated with the model object.

# Arguments
* `model::NumericModel`: Model object that describes the current model environment.
* `calibration::OrderedDict`: Contains the model calibration.
# Returns
* `residuals::Dict`: Contains the residuals of the arbitrage equations, and the residuals of the transition equations.
"""
function steady_state_residuals(model, calibration)
    m = calibration[:exogenous]
    s = calibration[:states]
    x = calibration[:controls]
    p = calibration[:parameters]
    res = Dolo.arbitrage(model, m, s, x, m, s, x, p)
    S = Dolo.transition(model, m, s, x, m, p)
    return Dict(:arbitrage=>res, :transition=>S-s)
end


function steady_state_residuals(model)
    return steady_state_residuals(model, model.calibration)
end

 

function find_deterministic_equilibrium(model, calibration)
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
    !converged(sol) &&  error("Nonlinear solver failed to find steady state")

    # otherwise set controls
    out = deepcopy(calibration)
    out[:states] = sol.zero[1:ns]
    out[:controls] = sol.zero[ns+1:end]
    out
end

function find_deterministic_equilibrium(model)
    return find_deterministic_equilibrium(model, model.calibration)
end