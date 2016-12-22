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
    res = Dolo.arbitrage(model,m,s,x,m,s,x,p)
    S = Dolo.transition(model,m,s,x,m,p)
    return Dict(:arbitrage=>res,:transition=>S-s)
end


function steady_state_residuals(model)
    return steady_state_residuals(model, model.calibration)
end
