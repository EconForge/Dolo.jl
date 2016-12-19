"""
Computes the steady state of a model for a user-defined calibration.

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


"""
Computes the steady state of a model, where the calibration is provided by the model object.
"""
function steady_state_residuals(model)
    return steady_state_residuals(model, model.calibration)
end
