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
