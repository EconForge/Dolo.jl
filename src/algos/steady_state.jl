function calibrated_steady_state(model::AModel)

    c = model.calibration
    s = NamedTuple(k=>c[k] for k in variables(model.states))
    x = NamedTuple(k=>c[k] for k in variables(model.controls))
    return (;s,x)

end
