using Dolo


model = Model("examples/models/rbc.yaml")
println(model.exogenous)

Dolo.set_calibration!(model; sig_z=0.01)
println(model.exogenous)

Dolo.set_calibration!(model; sig_z=0.01)
println(model.exogenous)