import Dolo

model = Dolo.yaml_import("./model_tax.yaml")

model.calibration

m0, s0, x0, p = model.calibration[:exogenous, :states, :controls, :parameters]

Dolo.transition(model, m0, s0, x0, m0, p)



model = Dolo.yaml_import("examples/models/rbc.yaml")



model.calibration

m0, s0, x0, p = model.calibration[:exogenous, :states, :controls, :parameters]

Dolo.transition(model, m0, s0, x0, m0, p)


