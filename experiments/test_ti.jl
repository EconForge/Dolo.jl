using Dolo


model = Model("examples/models/rbc_mc.yaml")


Dolo.get_calibration(model; beta=0.95)

Dolo.set_calibration!(model; beta=0.94)


model.calibration.flat[:beta]



sol0 = time_iteration(model)
sol = time_iteration(model, sol0.dr)




ssol0 = improved_time_iteration(model)

ssol = improved_time_iteration(model, ssol0.dr)