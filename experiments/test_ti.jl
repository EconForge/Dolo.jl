using Dolo


model = Model("experiments/cons.yaml")



sol0 = time_iteration(model)
sol = time_iteration(model, sol0.dr)

improved_time_iteration(model, sol0.dr)



μ = Dolo.ergodic_distribution(model, sol)



using Plots