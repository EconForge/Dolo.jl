using Dolo

# model = yaml_import("examples/models/rbc.yaml")


model = yaml_import("experiments/ayiagari.yaml")



Dolo.get_domain(model.exogenous)

gr = Dolo.discretize(model)


sol = improved_time_iteration(model)

# sim = simulate(model, sol.dr)

tab = tabulate(model, sol.dr, :a)


using SimplePlots


plot(tab[:a], tab[:c])