using Dolo

model = yaml_import("experiments/ayiagari.yaml")

Dolo.get_domain(model.exogenous)

Dolo.discretize(model)


sol = improved_time_iteration(model)

tab = tabulate(model, sol.dr, :a)
