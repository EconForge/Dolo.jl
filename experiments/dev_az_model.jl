using Dolo

model = yaml_import("examples/models/az_model.yaml")

gr = Dolo.discretize(model)

sol = improved_time_iteration(model)

# sim = simulate(model, sol.dr)

tab = tabulate(model, sol.dr, :k)

using SimplePlots

plot(tab[:k], tab[:i])