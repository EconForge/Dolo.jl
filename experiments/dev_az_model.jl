using Dolo

model = yaml_import("examples/models/az_model.yaml")

dr_global = time_iteration(model)

