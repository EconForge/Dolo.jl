using Dolo

model_rbc_mc = yaml_import("examples/models/rbc_mc.yaml")
model_rbc = yaml_import("examples/models/rbc.yaml")
model_rbc_iid = yaml_import("examples/models/rbc_iid.yaml")

F_mc = Dolo.Euler(model_rbc_mc)
F = Dolo.Euler(model_rbc)
F_iid = Dolo.Euler(model_rbc_iid)

function multiply_tuple_elements(some_tuple)
    tuple_product = 1
    for k in 1:length(some_tuple)
        tuple_product *= some_tuple[k]
    end
    return tuple_product
end

# test the endo and exo grids of the model rbc
length(F.grid.exo) == Dolo.get_options(model_rbc)[:discretization][:exo][:n]
length(F.grid.endo) == multiply_tuple_elements(Dolo.get_options(model_rbc)[:discretization][:endo][:n])

# test the endo grid of the model rbc idd (no exogenous grid is specified in the options)
length(F_iid.grid.endo) == multiply_tuple_elements(Dolo.get_options(model_rbc_iid)[:discretization][:endo][:n])

# test the endo grid of the model rbc mc ("...")
length(F_mc.grid.exo) == Dolo.get_options(model_rbc_mc)[:discretization][:exo][:n]
length(F_mc.grid.endo) == multiply_tuple_elements(Dolo.get_options(model_rbc_mc)[:discretization][:endo][:n])
# Problem : apparently, !Cartesian is not recognized by Dolang.
# But, by checking the yaml file, we can still observe that the size of the grid in F is not correct.
