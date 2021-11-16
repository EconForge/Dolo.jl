using Dolo


model = Model("examples/models/rbc.yaml")

Dolo.new_time_iteration(model);
Dolo.improved_time_iteration(model);


# grid, exo = Dolo.discretize(model)

println("Endo Grid size: $(Dolo.n_nodes(grid.endo))")
println("Exo Grid size: $(Dolo.n_nodes(grid.exo))")
println("Exo Grid I size: $(Dolo.n_inodes(exo, 1))")


grid, exo = Dolo.discretize(model;)

println("Endo Grid size: $(Dolo.n_nodes(grid.endo))")
println("Exo Grid size: $(Dolo.n_nodes(grid.exo))")
println("Exo Grid I size: $(Dolo.n_inodes(exo, 1))")

sol0 = time_iteration(model)



F = Dolo.Euler(model)

@time Dolo.loop_iti(F, F.x0)

@time Dolo.loop_ti(F, F.x0)


# μ = Dolo.ergodic_distribution(model, sol)



# using Plots