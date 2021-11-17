using Dolo


model = Model("examples/models/rbc.yaml")


model = Model("consumption_savings.yaml")

# F = Dolo.Euler(model);

# F(F.x0, F.x0)

# fun(u) = F(u, F.x0, false)

# fun(F.x0)

# R_i, D_i = Dolo.DiffFun(fun, F.x0)

# r = F(F.x0, F.x0, false)
# j = Dolo.df_A(F, F.x0, F.x0)

# Dolo.newton(u->F(u,F.x0), F.x0)



sol = Dolo.time_iteration(model, maxit=1000);
tab = tabulate(model, sol.dr, :m)
pl = plot(tab[:m], tab[:m])
plot!(pl, tab[:m], tab[:c])



Dolo.time_iteration(model; ignore_constraints=true);


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
