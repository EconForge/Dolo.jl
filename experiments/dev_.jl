using Dolo

model = Model("examples/models/rbc.yaml")


F = Dolo.Euler(model, ignore_constraints=false);

sol = improved_time_iteration(model, details=false)

# P, μ = Dolo.ergodic_distribution(model, sol.dr, F.grid_exo, F.grid_endo, F.dprocess);

x0 = Dolo.MSM([sol.dr(i,F.s0) for i in 1:length(sol.dr.grid_exo.nodes)])

@time Π, dΠ = Dolo.new_transition_dev(model, F.dprocess, x0, F.grid_exo, F.grid_endo);



using StaticArrays

m = @SMatrix [1.0 2.0 ; 3.0 4.0]
x = @SVector [0.4, 0.5]

1 ./ x .* m 


m