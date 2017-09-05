include("tmp_module.jl")

include("test_processes.jl")
# this loads:
# - mc: markov chain
# - mvn: multivariate normal law


grid = CartesianGrid([0.0,1.0], [0.1, 0.4], [5,4])
gg = nodes(grid)
n_x = 3
values = rand(size(gg,1), n_x)
vv = [values for i=1:size(mc.values,1)]
dr = DecisionRule(mc, grid, vv)
set_values!(dr,vv) # filter coefficients

dr(1, [0.1 0.5])
s0 = [0.2, 0.5]'
res =  dr(1, s0)


grid = CartesianGrid([0.0,1.0], [0.1, 0.4], [5,4])
gg = nodes(grid)
n_x = 3
values = rand(size(gg,1), n_x)
vv = [values]
dr = DecisionRule(mvn, grid, vv)
set_values!(dr, values) # filter coefficients
s0 = [0.2, 0.5]'
res =  dr(1, s0)
dr(1, [0.1 0.5])

@testset "DR types" begin
    grid_smol = Dolo.SmolyakGrid{2}([0.0, 0.0], [1.0, 1.0,], [3,3])
    grid_exo = Dolo.EmptyGrid()
    grid_exo2 = Dolo.UnstructuredGrid{2}(rand(1, 2))
    grid_rand = Dolo.RandomGrid{2}([0.0, 0.0], [1.0, 1.0,], 50)
    pts = rand(10000, 2)
    pts_vec = reinterpret(Point{2}, pts', (size(pts, 1),))

    @testset for (grid_endo, drT) in [(grid_smol, Dolo.SmolyakDR),
                                      (grid_rand, Dolo.CompletePolyDR)]


        dr = drT(grid_exo, grid_endo, Val{2})
        g = nodes(grid_endo)
        sg_vals = [sin.(g[:, 1]) cos.(g[:, 2])]
        sg_vals_vec = reinterpret(Point{2}, sg_vals', (size(g,1),))
        set_values!(dr, [sg_vals])
        out = evaluate(dr, pts)
        @test maximum(abs, out - evaluate(dr, pts_vec)) < 1e-14
        set_values!(dr, [sg_vals_vec])
        @test maximum(abs, out - evaluate(dr, pts_vec)) < 1e-14
        @test maximum(abs, out - evaluate(dr, pts)) < 1e-14
        @test_throws MethodError evaluate(dr, pts_vec[1])
        @test_throws ErrorException set_values!(dr, [sg_vals_vec, sg_vals_vec])

        dr2 = drT(grid_exo2, grid_endo, Val{2})
        set_values!(dr2, [sg_vals])
        out = evaluate(dr2, 1, pts)
        @test maximum(abs, out - evaluate(dr2, 1, pts_vec)) < 1e-14
        set_values!(dr2, [sg_vals_vec])
        @test maximum(abs, out - evaluate(dr2, 1, pts_vec)) < 1e-14
        @test maximum(abs, out - evaluate(dr2, 1, pts)) < 1e-14
        @test_throws MethodError evaluate(dr2, pts)
        @test_throws MethodError evaluate(dr2, pts_vec)
        @test_throws MethodError evaluate(dr2, 1, pts_vec[1])
        @test_throws BoundsError evaluate(dr2, 2, pts_vec)
        @test_throws ErrorException set_values!(dr2, [sg_vals_vec, sg_vals_vec])
    end
end
