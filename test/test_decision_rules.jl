@testset "DR types" begin
    grid_smol = Dolo.SmolyakGrid{2}([0.0, 0.0], [1.0, 1.0,], [3,3])
    grid_exo = Dolo.EmptyGrid()
    grid_exo2 = Dolo.UnstructuredGrid{2}(rand(1, 2))
    grid_rand = Dolo.RandomGrid{2}([0.0, 0.0], [1.0, 1.0,], 50)
    grid_cart = Dolo.CartesianGrid{2}([0.0, 0.0], [1.0, 1.0], [10, 15])
    pts = rand(200, 2)
    pts_vec = reinterpret(Point{2}, pts', (size(pts, 1),))

    function test_dr_empty(dr, g_vals, g_vals_vec)
        set_values!(dr, [g_vals])
        out = evaluate(dr, pts)
        @test maximum(abs, out - evaluate(dr, pts_vec)) < 1e-14
        set_values!(dr, [g_vals_vec])
        @test maximum(abs, out - evaluate(dr, pts_vec)) < 1e-14
        @test maximum(abs, out - evaluate(dr, pts)) < 1e-14
        @test_throws MethodError evaluate(dr, pts_vec[1])
        @test_throws ErrorException set_values!(dr, [g_vals_vec, g_vals_vec])
    end

    function test_dr_unstructured(dr, g_vals, g_vals_vec)
        set_values!(dr, [g_vals])
        out = evaluate(dr, 1, pts)
        @test maximum(abs, out - evaluate(dr, 1, pts_vec)) < 1e-14
        set_values!(dr, [g_vals_vec])
        @test maximum(abs, out - evaluate(dr, 1, pts_vec)) < 1e-14
        @test maximum(abs, out - evaluate(dr, 1, pts)) < 1e-14
        @test_throws MethodError evaluate(dr, 1, pts_vec[1])
        @test_throws ErrorException set_values!(dr, [g_vals_vec, g_vals_vec])
    end

    test_sets = [
        (grid_smol, Dolo.SmolyakDR, Dolo.Smolyak),
        (grid_rand, Dolo.CompletePolyDR, Dolo.CompletePolnomial{3}),
        (grid_rand, Dolo.CompletePolyDR, Dolo.CompletePolnomial),
        (grid_rand, Dolo.CompletePolyDR, Dolo.CompletePolnomial{2}),
        (grid_smol, Dolo.CompletePolyDR, Dolo.CompletePolnomial{3}),
        (grid_smol, Dolo.CompletePolyDR, Dolo.CompletePolnomial),
        (grid_smol, Dolo.CompletePolyDR, Dolo.CompletePolnomial{2}),
        (grid_exo2, Dolo.CompletePolyDR, Dolo.CompletePolnomial{3}),
        (grid_exo2, Dolo.CompletePolyDR, Dolo.CompletePolnomial),
        (grid_exo2, Dolo.CompletePolyDR, Dolo.CompletePolnomial{2}),
        (grid_cart, Dolo.CompletePolyDR, Dolo.CompletePolnomial{3}),
        (grid_cart, Dolo.CompletePolyDR, Dolo.CompletePolnomial),
        (grid_cart, Dolo.CompletePolyDR, Dolo.CompletePolnomial{2}),
        (grid_cart, Dolo.BasisMatricesDR, Dolo.Chebyshev),
    ]

    @testset for (grid_endo, dr_name, drT) in test_sets

        dr1 = dr_name(grid_exo, grid_endo, Val{2})
        g = nodes(dr1)
        g_vals = [sin.(g[:, 1]) cos.(g[:, 2])]
        g_vals_vec = reinterpret(Point{2}, g_vals', (size(g,1),))
        test_dr_empty(dr1, g_vals, g_vals_vec)

        dr2 = Dolo.DecisionRule(grid_exo, grid_endo, Val{2}, drT)
        test_dr_empty(dr2, g_vals, g_vals_vec)

        dru = dr_name(grid_exo2, grid_endo, Val{2})
        test_dr_unstructured(dru, g_vals, g_vals_vec)

    end

    # test BSpline
    function test_bspline_dr(order)
        dr1 = Dolo.DecisionRule(grid_exo, grid_cart, Val{2}, Dolo.BSpline{order})
        g = nodes(dr1)
        g_vals = [sin.(g[:, 1]) cos.(g[:, 2])]
        g_vals_vec = reinterpret(Point{2}, g_vals', (size(g,1),))
        test_dr_empty(dr1, g_vals, g_vals_vec)

        dr2 = Dolo.DecisionRule(grid_exo2, grid_cart, Val{2}, Dolo.BSpline{order})
        test_dr_unstructured(dr2, g_vals, g_vals_vec)
    end
    @testset for order in 1:4
        test_bspline_dr(order)
    end

    @testset "PWLinear" begin
        dr1 = Dolo.DecisionRule(grid_exo, grid_cart, Val{2}, Dolo.PWLinear)
        g = nodes(dr1)
        g_vals = [sin.(g[:, 1]) cos.(g[:, 2])]
        g_vals_vec = reinterpret(Point{2}, g_vals', (size(g,1),))
        test_dr_empty(dr1, g_vals, g_vals_vec)

        dr2 = Dolo.DecisionRule(grid_exo2, grid_cart, Val{2}, Dolo.PWLinear)
        test_dr_unstructured(dr2, g_vals, g_vals_vec)
    end

end
